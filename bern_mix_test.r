library(tidyverse)
library(rstan)

args <- list(post_thres = 0.8)

dat <- read_tsv("sites_dist.tsv.gz",
                col_types = cols(site_id = col_character()))


sites <- dat$site_id

nsites <- 10000
set.seed(38924)
sites <- sample(sites, size = nsites, replace = FALSE)
sites <- tibble(site_id = sites,
                id = 1:length(sites))

stan_data <- dat %>%
  filter(site_id %in% sites$site_id) %>%
  select(site_id, n_decrease, n_equal, n_increase) %>%
  pmap_dfr(function(site_id, n_decrease, n_equal, n_increase){
    tibble(site_id = site_id,
           x = rep(c(-1,0,1), times = c(n_decrease, n_equal, n_increase)))
  }) %>%
  left_join(sites, by = "site_id")
stan_data  
stan_data <- list(x = stan_data$x,
     id = stan_data$id,
     nobs = length(stan_data$x),
     nsites = max(stan_data$id))


m1.model <- stan_model("bernoulli_mix_multisite.stan",
                       model_name = "bern_change")
# 100 sites 8 seconds
# 1000 sites 192 seconds
# 10000 sites 4832 seconds (~1.3 hrs)
# 100000 if continueas linearlu, ~53 hrs or  ~2.5 days
m1.stan <- sampling(m1.model,
                    data = stan_data,
                    chains = 4,
                    iter = 2000,
                    warmup = 1000,
                    thin = 1,
                    cores = 4)
print(m1.stan)

post <- rstan::extract(m1.stan)
# First find the prob that the probability of change (p) greater than average (P)
P_p_diff <- post$p - as.vector(post$P)
# dim(P_p_diff)
P_p_diff <- colSums(P_p_diff > 0) / nrow(P_p_diff)

# Second find that the probability of increase given change (q) greater than average (Q)
Q_q_diff <- post$q - as.vector(post$Q)
# dim(Q_q_diff)
Q_q_diff <- colSums(Q_q_diff > 0) / nrow(Q_q_diff)

ids <- which(P_p_diff > 0.8 & Q_q_diff > 0.8)
site_ids <- sites %>%
  filter(id %in% ids) %>%
  select(site_id) %>% unlist %>% as.character
dat %>%
  filter(site_id %in% c(site_ids)) %>%
  print(n = 100)


res <- rep(0, stan_data$nsites)
for(iter in 1:length(post$P)){
  # i <- 1
  # post <- matrix(1:20, nrow = 5)
  # post
  # post - 1:5
  
  res <- res + 1*( (post$p[iter,] - post$P[iter]) > 0 & (abs(post$q[iter,] - post$Q[iter]) > 0.2) )
}
res <- tibble(id = 1:length(res),
       p_directional = res / length(post$P)) %>%
  arrange(desc(p_directional)) %>%
  print(n = 100)

res <- dat %>%
  inner_join(sites %>% left_join(res, by = "id"),
            by = "site_id")
write_tsv(res, "p_directional.tsv.gz")

res %>%
  arrange(desc(p_directional))
res %>%
  filter(pval_inc < 0.05 & pval_dec < 0.05)


sites %>%
  filter(site_id == "829901")
id_num <- 388

pars <- names(m1.stan)


selected_pars <- pars[c(1,2, 2 + id_num, 2 + stan_data$nsites + id_num)]
selected_pars

bayesplot::mcmc_areas(m1.stan, selected_pars)
bayesplot::mcmc_pairs(m1.stan, selected_pars)
bayesplot::mcmc_intervals(m1.stan, selected_pars)
bayesplot::mcmc_acf(m1.stan, selected_pars)
