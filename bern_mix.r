#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
# 
# This file is part of hct.
# 
# hct is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hct is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hct.  If not, see <http://www.gnu.org/licenses/>.

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Fit bernoulli mix model on site/patient data"))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Table of site/patient direction of change",
                                 "data"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--q_thres",
                     help = paste("Threshold of differece for q_i - Q"),
                     type = "numeric",
                     default = 0.1)
  p <- add_argument(p, "--min_patients",
                    help = "Minimum number of patients",
                    type = "numeric",
                    default = 5)
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--iter",
                    help = "Number of Stan iterations per chain.",
                    type = "numeric",
                    default = 3000)
  p <- add_argument(p, "--warmup",
                    help = "Number of Stan warmup iterations per chain.",
                    type = "numeric",
                    default = 2000)
  p <- add_argument(p, "--chains",
                    help = paste("Number of Stan chains. It will be equal",
                                 "to the number of required threads.")
                    type = "numeric",
                    default = 4)
  p <- add_argument(p, "--vp",
                    help = paste("Variance parameter for distribution of p_i"),
                    type = "numeric",
                    default = 5)
  p <- add_argument(p, "--vq",
                    help = paste("Variance parameter for distribution of q_i"),
                    type = "numeric",
                    default = 5)
                    
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$iter < 0){
    stop("ERROR: iter must be positive", call. = TRUE)
  }
  if(args$warmup < 0){
    stop("ERROR: warmup must be positive", call. = TRUE)
  }
  if(args$warmup >= args$iter){
    stop("ERROR: warmup must be less than iter", call. = TRUE)
  }
  if(args$chains < 1){
    stop("ERROR: chains must be at least 1", call. = TRUE)
  }
  if(args$vp <= 0){
    stop("ERROR: vp must be positive", call. = TRUE)
  }
  if(args$vq <= 0){
    stop("ERROR: vq must be positive", call. = TRUE)
  }
  
  return(args)
}

args <- process_arguments()
# args <- list(input = "sites_dist.tsv.gz",
#              q_thres = 0.1,
#              min_patients = 5,
#              outdir = "output",
#              iter = 3000,
#              warmup = 2000,
#              chains = 4,
#              vp = 10,
#              vq = 10)

library(tidyverse)
library(rstan)

dat <- read_tsv(args$input,
                col_types = cols(site_id = col_character()))
if(max(dat$n_patients) < args$min_patients){
  cat("Not enough patients.\n")
  q()
}

# Creating site index, no need to randomly subset anymore
# # subseting randomly
# nsites <- 1000
# set.seed(38924)
# dat <- dat %>%
#   filter(site_id %in% sample(unique(dat$site_id), size = nsites))
# sites <- sample(sites, size = nsites, replace = FALSE)
sites <- tibble(site_id = dat$site_id,
                id = 1:length(dat$site_id))

# Prepare data for stan
cat("Preparing data for stan...\n")
stan_data <- dat %>%
  filter(site_id %in% sites$site_id) %>%
  select(site_id, n_decrease, n_equal, n_increase) %>%
  pmap_dfr(function(site_id, n_decrease, n_equal, n_increase){
    tibble(site_id = site_id,
           x = rep(c(-1,0,1), times = c(n_decrease, n_equal, n_increase)))
  }) %>%
  left_join(sites, by = "site_id")
stan_data <- list(x = stan_data$x,
     id = stan_data$id,
     nobs = length(stan_data$x),
     nsites = max(stan_data$id),
     vp = args$vp,
     vq = args$vq)

cat("Compiling stan model...\n")
# First find the location of the file via introspection
stan_file <- commandArgs(trailingOnly = FALSE)
i <- which(stan_file %>%
             str_detect("--file"))
stan_file <- stan_file[i] %>%
  str_remove("^--file=") %>%
  dirname()
stan_file <- file.path(stan_file, "stan", "bernoulli_mix_multisite.stan")
# stan_file <- "~/micropopgen/src/hct/stan/bernoulli_mix_multisite.stan"
# stan_file <- "~/micropopgen/src/hct/stan/bernoulli_mix_multisite_hyper.stan"
m1.model <- stan_model(stan_file,
                       model_name = "bern_change")

cat("Running Stan...\n")
# 100 sites 8 seconds
# 1000 sites 192 seconds, 241 seconds for hyper, 500 seconds for double hyper
# 10000 sites 4832 seconds (~1.3 hrs)
# 100000 if trend continues linearlu, ~53 hrs or  ~2.5 days
m1.stan <- sampling(m1.model,
                    data = stan_data,
                    chains = args$chains,
                    iter = args$iter,
                    warmup = args$warmup,
                    thin = 1,
                    cores = args$chains)
# load("output/m1.stan.rdat")
load("~/micropopgen/exp/2021/today/output/m1.stan.rdat")
pars <- c("P", "Q")
print(m1.stan, pars = pars)
# bayesplot::mcmc_pairs(m1.stan, pars = pars)
# bayesplot::mcmc_trace(m1.stan, pars = pars)
# bayesplot::mcmc_acf(m1.stan, pars = pars)

# mu <- 0.2
# v <- 10
# hist(rbeta(n = 1000, shape1 = mu * v,  shape2 = (1-mu) * v), breaks = 20)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

filename <- file.path(args$outdir, "model_summaries.tsv.gz")
write_tsv(summary(m1.stan)$summary %>%
            as.data.frame() %>%
            rownames_to_column(var = "var"),
          filename)

# Save model
cat("Saving stan model...\n")
filename <- file.path(args$outdir, "m1.stan.rdat")
save(m1.stan, file = filename)

cat("Extracting posterior...\n")
# Do I need post? or can iterate over stan object?
post <- rstan::extract(m1.stan)
rm(m1.stan)
gc()
# # First find the prob that the probability of change (p) greater than average (P)
# P_p_diff <- post$p - as.vector(post$P)
# # dim(P_p_diff)
# P_p_diff <- colSums(P_p_diff > 0) / nrow(P_p_diff)
# 
# # Second find that the probability of increase given change (q) greater than average (Q)
# Q_q_diff <- post$q - as.vector(post$Q)
# # dim(Q_q_diff)
# Q_q_diff <- colSums(Q_q_diff > 0) / nrow(Q_q_diff)
# 
# ids <- which(P_p_diff > 0.8 & Q_q_diff > 0.8)
# site_ids <- sites %>%
#   filter(id %in% ids) %>%
#   select(site_id) %>% unlist %>% as.character
# dat %>%
#   filter(site_id %in% c(site_ids)) %>%
#   print(n = 100)

cat("Calculating p_directional...\n")
res <- rep(0, stan_data$nsites)
for(iter in 1:((args$iter - args$warmup) * args$chains)){
  # i <- 1
  # post <- matrix(1:20, nrow = 5)
  # post
  # post - 1:5
  
  res <- res + 1*( (post$p[iter,] - post$P[iter]) > 0 & (abs(post$q[iter,] - post$Q[iter]) > args$q_thres) )
}

cat("Joining results and wirting output...\n")
res <- tibble(id = 1:length(res),
       p_directional = res / length(post$P))
res <- dat %>%
  inner_join(sites %>% left_join(res, by = "id"),
            by = "site_id")
filename <- file.path(args$outdir, "p_directional.tsv.gz")
write_tsv(res, filename)
 
# res %>%
#   arrange(desc(p_directional))

