#!/usr/bin/env Rscript

args <- list(pdir = "p_directional/MGYG-HGUT-00074.tsv.gz",
             s_coef = "s_coef/MGYG-HGUT-00074.tsv",
             fit = "FIT/MGYG-HGUT-00074.tsv",
             type = 'single',
             outdir = "output")
print(args)

library(tidyverse)
library(ggforce)

# Read results tables
cat("Reading selection test results...\n")

read_sel_tests <- function(pdir = NULL, s_coef = NULL, fit = NULL){
  res <- list(pdir = NULL, s_coef = NULL, fit = NULL)
  
  if(!is.null(pdir)){
    cat("\tReading & processing p_directional...\n")
    res$pdir <- read_tsv(args$pdir,
                         col_types = cols(site_id = col_character()))
    
    # Get signed odds ratio for p_directional. This needs to be changed to
    # inclide q estimates (probably obtained from posterior)
    res$pdir <- res$pdir %>%
      mutate(signed_or = sign(n_increase - n_decrease) * p_directional / (1 - p_directional))
  }
  if(!is.null(s_coef)){
    cat("\tReading & processing selection coefficient (s)...\n")
    res$s_coef <- read_tsv(args$s_coef,
                           col_types = cols(site_id = col_character(),
                                            ref_id = col_character()))
    
    # Process s coefs. If multiple individuals, do a t-test. Otherwise
    # keep without p-val or variance
    res$s_coef <- res$s_coef %>%
      group_by(site_id) %>%
      summarise(s.mean = mean(s),
                s.sd = sd(s),
                npts = length(s),
                .groups = 'drop') %>%
      mutate(t.s = s.mean / (s.sd / sqrt(npts))) %>%
      mutate(pval = 2 * pt(abs(t.s), df = npts - 1, lower.tail = FALSE))
    
    # # Check that s.mean and t.s give equivalent sign
    # table(sign(Scoef$s.mean) == sign(Scoef$t.s), useNA = 'always')
    
  }
  if(!is.null(fit)){
    cat("\tReading & processing FIT...\n")
    res$fit <- read_tsv(args$fit,
                        col_types = cols(site_id = col_character(),
                                         ref_id = col_character()))
    
    # Process FIT test. Use Fisher method  for combining p-values
    # Probably Fisher method is wrong because it does not consider direction?
    # Fit <- Fit %>%
    #   filter(!is.na(pval)) %>%
    #   group_by(site_id) %>%
    #   summarise(n = length(pval),
    #             fisher.chi = -2 * sum(log(pval)),
    #             .groups = 'drop') %>%
    #   mutate(pval = pchisq(fisher.chi, df = 2 * n, lower.tail = FALSE))
    # Should try inverse variance method
    res$fit <- res$fit %>%
      filter(!is.na(pval)) %>%
      group_by(site_id) %>%
      summarise(y.mean = sum(Y.mean/Y.s2) / sum(1/Y.s2),
                y.var = 1 / sum(1/Y.s2),
                npts = length(Y.mean)) %>%
      mutate(z.score = y.mean / sqrt(y.var)) %>%
      mutate(pval = 2 * pnorm(abs(z.score), mean = 0, sd = 1, lower.tail = FALSE)) %>%
      arrange(pval)
    
    
  }
  

  return(res)
}


cat("Processing selection test results...\n")





# Get overall number of tests
cat("Calculate number of tests...\n")
Res <- tibble(method = c("p_directional", "s_coef", "FIT"),
              stat = rep("n_sites", 3),
              value = c(sum(!is.na(Pdir$p_directional)),
                        sum(!is.na(Scoef$s.mean)),
                        sum(!is.na(Fit$z.score)))) %>%
  bind_rows(tibble(method = c("p_directional", "s_coef", "FIT"),
                stat = rep("tested_sites", 3),
                value = c(sum(!is.na(Pdir$p_directional)),
                          sum(!is.na(Scoef$pval)),
                          sum(!is.na(Fit$z.score)))))


# Prepare output dir
if(!dir.exists(args$outdir)){
  args$fit_dir <- file.path(args$outdir, "FIT")
  args$s_dir <- file.path(args$outdir, "s_coef")
  dir.create(args$outdir)
  dir.create(args$fit_dir)
  dir.create(args$s_dir)
}

### P_directional vs s_coef
cat("Compare p_directional versus selection coefficient (s)...\n")
# Comapre p_directional vs s_coef mean
Dat <- Pdir %>%
  select(site_id, p_directional,signed_or, pdir.pts = n_patients) %>%
  inner_join(Scoef %>%
              select(site_id, s.mean, s.pts = npts, s.pval = pval),
            by = "site_id")

# Compare numbers
p1 <- ftable(pdir.pts ~ s.pts, data = Dat) %>%
  as_tibble() %>%
  mutate(s.pts = as.numeric(as.character(s.pts)),
         pdir.pts = as.numeric(as.character(pdir.pts))) %>%
  ggplot(aes(x = pdir.pts, y = s.pts)) +
  geom_abline(slope = 1) +
  geom_point(aes(size = Freq)) +
  scale_size_continuous(range=c(-1, 5)) +
  AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "npts.png")
ggsave(filename, p1, width = 5, height = 5, dpi = 150)

# Compare unsigned
p1 <- Dat %>%
  filter(!is.na(s.pval)) %>%
  ggplot(aes(x = p_directional / (1 - p_directional), y = -log10(s.pval))) +
  geom_point(aes(col = s.pts)) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_smooth(method = "lm", formula = y ~ x) +
  AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "unsigned.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 150)

# Compare signed all
p1 <-  Dat %>%
    ggplot(aes(x = signed_or, y = s.mean)) +
    geom_point(aes(col = s.pts)) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_smooth(method = "lm", formula = y ~ x) +
    AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "signed_all.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 150)

# Compare signed tested
p1 <-  Dat %>%
    filter(!is.na(s.pval)) %>%
    ggplot(aes(x = signed_or, y = -log10(s.pval) * sign(s.mean))) +
    geom_point(aes(col = s.pts)) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_smooth(method = "lm", formula = y ~ x) +
    AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "signed_tested.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 150)

cat("\tcalculating correlations...\n")
thres <- c(0, 0.5, 0.6, 0.7, 0.8)
rs <- NULL
for(u in thres){
  rs <- c(rs,
          cor(Dat$signed_or[Dat$p_directional >= u & !is.na(Dat$s.pval)],
              (-log10(Dat$s.pval) * sign(Dat$s.mean))[Dat$p_directional >= u & !is.na(Dat$s.pval)],
              method = 'spearman'))
}
Cors <- tibble(pdir_thres = thres,
               r = rs,
               method = 's_coef')

### P_directional vs FIT
cat("Compare p_directional versus FIT...\n")
# Comapre p_directional vs FIT
Dat <- Pdir %>%
  select(site_id, p_directional,signed_or, pdir.pts = n_patients) %>%
  inner_join(Fit %>%
              select(site_id, z.score, fit.pts = npts, fit.pval = pval),
            by = "site_id")

# Compare numbers
p1 <- ftable(pdir.pts ~ fit.pts, data = Dat) %>%
  as_tibble() %>%
  mutate(fit.pts = as.numeric(as.character(fit.pts)),
         pdir.pts = as.numeric(as.character(pdir.pts))) %>%
  ggplot(aes(x = pdir.pts, y = fit.pts)) +
  geom_abline(slope = 1) +
  geom_point(aes(size = Freq)) +
  scale_size_continuous(range=c(-1, 5)) +
  AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "npts.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 150)

# Compare unsigned
p1 <- Dat %>%
  # filter(!is.na(s.pval)) %>%
  ggplot(aes(x = p_directional / (1 - p_directional), y = abs(z.score))) +
  geom_point(aes(col = fit.pts)) +
  scale_color_gradient(low = "blue", high = "red") +
  geom_smooth(method = "lm", formula = y ~ x) +
  # ggforce::facet_zoom(y = abs(z.score) < 10) +
  AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "unsigned.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 150)

# Compare signed all (=tested)
p1 <-  Dat %>%
    ggplot(aes(x = signed_or, y = z.score)) +
    geom_point(aes(col = fit.pts)) +
    scale_color_gradient(low = "blue", high = "red") +
    geom_smooth(method = "lm", formula = y ~ x) +
    # ggforce::facet_zoom(y = abs(z.score) < 10) +
    AMOR::theme_blackbox()
filename <- file.path(args$s_dir, "signed_all.png")
ggsave(filename, p1, width = 6, height = 5, dpi = 150)

# Cors
cat("\tcalculating correlations...\n")
thres <- c(0, 0.5, 0.6, 0.7, 0.8)
rs <- NULL
for(u in thres){
  rs <- c(rs,
          cor(Dat$signed_or[Dat$p_directional >= u],
              Dat$z.score[Dat$p_directional >= u],
              method = 'spearman'))
}
Cors <- bind_rows(Cors,
  tibble(pdir_thres = thres,
               r = rs,
               method = 'FIT'))

cat("Writing summaries...\n")
write_tsv(Cors, file.path(args$outdir, "cors.tsv"))
write_tsv(Res, file.path(args$outdir, "ntests.tsv"))
