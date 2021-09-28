#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
# This file is part of This program.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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

count_tests <- function(Dat){
  # Get overall number of tests
  Res <- tibble(method = c("p_directional", "s_coef", "FIT"),
                stat = rep("n_sites", 3),
                value = c(sum(!is.na(Dat$pdir$p_directional)),
                          sum(!is.na(Dat$s_coef$s.mean)),
                          sum(!is.na(Dat$fit$z.score)))) %>%
    bind_rows(tibble(method = c("p_directional", "s_coef", "FIT"),
                     stat = rep("tested_sites", 3),
                     value = c(sum(!is.na(Dat$pdir$p_directional)),
                               sum(!is.na(Dat$s_coef$pval)),
                               sum(!is.na(Dat$fit$z.score)))))
  
  return(Res)
}

pdir_vs_scoef <- function(Pdir, Scoef, outdir = "./", plot = TRUE){
  ### P_directional vs s_coef
  # Comapre p_directional vs s_coef mean
  Dat <- Pdir %>%
    select(site_id, p_directional,signed_or, pdir.pts = n_patients) %>%
    inner_join(Scoef %>%
                 select(site_id, s.mean, s.pts = npts, s.pval = pval),
               by = "site_id")
  
  if(plot){
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
    filename <- file.path(outdir, "npts.png")
    ggsave(filename, p1, width = 5, height = 5, dpi = 150)
    
    # Compare unsigned
    p1 <- Dat %>%
      filter(!is.na(s.pval)) %>%
      ggplot(aes(x = p_directional / (1 - p_directional), y = -log10(s.pval))) +
      geom_point(aes(col = s.pts)) +
      scale_color_gradient(low = "blue", high = "red") +
      geom_smooth(method = "lm", formula = y ~ x) +
      AMOR::theme_blackbox()
    filename <- file.path(outdir, "unsigned.png")
    ggsave(filename, p1, width = 6, height = 5, dpi = 150)
    
    # Compare signed all
    p1 <-  Dat %>%
      ggplot(aes(x = signed_or, y = s.mean)) +
      geom_point(aes(col = s.pts)) +
      scale_color_gradient(low = "blue", high = "red") +
      geom_smooth(method = "lm", formula = y ~ x) +
      AMOR::theme_blackbox()
    filename <- file.path(outdir, "signed_all.png")
    ggsave(filename, p1, width = 6, height = 5, dpi = 150)
    
    # Compare signed tested
    p1 <-  Dat %>%
      filter(!is.na(s.pval)) %>%
      ggplot(aes(x = signed_or, y = -log10(s.pval) * sign(s.mean))) +
      geom_point(aes(col = s.pts)) +
      scale_color_gradient(low = "blue", high = "red") +
      geom_smooth(method = "lm", formula = y ~ x) +
      AMOR::theme_blackbox()
    filename <- file.path(outdir, "signed_tested.png")
    ggsave(filename, p1, width = 6, height = 5, dpi = 150)
  }
  
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
  
  return(list(dat = Dat, cors = Cors))
}



pdir_vs_fit <- function(Pdir, Fit, outdir = "./", plot = TRUE){
  # Comapre p_directional vs FIT
  Dat <- Pdir %>%
    select(site_id, p_directional,signed_or, pdir.pts = n_patients) %>%
    inner_join(Fit %>%
                 select(site_id, z.score, fit.pts = npts, fit.pval = pval),
               by = "site_id")
  
  if(plot){
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
    filename <- file.path(outdir, "npts.png")
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
    filename <- file.path(outdir, "unsigned.png")
    ggsave(filename, p1, width = 6, height = 5, dpi = 150)
    
    # Compare signed all (=tested)
    p1 <-  Dat %>%
      ggplot(aes(x = signed_or, y = z.score)) +
      geom_point(aes(col = fit.pts)) +
      scale_color_gradient(low = "blue", high = "red") +
      geom_smooth(method = "lm", formula = y ~ x) +
      # ggforce::facet_zoom(y = abs(z.score) < 10) +
      AMOR::theme_blackbox()
    filename <- file.path(outdir, "signed_all.png")
    ggsave(filename, p1, width = 6, height = 5, dpi = 150)
  }
  
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
  
  return(list(dat = Dat, cors = Cors))
}

compare_sel_results <- function(pdir = NULL, s_coef = NULL, fit = NULL,
                                outdir = "output/"){
  # Read results tables
  cat("Reading selection test results...\n")
  Sel <- read_sel_tests(pdir = pdir, s_coef = s_coef, fit = fit)
  
  cat("Calculate number of tests...\n")
  Ntests <- count_tests(Dat = Sel)
  write_tsv(Ntests, file.path(outdir, "ntests.tsv"))
  
  # Prepare output dir
  cat("Prepare output directory...\n")
  if(!dir.exists(outdir)){
    fit_dir <- file.path(outdir, "FIT")
    s_dir <- file.path(outdir, "s_coef")
    dir.create(outdir)
    dir.create(fit_dir)
    dir.create(s_dir)
  }
  
  # Compare tests
  Cors <- PvS <- PvF <- NULL
  if(!is.null(Sel$pdir) && !is.null(Sel$s_coef)){
    cat("Compare p_directional versus selection coefficient (s)...\n")
    PvS <- pdir_vs_scoef(Pdir = Sel$pdir, Scoef = Sel$s_coef,
                         outdir = s_dir, plot = TRUE)
    
    Cors <- Cors %>%
      bind_rows(PvS$cors)
  }
  if(!is.null(Sel$pdir) && !is.null(Sel$fit)){
    cat("Compare p_directional versus FIT...\n")
    PvF <- pdir_vs_fit(Pdir = Sel$pdir, Fit = Sel$fit,
                       outdir = fit_dir, plot = TRUE)
    Cors <- Cors %>%
      bind_rows(PvF$cors)
  }
  
  if(!is.null(Cors)){
    cat("Writing summaries...\n")
    write_tsv(Cors, file.path(outdir, "cors.tsv"))
  }
  
  return(list(Pvf = PvF, PvS = PvS))
}

args <- list(pdir = "p_directional/MGYG-HGUT-00074.tsv.gz",
             s_coef = "s_coef/MGYG-HGUT-00074.tsv",
             fit = "FIT/MGYG-HGUT-00074.tsv",
             type = 'single',
             outdir = "output")
print(args)

library(tidyverse)
library(ggforce)


if(args$type == 'single'){
  Res <- compare_sel_results(pdir = args$pdir,
                             s_coef = args$s_coef,
                             fit = args$fit,
                             outdir = args$outdir)
}else if(args$type == 'multi'){
  
}else{
  stop("ERROR: only type 'single' and 'multi' are supported")
}








