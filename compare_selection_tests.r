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

library(argparser)

read_sel_tests <- function(pdir = NA, s_coef = NA, fit = NA){
  res <- list(pdir = NULL, s_coef = NULL, fit = NULL)
  
  if(!is.na(pdir)){
    cat("\tReading & processing p_directional...\n")
    res$pdir <- read_tsv(pdir,
                         col_types = cols(site_id = col_character()))
    
    # Get signed odds ratio for p_directional. This needs to be changed to
    # inclide q estimates (probably obtained from posterior)
    res$pdir <- res$pdir %>%
      mutate(signed_or = sign(n_increase - n_decrease) * p_directional / (1 - p_directional))
  }
  if(!is.na(s_coef)){
    cat("\tReading & processing selection coefficient (s)...\n")
    res$s_coef <- read_tsv(s_coef,
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
  if(!is.na(fit)){
    cat("\tReading & processing FIT...\n")
    res$fit <- read_tsv(fit,
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
    plot_pvs(Dat = Dat, outdir = outdir)
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

plot_pvs <- function(Dat, outdir){
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

pdir_vs_fit <- function(Pdir, Fit, outdir = "./", plot = TRUE){
  # Comapre p_directional vs FIT
  Dat <- Pdir %>%
    select(site_id, p_directional,signed_or, pdir.pts = n_patients) %>%
    inner_join(Fit %>%
                 select(site_id, z.score, fit.pts = npts, fit.pval = pval),
               by = "site_id")
  
  if(plot){
    plot_pvf(Dat = Dat, outdir = outdir)
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
  
  Cors <- tibble(pdir_thres = thres,
                 r = rs,
                 method = 'FIT')
  
  return(list(dat = Dat, cors = Cors))
}

plot_pvf <- function(Dat, outdir){
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

compare_sel_results <- function(pdir = NULL, s_coef = NULL, fit = NULL,
                                outdir = "output/"){
  # Read results tables
  cat("Reading selection test results...\n")
  Sel <- read_sel_tests(pdir = pdir, s_coef = s_coef, fit = fit)
  
  # Prepare output dir
  cat("Prepare output directory...\n")
  if(!dir.exists(outdir)){
    fit_dir <- file.path(outdir, "FIT")
    s_dir <- file.path(outdir, "s_coef")
    dir.create(outdir)
    dir.create(fit_dir)
    dir.create(s_dir)
  }
  
  cat("Calculate number of tests...\n")
  Ntests <- count_tests(Dat = Sel)
  write_tsv(Ntests, file.path(outdir, "ntests.tsv"))
  
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
  
  return(list(PvF = PvF$dat, PvS = PvS$dat))
}

process_arguments <- function(){
  p <- arg_parser(paste("Compare multiple directional selection tests"))
  
  # Positional arguments
  p <- add_argument(p, "pdir",
                    help = paste("Results from p_directional"),
                    type = "character")
  p <- add_argument(p, "s_coefr",
                    help = paste("Results from selection coefficient (s)"),
                    type = "character")
  p <- add_argument(p, "fit",
                    help = paste("Results from Frequency Increase (FIT) test"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--type",
                     help = paste("If single, 'pdir', 's_coef' and 'fit' must",
                                  "be single files. If 'multi' then they",
                                  "must be directories."),
                     type = "character",
                     default = "single")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
# args <- list(pdir = "p_directional/MGYG-HGUT-00074.tsv.gz",
#              s_coef = "s_coef/MGYG-HGUT-00074.tsv",
#              fit = "FIT/MGYG-HGUT-00074.tsv",
#              type = 'single',
#              outdir = "output")
print(args)

library(tidyverse)
library(ggforce)


if(args$type == 'single'){
  Res <- compare_sel_results(pdir = args$pdir,
                             s_coef = args$s_coef,
                             fit = args$fit,
                             outdir = args$outdir)
}else if(args$type == 'multi'){
  write("Reading files in input dirs...\n")
  pdirs <- list.files(args$pdir, full.names = TRUE, recursive = FALSE)
  pdirs <- tibble(spec = pdirs %>% basename %>% str_remove(pattern = "[.]tsv[.]gz$"),
         pdir = pdirs)
  
  scoefs <- list.files(args$s_coef, full.names = TRUE, recursive = FALSE) 
  scoefs <- tibble(spec = scoefs %>% basename %>% str_remove(pattern = "[.]tsv$"),
                   s_coef = scoefs)
  
  fits <- list.files(args$fit, full.names = TRUE, recursive = FALSE)
  fits <- tibble(spec = fits %>% basename %>% str_remove(pattern = "[.]tsv$"),
                 fit = fits)
  
  # Cat matching files
  Files <- pdirs %>%
    left_join(scoefs, by = "spec") %>%
    left_join(fits, by = "spec")
  
  cat("\t", nrow(Files), " p_directional files...\n", sep = "")
  cat("\t", sum(!is.na(Files$s_coef)), " also have s_coef file...\n", sep = "")
  cat("\t", sum(!is.na(Files$fit)), " also have FIT file...\n", sep = "")
  
  # print(Files)
  
  Files <- Files %>%
    filter(!(is.na(s_coef) & is.na(fit)))
  cat("\t", nrow(Files), " specs to analyze...\n", sep = "")
  
  # Prepare output dir
  if(!dir.exists(args$outdir)){
    dir.create(args$outdir)
  }
  
  PvS <- PvF <- NULL
  for(i in 1:nrow(Files)){
    spec <- Files$spec[i]
    cat(spec, "\n")
    Res <- compare_sel_results(pdir = Files$pdir[i],
                               s_coef = Files$s_coef[i],
                               fit = Files$fit[i],
                               outdir = file.path(args$outdir, spec))
    
    if(!is.null(Res$PvS)){
      PvS <- PvS %>%
        bind_rows(Res$PvS %>% mutate(spec = spec))
    }
    if(!is.null(Res$PvF)){
      PvF <- PvF %>%
        bind_rows(Res$PvF %>% mutate(spec = spec))
    }
  }
  
  # Overall plots
  cat("Overall analysis...\n")
  ov_outdir <- file.path(args$outdir, "overall")
  ov_sdir <- file.path(ov_outdir, "s_coef")
  ov_fitdir <- file.path(ov_outdir, "FIT")
  dir.create(ov_outdir)
  dir.create(ov_sdir)
  dir.create(ov_fitdir)
  
  # Plots
  cat("\tPlotting...\n")
  plot_pvs(PvS, ov_sdir)
  plot_pvf(PvF, ov_fitdir)
  
  # Overall cors
  cat("\tCorrelations...\n")
  thres <- c(0, 0.5, 0.6, 0.7, 0.8)
  rs <- NULL
  for(u in thres){
    rs <- c(rs,
            cor(PvS$signed_or[PvS$p_directional >= u & !is.na(PvS$s.pval)],
                (-log10(PvS$s.pval) * sign(PvS$s.mean))[PvS$p_directional >= u & !is.na(PvS$s.pval)],
                method = 'spearman'))
  }
  Cors <- tibble(pdir_thres = thres,
                 r = rs,
                 method = 's_coef')
  rs <- NULL
  for(u in thres){
    rs <- c(rs,
            cor(PvF$signed_or[PvF$p_directional >= u],
                PvF$z.score[PvF$p_directional >= u],
                method = 'spearman'))
  }
  Cors <- Cors %>%
    bind_rows(tibble(pdir_thres = thres,
                     r = rs,
                     method = 'FIT'))
  write_tsv(Cors, file.path(ov_outdir, "cors.tsv"))
  
}else{
  stop("ERROR: only type 'single' and 'multi' are supported")
}








