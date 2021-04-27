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

# Monotonic model of site maf change. Assumes sites are replicates
# within a gene
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Make monotonic model of of site maf changes."))
  
  # Positional arguments
  p <- add_argument(p, "data",
                    help = paste("Site data, must have columns freq, pt",
                                 "day", "gene_id."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--min_sites",
                     help = paste("Min number of SNV sites per gene."),
                     type = "numeric",
                     default = 5)
  p <- add_argument(p, "--iter",
                    help = "Stan iterations",
                    type = "numeric",
                    default = 3000)
  p <- add_argument(p, "--warmup",
                    help = "Stan warmup iterations",
                    type = "numeric",
                    default = 2000)
  p <- add_argument(p, "--cores",
                    help = "Stan cores",
                    type = "numeric",
                    default = 3000)
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
# args <- list(data = "test_data.tsv",
#              min_sites = 5,
#              outdir = "output/",
#              iter = 3000,
#              warmup = 2000,
#              cores = 4)
library(tidyverse)
library(brms)

# Read data
cat("Reading site data...\n")
dat <- read_tsv(args$data,
                col_types = cols(site_id = col_character(),
                                 sample = col_character(),
                                 freq = col_number(),
                                 date = col_date(),
                                 day = col_number(),
                                 ref_pos = col_number(),
                                 gene_id = col_character(),
                                 recode = col_logical()))
# dat

# Keep only sites with multiple patients
# Could be moved to process_midas.r
cat("Remove sites without multiple patients...\n")
n_pts <- dat %>%
  group_by(site_id) %>%
  summarise(n_pt = length(unique(pt)),
            .groups = 'drop')
sites <- n_pts$site_id[ n_pts$n_pt > 1 ]
dat <- dat %>%
  filter(site_id %in% sites)

cat("Count sites per gene...\n")
gene_sites <- dat %>%
  group_by(gene_id) %>%
  summarise(n_sites = length(unique(site_id)),
            n_pt = length(unique(pt)),
            n_obs = length(freq),
            start_day = min(day),
            end_day = max(day),
            start_date = min(date),
            end_date = max(date),
            .groups = 'drop')
# gene_sites

cat("Selecting genes...\n")
genes <- gene_sites$gene_id[ gene_sites$n_sites >= args$min_sites ]
dat <- dat %>%
  filter(gene_id %in% genes)

cat("Prepare output directories")
# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
  dir.create(file.path(args$outdir, "slope_ranefs"))
  # dir.create(file.path(args$outdir, "conditional_day"))
  dir.create(file.path(args$outdir, "pt_preds"))
}

cat("Writing gene info...\n")
write_tsv(gene_sites, file.path(args$outdir, "gene_info.tsv"))

cat("Checking if any genes are kept...\n")
if(length(genes) == 0){
  cat("No genes left to model...\n")
  q()
}

# Select one gene, create model, will simply be used for update
cat("Creating dummy model...\n")
d <- dat %>%
  filter(gene_id == genes[1])
m1 <- brm(freq ~ mo(day) + (1 + mo(day) | pt),
          data = d %>%
            mutate(day = factor(day, ordered = TRUE)),
          cores = args$cores,
          iter = args$iter, warmup = args$warmup,
          control = list(adapt_delta = 0.99))

# Fit model on every gene
cat("Modelling all selected Genes")
Res <- (dat %>%
  split(.$gene_id))[1] %>%
  map_dfr(function(d, m1, cores = 4,
                   iter = 3000, warmup = 2000,
                   outdir = "output"){
    # Re-fit model with new data without recompiling
    m1 <- update(m1, newdata = d %>%
                   mutate(factor(day, ordered = TRUE)),
                 recompile = FALSE,
                 cores = cores,
                 iter = iter,
                 warmup = warmup,
                 control = list(adapt_delta = 0.99))
    
    # Get gene info
    gene <- unique(d$gene_id)

    # Plot intervals for pt slopes
    pars <- names(m1$fit) 
    pars <- pars[ (pars %>% str_detect(pattern = "^r_pt")) & (pars %>% str_detect(pattern = "day]$")) ]
    p1 <- bayesplot::mcmc_intervals(m1, pars = pars, prob = 0.5, prob_outer = 0.8)
    filename <- file.path(outdir, "slope_ranefs", paste0(gene, ".jpeg"))
    ggsave(filename, p1, width = 4, height = 6, dpi = 150)
    
    # Select pts with some slope differences
    interesting_pts <- p1$data %>%
      filter(l > 0 | h < 0) %>%
      select(parameter) %>%
      unlist %>%
      as.character() %>%
      str_remove("^r_pt\\[") %>%
      str_remove(",moday\\]$")
    
    # Plot conditional effects
    # For some reason it is failing ti create the plots but without errors !!!
    # filename <- file.path(outdir, "conditional_day", paste0(gene, ".jpeg"))
    # jpeg(filename, width = 900, height = 900)
    # conditional_effects(m1, prob = 0.8)
    # dev.off()
    
    # Plot patient-wise predictions
    # Using predict
    d_pred <- tibble(day = rep(unique(d$day), times = length(unique(d$pt))),
                     pt = rep(unique(d$pt), each = length(unique(d$day))))
    p1 <- d_pred %>%
      bind_cols(predict(m1,
                        newdata = d_pred,
                        robust = TRUE,
                        probs = c(0.1, 0.5, 0.9)) %>%
                  as_tibble()) %>%
      ggplot(aes(x = day, y = Estimate, group = pt,
                 col = pt %in% interesting_pts)) +
      geom_line() +
      geom_point() +
      theme_classic() +
      ggtitle(label = gene) +
      theme(legend.position = "bottom")
    filename <- file.path(outdir, "pt_preds", paste0(gene, ".jpeg"))
    ggsave(filename, p1, width = 6, height = 6, dpi = 150)
    
    # Store results
    res <- broom.mixed::tidyMCMC(m1$fit, rhat = TRUE, ess = TRUE,
                          robust = FALSE, conf.int = TRUE,
                          conf.level = 0.8,
                          conf.method = "HPDinterval",
                          index = FALSE) %>%
      mutate(gene_id = gene)
    
    res
  }, m1 = m1,
  cores = args$cores,
  iter = args$iter,
  warmup = args$warmup,
  outdir = args$outdir)
# Res
cat("Writing model results...\n")
write_tsv(Res, file.path(args$outdir, "model_results.tsv"))
