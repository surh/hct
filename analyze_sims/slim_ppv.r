#!/usr/bin/env Rscript

# (C) Copyright 2022 Sur Herrera Paredes
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
# along with This program.  If not, see <https://www.gnu.org/licenses/>.

#' # Setup
library(tidyverse)

args <- list(comps_dir = "2022-06-16.comparisons/comparisons/",
             meta = "pars_per_simulation.tsv",
             nsites = "sites_per_sim.tsv",
             outdir = "ppv_output/")

#' We choose the methods and colors to display
# selected_tests <- c("maf_FDR", "s_coef_FDR", "FIT_FDR",
#                     "P(directional)", "P(directional,-)")
# test_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")
selected_tests <- c("maf_FDR", "s_coef_FDR", "FIT_FDR",
                    "P(directional)")
test_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")

#' Read metadata of slimulations
Meta <- read_tsv(args$meta,
                 col_types = cols(sim_id = col_character(),
                                  run_id_short = col_character(),
                                  .default = col_number()))
Meta

Nsites <-read_tsv(args$nsites,
                  col_types = cols(sim_id = col_character(),
                                   
                                   
                                   .default = col_number()))
Nsites
Meta <- Meta %>%
  left_join(Nsites %>%
              select(sim_id, n_sites, n_standing, n_selected),
            by = "sim_id")
Meta

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

#' Reading the ppvs
PPVs <- list.dirs(args$comps_dir, recursive = FALSE, full.names = TRUE) %>%
  map_dfr(function(slim_dir){
    sim_id <- basename(slim_dir)
    
    subdir <- file.path(slim_dir, "ppv")
    filename <- file.path(subdir, "ppv_selected_thres.tsv")
    
    if(file.exists(filename)){
      dat <- read_tsv(filename, col_types = cols(test = col_character(),
                                                 .default = col_number())) %>%
        mutate(sim_id = sim_id)
    }else{
      dat <- NULL
    }
    
    return(dat)
  })
PPVs


#' Calculate FPR
PPVs <- PPVs %>%
  mutate(fpr = (n_pos - n_true_pos) / n_sites)
PPVs


#' Plot PPVs
p1 <- PPVs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  # filter(n_true_pos > 0) %>%
  ggplot(aes(x = thres, y = ppv)) +
  facet_wrap(~ test, scales = "free_x") +
  # ggbeeswarm::geom_beeswarm(aes(fill = thres), size = 2, col = "black", shape = 21, cex = 0.8) +
  ggbeeswarm::geom_beeswarm(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("#bae4bc", "#7bccc4", "#43a2ca", "#0868ac",
                                "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  geom_hline(yintercept = 0.05, col = "black") +
  ylim(c(0, 1)) +
  AMOR::theme_blackbox() 
p1
filename <- file.path(args$outdir, "ppv_selected.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "ppv_selected.svg")
ggsave(filename, p1, width = 6, height = 4)

#' Plot FPR for selected
p1 <- PPVs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  # filter(n_true_pos > 0) %>%
  ggplot(aes(x = thres, y = fpr)) +
  facet_wrap(~ test, scales = "free_x") +
  # ggbeeswarm::geom_beeswarm(aes(fill = thres), size = 2, col = "black", shape = 21, cex = 0.8) +
  ggbeeswarm::geom_beeswarm(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("#bae4bc", "#7bccc4", "#43a2ca", "#0868ac",
                                "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  geom_hline(yintercept = 0.05, col = "black") +
  ylim(c(0, 1)) +
  AMOR::theme_blackbox() 
p1
filename <- file.path(args$outdir, "fpr_selected.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "fpr_selected.svg")
ggsave(filename, p1, width = 6, height = 4)

#' We see that PPV might be a bit better for s_coef at the selected coeffients,
#' but FPR is now worse in contrast to neutral slimulations. While P(directional)
#' FPR rate is not affected by the presence of sites under selection

#' Now plot number of discoveries
meds <- PPVs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  group_by(test, thres) %>%
  summarise(n_true_pos = median(n_true_pos),
            .groups = 'drop')
meds

p1 <- PPVs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  # filter(n_true_pos > 0) %>%
  ggplot(aes(x = thres, y = n_true_pos + 1)) +
  facet_wrap(~ test, scales = "free_x") +
  # ggbeeswarm::geom_beeswarm(aes(fill = thres), size = 2, col = "black", shape = 21, cex = 0.8) +
  ggbeeswarm::geom_beeswarm(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("#bae4bc", "#7bccc4", "#43a2ca", "#0868ac",
                                "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  # geom_hline(yintercept = 0.05, col = "black") +
  # geom_crossbar(data = . %>%
  #                 group_by(test, thres) %>%
  #                 summarise(n_true_pos = median(n_true_pos),
  #                           .groups = 'drop'),
  #               aes(x = thres, ymin = n_true_pos, ymax = n_true_pos)) +
  # geom_line(stat = "hline", yintercept = "median") +
  # geom_errorbar(stat = "hline", yintercept = "mean") +
  geom_crossbar(data = meds, aes(ymin = n_true_pos + 1, ymax = n_true_pos + 1), fatten = 0, size = 1) +
  scale_y_log10() +
  AMOR::theme_blackbox() 
p1
filename <- file.path(args$outdir, "ndiscoveries_selected.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "ndiscoveries_selected.svg")
ggsave(filename, p1, width = 6, height = 4)

#' We see that in more than 50% of cases, s_coef does not discover any
#' selected sites, while P(directional) makes some discoveries
#' in more than half of the cases

#' Overall PPV is the main weakness but the method outperforms or matches
#' s_coef in all the other metrics, other methods have more serious issiues


#' # Session info
sessionInfo()