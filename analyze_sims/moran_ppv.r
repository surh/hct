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

library(tidyverse)
source(file.path(this.path::this.dir(), "functions.r"))
args <- list(comps_dir = "old_moran_sim_selmethods/2022-05-17.comparisons/",
             meta = "pars_per_simulation_moran.tsv",
             outdir = "ppv_output/")

#' Read metadata of slimulations
Meta <- read_tsv(args$meta,
                 col_types = cols(sim_id = col_character(),
                                  .default = col_number()))
#' Number of sites under selection per simulation
Meta <- Meta %>%
  mutate(n_selected = nsites * prop_selected) 
Meta

#' Prepare output dir
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
#' FPR = FP / (FP + TN)
PPVs <- PPVs %>%
  left_join(Meta %>%
              select(sim_id, n_selected)) %>%
  mutate(fpr = (n_pos - n_true_pos) / (n_sites - n_pos))
PPVs
# PPVs <- PPVs %>%
#   mutate(fpr = (n_pos - n_true_pos) / n_sites )
# PPVs

#' Plot PPVs
p1 <- PPVs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  ggplot(aes(x = thres, y = ppv)) +
  facet_wrap(~ test, scales = "free_x") +
  # ggbeeswarm::geom_beeswarm(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5) +
  ggbeeswarm::geom_quasirandom(aes(col = thres), size = 2, cex = 0.8, alpha = 0.3) +
  # geom_point(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5,
  #            position = position_jitter(width = 0.3)) +
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
  ggbeeswarm::geom_quasirandom(aes(col = thres), size = 2, cex = 0.8, alpha = 0.3) +
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

#' s_coef & P(directional) essentially perform on par here. Though at OR of 1
#' P(directional is probably worse). Also s_coef has a very strong bimodal
#' distribution on PPV. 







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

  ggplot(aes(x = thres, y = n_true_pos + 1)) +
  facet_wrap(~ test, scales = "free_x") +
  # ggbeeswarm::geom_beeswarm(aes(fill = thres), size = 2, col = "black", shape = 21, cex = 0.8) +
  # ggbeeswarm::geom_beeswarm(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5) +
  ggbeeswarm::geom_quasirandom(aes(col = thres), size = 2, cex = 0.8, alpha = 0.3) +
  scale_color_manual(values = c("#bae4bc", "#7bccc4", "#43a2ca", "#0868ac",
                                "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  geom_hline(yintercept = 10, col = "black") +
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

#' s_coef does a bit better here. Specially at higher significance thresholds

#' Overall PPV is the main weakness but the method outperforms or matches
#' s_coef in all the other metrics, other methods have more serious issues

#' # Session info
sessionInfo()