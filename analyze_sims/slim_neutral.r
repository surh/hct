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

args <- list(comps_dir = "2022-06-03.selection_methods/comparisons/",
             meta = "pars_per_simulation.tsv",
             nsites = "sites_per_sim.tsv",
             outdir = "neutral_output/")

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


#' Select slimulations without sites under selection
Meta <- Meta %>%
  filter(n_selected == 0)
Meta

#' The code below calculates number of false positive as desired thresholds
#' It is for demonstration purposes, need to be incorporated at compare script:w

Scoef <- read_tsv("2022-06-03.selection_methods/s_coef/slim_1.tsv",
                  col_types = cols(site_id = col_character())) %>%
  mutate(qval = p.adjust(pval, 'fdr')) %>%
  filter(!is.na(pval))
Scoef %>%
  arrange(pval)
Scoef$qval

alpha_thres <- c(0.1, 0.05, 0.01, 0.001)
or_thres <- c(1, 2, 3 , 4)

# n_fp is the critical part here
tibble(test = "s_coef_FDR",
       ntests = nrow(Scoef),
       thres = alpha_thres,
       n_fp = colSums(Scoef$qval < matrix(rep(alpha_thres, each = nrow(Scoef)),
                                          ncol = length(alpha_thres))))






