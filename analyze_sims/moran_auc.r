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

# setwd("~/micropopgen/exp/2022/today/")
library(tidyverse)
source(file.path(this.path::this.dir(), "functions.r"))
args <- list(comps_dir = "old_moran_sim_selmethods/2022-05-17.comparisons/",
             meta = "pars_per_simulation_moran.tsv",
             # nsites = "sites_per_sim.tsv",
             outdir = "auc_output/")

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

#' # AUC summaries
#' Read all AUC files for when there was selection
AUC <- list.dirs(args$comps_dir, recursive = FALSE, full.names = TRUE) %>%
  map_dfr(function(slim_dir){
    sim_id <- basename(slim_dir)
    cat(sim_id, "\n")
    
    roc_dir <- file.path(slim_dir, "roc")
    if(dir.exists(roc_dir)){
      auc_file <- file.path(roc_dir, "roc_auc.tsv")
      AUC <- read_tsv(auc_file, col_types = cols(test = col_character(),
                                                 AUC = col_number())) %>%
        mutate(sim_id = sim_id)
    }else{
      AUC <- NULL
    }
    
    return(AUC)
  })
AUC

#' Now we plot the AUCs
#' # NOTE: This uses a slightly different graphical representation than
#' the SLiMulation. That is because the beeswarm looks really ugly with so many
#' overlapping points near the top
table(AUC$test) %>%
  as.data.frame() %>% as_tibble() %>%
  rename(test = Var1) %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) 

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  ggplot(aes(x = test, y = AUC)) +
  # geom_boxplot(outlier.colour = NA) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.5)) +
  # ggbeeswarm::geom_beeswarm(priority = "none") +
  geom_point(position = position_jitter(width = 0.3), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_overall_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_overall_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' #STILL NEED TO FIX ROC CALCULATIONS 
#' Value below should be zero
1 - max(AUC$AUC)
warning("FIX ROC")





