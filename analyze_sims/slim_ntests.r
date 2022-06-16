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
             outdir = "ntests_output/")


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

#' Merge metadata with number of sites per slimulation &
#' Clean metadata to only changing variables
Nsites <-read_tsv(args$nsites,
                  col_types = cols(sim_id = col_character(),
                                   
                                   
                                   .default = col_number()))
Nsites
Meta <- Meta %>%
  left_join(Nsites %>%
              select(sim_id, n_sites, n_standing, n_selected),
            by = "sim_id")
Meta <- Meta %>%
  select(sim_id, Mu, Rho, tractlen, scoef, 
         prop_selection, n_sites, n_standing, n_selected)


#' Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}



#' # AUC summaries
#' Read all ntests_summaries
Ntests <- list.dirs(args$comps_dir, recursive = FALSE, full.names = TRUE) %>%
  map_dfr(function(slim_dir){
    sim_id <- basename(slim_dir)
    
    subdir <- file.path(slim_dir, "summaries")
    if(dir.exists(subdir)){
      filename <- file.path(subdir, "ntests_summary.tsv")
      dat <- read_tsv(filename, col_types = cols(test = col_character(),
                                                    .default = col_number())) %>%
        mutate(sim_id = sim_id)
    }else{
      dat <- NULL
    }
    
    return(dat)
  })
Ntests <- Ntests %>%
  mutate(test = replace(test, test == "maf", "maf_FDR")) %>%
  mutate(test = replace(test, test == "s_coef", "s_coef_FDR")) %>%
  mutate(test = replace(test, test == "FIT", "FIT_FDR")) 
Ntests

#' This table indicate the total number of sites that were tested per method
#' per slimulation, the total number of sites that were under selection
#' AND could be tested (n_selected_s), the total number of sites under selection
#' in the slimulation (after accounting for sampling, n_selected_total),
#' the total number of sites in the slimulation (n_sites_total, after
#' accounting for sampling), the number of sites in the slimulation
#' (after sampling) that represented standing variation (n_sites_standing)

p1 <- Ntests %>%
  left_join(Meta %>%
              select(sim_id, Mu),
            by = "sim_id") %>%
  ggplot(aes(x = test, y = n_tested)) +
  facet_wrap(~factor(2 * 1000 * Mu), scales = "free_y",
             labeller = function(x){
               data.frame(paste0("2 * Ne * Mu=", as.character(x[,1])))
             }) +
  # geom_point(aes(col = n_selected_total > 0)) +
  ggbeeswarm::geom_beeswarm(aes(fill = n_selected_total > 0),
                            size = 2, cex = 1, col = "black", shape = 21) +
  scale_fill_discrete(name = "Selection?") +
  AMOR::theme_blackbox() +
  scale_y_log10()
p1 
filename <- file.path(args$outdir, "n_tested_all_alt.png")
ggsave(filename, p1, width = 10, height = 4)
filename <- file.path(args$outdir, "n_tested_all_alt.svg")
ggsave(filename, p1, width = 10, height = 4)


p1 <- Ntests %>%
  left_join(Meta %>%
              select(sim_id, Mu),
            by = "sim_id") %>%
  mutate(test = factor(test, selected_tests)) %>%
  ggplot(aes(x = n_selected_total > 0, y = n_tested)) +
  facet_wrap(~factor(2 * 1000 * Mu), scales = "free_y",
             labeller = function(x){
               data.frame(paste0("2 * Ne * Mu=", as.character(x[,1])))
             }) +
  # geom_point(aes(col = n_selected_total > 0)) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  scale_color_manual(values = test_colors) +
  AMOR::theme_blackbox() +
  scale_y_log10()
p1
filename <- file.path(args$outdir, "n_tested_all.png")
ggsave(filename, p1, width = 10, height = 4)
filename <- file.path(args$outdir, "n_tested_all.svg")
ggsave(filename, p1, width = 10, height = 4)

#' Overall, we see two orders of magnitude difference in the number of 
#' sites that can be tested by the other tests vs P(directional)

p1 <- Ntests %>%
  left_join(Meta %>%
              select(sim_id, Mu),
            by = "sim_id") %>%
  mutate(test = factor(test, selected_tests)) %>%
  filter(n_selected_total > 0) %>%
  ggplot(aes(x = factor(2 * 1000 * Mu), y = n_tested_s)) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2, groupOnX = TRUE) +
  scale_color_manual(values = test_colors) +
  AMOR::theme_blackbox() +
  scale_y_log10()
p1
filename <- file.path(args$outdir, "n_tested_sel.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "n_tested_sel.svg")
ggsave(filename, p1, width = 6, height = 4)

#' Confirming that MAF and P(directional) always make the same number of tests)
Ntests %>%
  group_by(sim_id) %>%
  summarise(diff = n_tested_s[ test == "maf_FDR"] == n_tested_s[ test == "P(directional)"]) %>%
  filter(!diff)

#' Printing cases where s_coef end up with a a few fewer tests
Ntests %>%
  filter(sim_id %in% (Ntests %>%
           group_by(sim_id) %>%
           summarise(diff = n_tested_s[ test == "s_coef_FDR"] == n_tested_s[ test == "P(directional)"]) %>%
           filter(!diff) %>%
           print() %>%
           select(sim_id) %>% unlist() %>% as.character())
  ) %>%
  print(n = 1000)

#' Overall most selected sites end up being included in the test, this is
#' probably because I chose selected sites to start with a minimum 5%
#' minor allele frequency in the starting population, so it is likely
#' that that is enough to ensure that the allele survives in enough population
#' to make the cut for the test most times.


