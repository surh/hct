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
             outdir = "ntests_output/")

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

warning("SUMMARIES NOT CREATED YET")
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
