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

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate the Frequency Increment Test (FIT) per site",
                        "from SliM simulations"))
  
  # Positional arguments
  p <- add_argument(p, "simdir",
                    help = paste("Directory with simulation output"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--n_timepoints",
                    help = paste("Number of tiempoints to extract for FIT.",
                                 "Timepoints will be evenly spaced."),
                    type = "numeric",
                    default = 3)
  p <- add_argument(p, "--output",
                    help = paste("File to write output"),
                    type = "character",
                    default = "FIT.tsv")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  return(args)
}


args <- process_arguments()
# args <- list(simdir = "/home/sur/micropopgen/exp/2022/today2/sim_x/",
#              n_timepoints = 3,
#              output = "alt_methods/FIT.tsv")

library(tidyverse)
library(HMVAR)


#' We read all the allele frequencies
Freqs <- list.dirs(args$simdir,
                   recursive = FALSE,
                   full.names = TRUE) %>%
  map(function(pop_dir, n_timepoints = 3){
    pop_id <- basename(pop_dir)
    
    dat <- HMVAR::read_midas_abun(file.path(pop_dir, "snps_freq.txt"))
    samples <- setdiff(colnames(dat), "site_id")
    gens <- str_remove(samples, "^gen_") %>% as.numeric()
    gens <- set_names(gens, nm = samples)
    gens <- sort(gens)
    n_gens <- length(gens)
    if(n_gens < n_timepoints)
      stop("ERROR: less samples than timepoints", call. = TRUE)
    
    samples <- floor(seq(from = 1, 
                         to = n_gens,
                         length.out = n_timepoints))
    samples <- names(gens[samples])
    cat("Selecting samples", samples, "\n")
    
    
    # meta <- tibble(sample = c("gen_1", "gen_20", "gen_50"),
    #                Group = pop_id)
    # samples <-c("gen_1", "gen_20", "gen_50")
    
    
    dat <- dat %>%
      select(site_id, all_of(samples))
    
    colnames(dat) <- c("site_id", paste0(samples, ".", pop_id))
    dat
  }, n_timepoints = args$n_timepoints)


#' We create a snps_freq like table
Dat <- Freqs[[1]] %>%
  full_join(Freqs[[2]],
            by = "site_id")
for(i in 3:length(Freqs)){
  Dat <- Dat %>%
    full_join(Freqs[[i]],
              by = "site_id")
}
Dat <- Dat %>%
  mutate(across(starts_with("gen_"), replace_na, 0))

#' Create a map
map <- tibble(sample = colnames(Dat)[-1]) %>%
  mutate(meta = sample) %>%
  separate(meta, sep = "[.]", into = c("gen", "pop_id")) %>%
  separate(gen, sep = "_", into = c("discard", "time")) %>%
  transmute(sample, time = as.numeric(time), pop_id)

#' Create long format freq
Dat <- Dat %>%
  pivot_longer(-site_id,
               values_to = "freq",
               names_to = "sample") %>%
  left_join(map, by = c("sample"))

#' Calculate FIT
Res <- FIT(Dat = Dat, aggregate = TRUE)
write_tsv(Res, args$output)
