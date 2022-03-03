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
  p <- arg_parser(paste("Calculate selection coefficient (s) per site",
                        "from Moran Process simulations"))
  
  # Positional arguments
  p <- add_argument(p, "simdir",
                    help = paste("Directory with simulation output"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--N",
                     help = paste("population size in Moran Process"),
                     type = "numeric",
                     default = 100)
  p <- add_argument(p, "--T",
                    help = paste("Number of Moran process steps"),
                    type = "numeric",
                    default = 100)
  p <- add_argument(p, "--x_0",
                    help = paste("Number of individuals with genotype A",
                                 "at the beginning of the simulation."),
                    type = "numeric",
                    default = 100)
  p <- add_argument("--output",
                    help = paste("File to write output"),
                    type = "character",
                    default = "s_coef.tsv")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  return(args)
}

args <- process_arguments()
# args <- list(simdir = "/home/sur/micropopgen/exp/2022/today2/sims/sim_1/",
#              # maf_changes = "output/maf_changes.tsv",
#              N = 100,
#              T = 100,
#              x_0 = 99,
#              output = "s_coef.tsv")
library(tidyverse)

# Read maf_changes
Dat <- read_tsv(file.path(args$simdir, "maf_changes.tsv"),
                          col_types = cols(site_id = col_character()))

# Calculate selection coefficient
cat("Calculating selection coefficient...\n")
s_coef <- Dat %>%
  mutate(v0 = args$x_0 / args$N) %>%
  mutate(v1 = v0 + maf_change) %>%
  mutate(s = (1 / args$T) * log( (v1 / (1 - v1)) * ((1 - v0) / v0) ))  %>%
  select(site_id, s) %>%
  filter(!is.infinite(s))

# Process s coefs. If multiple individuals, do a t-test. Otherwise
# keep without p-val or variance
cat("Aggregating selection coefficients...\n")
s_coef <- s_coef %>%
  group_by(site_id) %>%
  summarise(s.mean = mean(s),
            s.sd = sd(s),
            npts = length(s),
            .groups = 'drop') %>%
  mutate(t.s = s.mean / (s.sd / sqrt(npts))) %>%
  mutate(pval = 2 * pt(abs(t.s), df = npts - 1, lower.tail = FALSE)) 

write_tsv(s_coef, args$output)
