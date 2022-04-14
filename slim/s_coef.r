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
                        "from SliM simulations"))
  
  # Positional arguments
  p <- add_argument(p, "maf_changes",
                    help = paste("File with dereived allele changes"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--time",
                    help = paste("Number of SliM generations"),
                    type = "numeric",
                    default = 50)
  p <- add_argument(p, "--output",
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
# args <- list(maf_changes = "/home/sur/micropopgen/exp/2022/today2/sim_x/maf_changes.tsv",
#              time = 50,
#              output = "s_coef.tsv")
print(args)
library(tidyverse)
library(HMVAR)

# Read maf_changes
Dat <- read_tsv(args$maf_changes,
                col_types = cols(site_id = col_character(),
                                 pop_id = col_character())) %>%
  rename(v0 = t_0, v1 = t_n)

# Calculate selection coefficient
s_coef <- s_coefficient(Dat = Dat, aggregate = TRUE, time = args$time)
write_tsv(s_coef, args$output)
