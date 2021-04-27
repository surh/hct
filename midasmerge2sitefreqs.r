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

# Objective is simply to do some of the processing steps which take time
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Reads data from midas merge and writes an",
                        "extended format, keeping only sample/site",
                        "combinations with multiple patients"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("Directory with results from midas merge",
                                 "for one species."),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("Mapping file."))
  
  # Read arguments
  p <- add_argument(p, "--output",
                    help = "Name of output file",
                    type = "character",
                    default = "site_data.tsv")
   
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  return(args)
}

args <- process_arguments()
# args <- list(midas_dir = "MGYG-HGUT-00099/",
#              map = "hct_map_temp.tsv")

library(tidyverse)
library(HMVAR)


# Read metadata
meta <- read_tsv(args$map)

# Read midas data for one species
cat("Reading midas data...\n")
Dat <- read_midas_data(args$midas_dir,
                       map = meta %>%
                         select(sample, Group = pt),
                       cds_only = TRUE)

# Merge depth and freq from midas data
dat <- match_freq_and_depth(freq = Dat$freq, depth = Dat$depth,
                            info = Dat$info %>%
                              select(site_id, ref_id, ref_pos, gene_id),
                            map = meta, 
                            depth_thres = 1)

# Keep only site/patient combinations with more than one day of observations
# It is kind of slow, is there a more efficient approach? could vectorize
cat("Removing sites covered by only one patient...\n")
dat <- dat %>%
  # head(10) %>%
  split(list(.$pt, .$site_id)) %>%
  map_dfr(function(d){
    if(nrow(d) > 1){
      return(d)
    }else{
      return(NULL)
    }
  })

# Recode minor alleles that don't end as dominant
cat("Recoding...\n")
Recode <- dat %>%
  group_by(site_id) %>%
  summarise(recode = mean(freq[ day == max(day) ]) < 0.5,
            .groups = 'drop')
dat <- dat %>%
  left_join(Recode, by = "site_id")
dat$freq[ dat$recode ] <- 1 - dat$freq[ dat$recode ]

cat("Writing output...\n")
write_tsv(dat, args$output)
