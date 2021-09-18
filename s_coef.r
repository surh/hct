#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate selection coefficient (Feder et. al 2014)"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("Directory with output from midas_merge.py"),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("Metadata file. Needs 'sample', 'pt', and
                                 'date'"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--output",
                    help = paste("File to write results"),
                    type = "character",
                    default = "s_coef.tsv")
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
print(args)

library(tidyverse)
library(HMVAR)

#' Calculate selection coefficient
#' 
#' Selection coefficient for alleles that change frequency between two
#' time points.
#' 
#' According to https://www.genetics.org/content/196/2/509
#'
#' @param Dat A tibble or data frame. It must have
#'
#' @return
#' @export
s_coefficient <- function(Dat){
  
  Dat %>%
    dplyr::filter(freq > 0 & freq < 1) %>%
    # dplyr::filter(site_id %in% c("100035", "100073", "100087", "69024", "69024", "69027")) %>%
    split(list(.$site_id, .$pt)) %>%
    purrr::map_dfr(function(d){
      
      if(nrow(d) < 2)
        return(NULL)
      
      pt <- unique(d$pt)
      site_id <- unique(d$site_id)
      ref_id <- unique(d$ref_id)
      ref_pos <- unique(d$ref_pos)
      
      d %>%
        # filter(date == max(date) | date == min(date)) %>%
        summarise(v0 = freq[ date == min(date)],
                  v1 = freq[ date == max(date)],
                  delta.t = as.numeric(max(date) - min(date))) %>%
        mutate(s = (1 / delta.t) * log( (v1 / (1 - v1)) * ((1 - v0) / v0) )) %>%
        mutate(pt = pt, site_id = site_id, ref_id = ref_id, ref_pos = ref_pos)
    })
}


meta <- read_tsv(args$map,
                 col_types = cols(sample = col_character(),
                                  pt = col_character(),
                                  date = col_date(format = "%Y-%m-%d")))

cat("Reading data...\n")
Dat <- read_midas_data(args$midas_dir)
Dat <- match_freq_and_depth(freq = Dat$freq,
                            depth = Dat$depth,
                            info = Dat$info %>%
                              select(site_id, ref_id, ref_pos),
                            map = meta,
                            depth_thres = 1)

cat("Calculating selection coefficients (s)...")
Res <- s_coefficient(Dat)
if(nrow(Res) > 0){
  cat("Ordering & writing results...\n")
  Res %>%
    arrange(desc(abs(s))) %>%
    write_tsv(args$output)
}else{
  warning("WARNING: no site x pt combination had at least two timepoints...\n")
}

