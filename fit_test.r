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
  p <- arg_parser(paste("Perform FIT test (Feder et. al 2014)"))
  
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
                     default = "FIT.tsv")
                     
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

#' Frequency Increment Test (FIT)
#' 
#' From https://www.genetics.org/content/196/2/509
#'
#' @param Dat Tibble, must have site_id, freq, pt,
#' date, ref_id & ref_pos columns
#'
#' @return A tibble
#' @export
FIT <- function(Dat){
  Dat %>%
    dplyr::filter(freq != 1 & freq != 0) %>%
    # head(10000) %>%
    # filter(site_id %in% c("100035", "100073", "100087", "69024", "69024", "69027")) %>%
    # arrange(site_id, pt, date) %>%
    # group_by(site_id, pt) %>% 
    # summarise(Y = diff(freq))
    split(list(.$site_id, .$pt)) %>%
    purrr::map_dfr(function(d){
      
      if(nrow(d) < 3)
        return(NULL)
      
      pt <- unique(d$pt)
      site_id <- unique(d$site_id)
      ref_id <- unique(d$ref_id)
      ref_pos <- unique(d$ref_pos)
      
      d %>%
        dplyr::arrange(date) %>%
        dplyr::mutate(vprev = dplyr::lag(freq, n = 1, default = NA),
                      t = as.numeric(date - min(date))) %>%
        dplyr::mutate(tprev = dplyr::lag(t, n = 1, default = NA)) %>%
        dplyr::mutate(Y = (freq - vprev) / sqrt(2 * vprev * (1 - vprev) * (t - tprev))) %>%
        dplyr::summarise(Y.mean = mean(Y, na.rm = TRUE),
                         Y.s2 = var(Y, na.rm = TRUE),
                         L = length(Y) - 1) %>%
        dplyr::mutate(t.fi = Y.mean / sqrt(Y.s2 / L)) %>%
        dplyr::mutate(pval = 2 * (1 - pt(q = abs(t.fi), df = L - 1)),
                      pt = pt,
                      site_id = site_id,
                      ref_id = ref_id,
                      ref_pos = ref_pos)
      
    })
}


# Trying Fit test
# https://www.genetics.org/content/196/2/509


meta <- read_tsv(args$map,
                 col_types = cols(sample = col_character(),
                                  pt = col_character(),
                                  date = col_date(format = "%Y-%m-%d")))
# meta
cat("Reading data...\n")
Dat <- read_midas_data(args$midas_dir)
Dat <- match_freq_and_depth(freq = Dat$freq,
                            depth = Dat$depth,
                            info = Dat$info %>%
                              select(site_id, ref_id, ref_pos),
                            map = meta,
                            depth_thres = 1)

cat("Performing FIT test...\n")
Res <- FIT(Dat) %>%
  arrange(pval)

cat("Writing results...\n")
write_tsv(Res, args$output)