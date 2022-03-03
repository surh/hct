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
                        "from Moran Process simulations"))
  
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
# args <- list(simdir = "/home/sur/micropopgen/exp/2022/today2/sims/sim_1/",
#              # maf_changes = "output/maf_changes.tsv",
#              N = 100,
#              T = 100,
#              x_0 = 99,
#              output = "s_coef.tsv")
library(tidyverse)

#' Frequency Increment Test (FIT)
#' 
#' From https://www.genetics.org/content/196/2/509
#'
#' @param Dat Tibble, must have site_id, freq, pt,
#' date, ref_id & ref_pos columns
#'
#' @return A tibble. Empty if no site x pt combination has at least three
#' timepoints
#' @export
FIT <- function(Dat){
  Dat %>%
    dplyr::filter(freq != 1 & freq != 0) %>%
    split(list(.$site_id, .$pt)) %>%
    purrr::map_dfr(function(d){
      
      if(nrow(d) < 3)
        return(NULL)
      
      pt <- unique(d$pt)
      site_id <- unique(d$site_id)
      ref_id <- unique(d$ref_id)
      ref_pos <- unique(d$ref_pos)
      
      # cat(site_id, pt, "\n")
      
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

#' Title
#'
#' @param indir 
#' @param n_timepoints 
#'
#' @return
#' @export
#'
#' @examples
read_sim_timepoints <- function(indir, n_timepoints){
  
  Freq <- NULL
  for(d in list.dirs(indir, recursive = FALSE)){
    f <- HMVAR::read_midas_abun(file.path(d, "snp_freqs.txt"))
    
    timepoints <- setdiff(colnames(f), "site_id") %>% str_remove("^t") %>% as.numeric()
    timepoints
    if((max(timepoints) - min(timepoints)) == (length(timepoints) - 1) && all(timepoints == round(timepoints))){
      timepoints <- seq(from = min(timepoints), to = max(timepoints), length.out = n_timepoints)
      cat("\t>Selected timepoints", timepoints, "\n")
    }else{
      stop("ERROR: no method implemented yet.", call. = TRUE)
    }
    
    Freq <- Freq %>%
      bind_rows(f %>%
                  select(site_id, all_of(paste0("t", timepoints))) %>%
                  pivot_longer(-site_id, names_to = "date", values_to = "freq") %>%
                  mutate(date = str_remove(string = date, pattern = "^t") %>% as.numeric()) %>%
                  mutate(pt = basename(d)))
    
  }
  
  return(Freq)
}


cat("Reading timepoints...\n")
Freq <- read_sim_timepoints(indir = file.path(args$simdir, "freqs/"),
                            n_timepoints = args$n_timepoints)
cat("Calculating Frequency Increment Test (FIT)...\n")
Fit <- FIT(Dat = Freq %>%
             mutate(ref_id = NA,
                    ref_pos = NA))
Fit <- Fit %>% select(-ref_id, -ref_pos)

cat("Aggregating FIT results...\n")
Fit <- Fit %>%
  filter(!is.na(t.fi)) %>%
  filter(!is.infinite(t.fi)) %>%
  group_by(site_id) %>%
  summarise(y.mean = sum(Y.mean/Y.s2) / sum(1/Y.s2),
            y.var = 1 / sum(1/Y.s2),
            npts = length(Y.mean)) %>%
  mutate(z.score = y.mean / sqrt(y.var)) %>%
  mutate(pval = 2 * pnorm(abs(z.score), mean = 0, sd = 1, lower.tail = FALSE)) %>%
  arrange(desc(abs(z.score)))

write_tsv(Fit, args$output)