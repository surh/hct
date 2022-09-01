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
source(file.path(this.path::this.dir(), "sim_functions.r"))

process_arguments <- function(){
  p <- arg_parser(paste("Script to simulate populations of alleles under",
                        "a modified Moran Process where a extinct allele",
                        "can be re-seeded via immigration."))
  
  # Optional arguments
  p <- add_argument(p, "--nsites",
                     help = paste("Number of sites in each population."),
                     type = "numeric",
                     default = 10000)
  p <- add_argument(p, "--npops",
                    help = paste("Number of populations to simulate with the",
                                 "same parameters"),
                    type = "numeric",
                    default = 10)
  p <- add_argument(p, "--pop_size",
                    help = paste("Population size in Moran Process"),
                    type = "numeric",
                    default = 100)
  p <- add_argument(p, "--p_imm",
                    help = paste("Probability of re-seeding locally extinct",
                                 "genotype"),
                    type = "numeric",
                    default = 0.05)
  p <- add_argument(p, "--x_0",
                    help = paste("Starting number of individuals with genotype",
                                 "A"),
                    type = "numeric",
                    default = 99)
  p <- add_argument(p, "--f",
                    help = paste("Fitness at not selected sites & of genotype",
                                 "A at selected sites"),
                    type = "numeric",
                    default = 1)
  p <- add_argument(p, "--g",
                    help = paste("Fitness of genotype B at selected sites"),
                    type = "numeric",
                    default = 3)
  p <- add_argument(p, "--T",
                    help = paste("Steps in the Moran process"),
                    type = "numeric",
                    default = 100)
  p <- add_argument(p, "--prop_selected",
                    help = paste("Proportion of sites under selection"),
                    type = "numeric",
                    default = 0.01)
  p <- add_argument(p, "--seed",
                    help = paste("Seed for numeric simulations"),
                    type = "numeric",
                    default = NULL)
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
                    
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$nsites < 1 || (args$nsites - round(args$nsites) != 0) ){
    stop("ERROR: nsites must a positive integer", call. = TRUE)
  }
  if(args$npops < 1 || (args$npops - round(args$npops) != 0) ){
    stop("ERROR: npops must a positive integer", call. = TRUE)
  }
  if(args$pop_size < 1 || (args$pop_size - round(args$pop_size) != 0) ){
    stop("ERROR: pop_size must a positive integer", call. = TRUE)
  }
  if(args$T < 1 || (args$T - round(args$T) != 0) ){
    stop("ERROR: T must a positive integer", call. = TRUE)
  }
  if(args$p_imm < 0 || args$p_imm > 1){
    stop("ERROR: p_imm must be a probability [0,1]", call. = TRUE)
  }
  if(args$prop_selected < 0 || args$prop_selected > 1){
    stop("ERROR: prop_selected must be in the range [0,1]", call. = TRUE)
  }
  if(args$x_0 > args$pop_size){
    stop("ERROR: There cannot be more individuals than the population size",
         call. = TRUE)
  }
  if(is.na(args$seed)){
    args$seed <- NULL
  }
  
  return(args)
}


args <- process_arguments()
# args <- list(nsites = 100,
#              npops = 1,
#              pop_size = 100,
#              p_imm = 0.05,
#              x_0 = 99,
#              f = 1,
#              g = 3,
#              T = 100,
#              prop_selected = 0.01,
#              seed = NULL,
#              outdir = "output")
print(args)
library(tidyverse)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
  args$freqs_dir <- file.path(args$outdir, "freqs")
  dir.create(args$freqs_dir)
}


if(is.null(args$seed)){
  warning("WARNING: No seed was provided. It won't be possible to replicate this results")
}else{
  set.seed(args$seed)
}

maf_changes <- NULL
for(p in 1:args$npops){
  cat("Simulating pop", p, "\n")
  Dat <- sim_population(nsites = args$nsites,
                        pop_size = args$pop_size,
                        x_0 = args$x_0,
                        p_imm = args$p_imm,
                        f = args$f,
                        g = args$g,
                        prop_selected = args$prop_selected,
                        T = args$T)
  
  maf_changes <- maf_changes %>%
    rbind(get_site_maf_changes(Dat$freq, t_0 = 0, t_n = args$T) %>%
            mutate(pop_id = p))
  
  # Writing allele frequencies
  popdir <- file.path(args$freqs_dir, paste0("pop_", p))
  dir.create(popdir)
  write_tsv(Dat$freq, file.path(popdir, "snp_freqs.txt"))
}

filename <- file.path(args$outdir, "maf_changes.tsv")
write_tsv(maf_changes, filename)

filename <- file.path(args$outdir, "snps_info.txt")
write_tsv(Dat$info, filename)

Sites <- maf_changes %>%
  group_by(site_id) %>%
  summarise(n_decrease = sum(maf_change < 0),
            n_equal = sum(maf_change == 0),
            n_increase = sum(maf_change >0)) %>%
  mutate(n_patients = n_decrease + n_equal + n_increase)
filename <- file.path(args$outdir, "sites.tsv")
write_tsv(Sites, filename)



