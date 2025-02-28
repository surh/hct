#!/usr/bin/env Rscript

# Copyright (C) 2022 Sur Herrera Paredes
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Process standing variation & select sites for selection"))
  
  # Positional arguments
  p <- add_argument(p, "ms_file",
                    help = paste("File produced by MS with standing variation"),
                    type = "character")
  p <- add_argument(p, "genome_size",
                    help = "Genome size",
                    type = "numeric")
  
  # Optional arguments
  p <- add_argument(p, "--prop_selected",
                     help = paste("Proportion of segregating sites to be used",
                                  "to simulate selection"),
                     type = "numeric",
                     default = 0.01)
  p <- add_argument(p, "--min_maf",
                    help = paste("Minimul minor allele frequency for a site",
                                 "to be under selection"),
                    type = "numeric",
                    default = 0.05)
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
                    
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$min_maf < 0 || args$min_maf > 1){
    stop("ERROR: --min_maf must be in the interval [0,1]", call. = TRUE)
  }
  
  return(args)
}

params <- process_arguments()
# params <- list(ms_file = "standing_variation/standing_variation.ms",
#                genome_size = 1e6,
#                prop_selected = 0.01,
#                min_maf = 0.05,
#                outdir = "standing_variation/")
print(params)

library(tidyverse)
ms <- read_lines(params$ms_file)

# Get first two lines
ms_cmd <- ms[[1]]
ms_seeds <- ms[[2]]


Res <- NULL
Geno <- NULL
i <- 3
sample_i <- 0
while(i <= length(ms)){


  if( str_detect(ms[[i]], "^//") ){
    sample_i <- sample_i + 1

    if(sample_i > 1){
      stop("ERROR: script only prepared for ms files with 1 sample",
           call. = TRUE)
    }

  }else if( str_detect(ms[[i]], "^segsites:") ){
    n_sites <- str_remove(ms[[i]], "^segsites:")
    n_sites <- as.numeric(n_sites)
  }else if( str_detect(ms[[i]], "^positions:") ){
    positions <- str_remove(ms[[i]], "^positions:")
    positions <- str_split(trimws(positions, which = "both"), " ")[[1]]
    positions <- floor(as.numeric(positions) * params$genome_size)

    if(length(positions) != n_sites){
      stop("ERROR: number of positions doesn't match number of sites",
           call. = TRUE)
    }

    Res[[sample_i]] <- tibble(site_id = paste0("std.", positions),
                              ref_id = "chr1",
                              ref_pos = positions)
    Geno[[sample_i]] <- tibble(site_id = Res[[sample_i]]$site_id)

  }else if(trimws(ms[[i]]) == ""){
    cat("Skipping line", i, "\n")
  } else{
    genotype <- as.numeric(str_split(ms[[i]], "")[[1]])

    if(length(genotype) != n_sites){
      stop("ERROR: number of genotyped positions doesn't match number of sites",
           call. = TRUE)
    }

    Geno[[sample_i]][ paste0("g", ncol(Geno[[sample_i]])) ] <- genotype
  }

  i <- i + 1
}

#' Remove the list format. That is planned in case
#' we eventually want to read an ms file with multiple samples
#' but it is not needed here
Res <- Res[[1]]
Geno <- Geno[[1]]


#' Tested some code below for speed on calculating derived allele
#' frequencies (maf).
# system.time(Geno[[1]] %>%
#               select(-site_id) %>%
#               summarise_each(mean)
# )
# Geno[[1]] %>%
#   select(-site_id) %>%
#   summarise_each(mean)
#
# system.time(Geno[[1]] %>%
#               pivot_longer(-site_id, names_to = "genome", values_to = "genotype") %>%
#               group_by(site_id) %>%
#               summarise(maf = mean(genotype),
#                         .groups = 'drop'))

# Much faster
# Calculate maf
# Geno %>%
#   pivot_longer(-site_id, names_to = "genome", values_to = "genotype") %>%
#   group_by(site_id) %>%
#   summarise(maf = mean(genotype),
#             .groups = 'drop')
#' Adding derived allele frequency to info table
Res <- Res %>%
  left_join(Geno %>%
              pivot_longer(-site_id, names_to = "genome", values_to = "genotype") %>%
              group_by(site_id) %>%
              summarise(maf = mean(genotype),
                        .groups = 'drop'),
            by = "site_id")
Res

#' We plot the site frequency spectrum
Res %>%
  ggplot(aes(x = maf)) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox()

#' We pring overall minor allele frequency distribution statistics
summary(Res$maf)


#' Now we pick sites to be under selection
#' We use `ceiling` to make sure at least one site
#' is picked
n_sel_sites <- ceiling(n_sites * params$prop_selected)
n_sel_sites


#' We pick among the sites with maf above threshold
canditate_sites <- Res$site_id[ Res$maf >= params$min_maf ]
sel_sites <- sample(canditate_sites, n_sel_sites, replace = FALSE)
print(sel_sites)
Res$sel <- 1*(Res$site_id %in% sel_sites)
Res

#' Plot positions
Res %>%
  ggplot(aes(x = ref_pos, y = maf)) +
  geom_point(aes(col = as.factor(sel))) +
  AMOR::theme_blackbox()


#' Write output
# Prepare output dir
if(!dir.exists(params$outdir)){
  dir.create(params$outdir)
}

write_tsv(Res, file.path(params$outdir, "snps_info.txt"))
write_tsv(Geno, file.path(params$outdir, "snps_freq.txt"))



#' # SessionInfo
sessionInfo()
