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


library(tidyverse)

#' Read genotype data in ms format
#'
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
read_ms <- function(file){
  
  ms <- read_lines(file)
  
  # Get first two lines
  # Not actually mandatory in SLiM output
  # ms_cmd <- ms[[1]]
  # ms_seeds <- ms[[2]]
  
  Res <- NULL
  Geno <- NULL
  i <- 1
  sample_i <- 0
  while(i <= length(ms)){
    # cat(i, "\n")
    if( str_detect(ms[[i]], "^//") ){
      sample_i <- sample_i + 1
      
      if(sample_i > 1){
        stop("ERROR: script only prepared for ms files with 1 sample",
             call. = TRUE)
      }
      
      # Finding sample beginning
      i <- i + 1 
      next
    }
    
    if(sample_i > 0){
      # Only if we are withing the samples region of the file
      if( str_detect(ms[[i]], "^segsites:") ){
        n_sites <- str_remove(ms[[i]], "^segsites:")
        n_sites <- as.numeric(n_sites)
      }else if( str_detect(ms[[i]], "^positions:") ){
        positions <- str_remove(ms[[i]], "^positions:")
        positions <- str_split(trimws(positions, which = "both"), " ")[[1]]
        positions <- as.numeric(positions)
        
        if(length(positions) != n_sites){
          stop("ERROR: number of positions doesn't match number of sites",
               call. = TRUE)
        }
        
        Res[[sample_i]] <- tibble(site_id = as.character(1:n_sites),
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
    }else{
      cat("Skipping line", i, "\n")
    }
    
    
    i <- i + 1
  }
  
  return(list(info = Res[[1]], freq = Geno[[1]]))
}

#' Process simulation of one population
#'
#' @param pop_dir 
#' @param n_genomes 
#' @param genome_size 
#'
#' @return
#' @export
#'
#' @examples
process_popdir <- function(pop_dir,
                           n_genomes = 10,
                           genome_size = 1e6,
                           filter = TRUE){
  # pop_dir <- file.path(args$sim_dir, "pop_1")
  
  pop_id <- basename(pop_dir)
  cat("\tProcessing population ", pop_id, "\n")
  
  # Process population directory filenames
  pop_files <- list.files(pop_dir, recursive = FALSE, full.names = T)
  ms_files <- pop_files[ basename(pop_files) %>% str_detect("[.]ms$") ]
  generations <- basename(ms_files) %>% str_remove("^gen_") %>% str_remove("[.]ms$") %>% as.numeric()
  
  
  Pop <- generations %>%
    # set_names() %>%
    map_dfr(function(g, pop_dir, n_genomes = 10, genome_size = 1e6){
      # g <- 10
      # Reading ms file for current generation sample
      cat("\t>Generation ", g, "\n")
      dat <- read_ms( file.path(pop_dir, paste0("gen_", g, ".ms")) )
      
      # Calculate derived allele frequency from sample and merge with
      # info
      dat$info <- dat$info %>%
        left_join(dat$freq[1:(n_genomes + 1)] %>%
                    pivot_longer(-site_id, names_to = "genome",
                                 values_to = "genotype") %>%
                    group_by(site_id) %>%
                    summarise(maf = mean(genotype),
                              .groups = 'drop'),
                  by = "site_id") %>%
        mutate(ref_pos = floor(genome_size * ref_pos))
      # dat$info
      
      gen_label <- paste0("gen_", g)
      
      # Read info from same generation and homogenize
      info <- HMVAR::read_midas_info(file.path(pop_dir,
                                               paste0("gen_", g, "_snp_info.txt")));
      # info
      
      # Produce output
      info %>%
        full_join(dat$info %>%
                    select(ref_id, ref_pos, maf)) %>%
        mutate(gen = gen_label,
               s_coef = as.numeric(s_coef))
    }, pop_dir = pop_dir,
    n_genomes = n_genomes,
    genome_size = genome_size)
  # Pop
  
  # Rearrange into table
  cat("\tRearranging allele frequenciess..\n")
  Pop <- Pop %>% 
    pivot_wider(id_cols = c(snp_id, ref_id, ref_pos, s_coef, m_type),
                values_from = "maf", names_from = "gen")
  # Pop
  
  if(filter)
    Pop <- remove_undetected(Pop = Pop)
  
  # # Check repeated positions
  # Pop %>%
  #   filter(ref_pos %in% as.numeric(names(which(table(Pop$ref_pos) > 1))))
  # 
  # # Check those with selection
  # Pop %>%
  #   filter(s_coef != 0)
  
  return(Pop)
}

remove_undetected <- function(Pop){
  # Discard absent sites
  cat("Discarding undetected...\n")
  undetected_ii <- Pop %>%
    select(starts_with("gen_")) %>%
    is.na %>%
    apply(1,all)
  table(undetected_ii)
  Pop %>%
    filter(!undetected_ii)
}

args <- list(start_dir = "standing_variation/",
             sim_dir = "sim_x/",
             n_genomes = 10,
             genome_size = 1e6,
             outdir = "output/",
             seed = 2308123)


set.seed(args$seed)
dat <- read_ms("standing_variation/standing_variation.ms")
# genomes <- sample(colnames(dat$freq)[-1], size = args$n_genomes)
genomes <- colnames(dat$freq)[-1]


#' Calculate Starting derived allele frequency for all sites
Info <- dat$info %>%
  left_join(dat$freq %>%
              select(site_id, all_of(genomes)) %>%
              pivot_longer(-site_id) %>%
              group_by(site_id) %>%
              summarise(maf = mean(value),
                        .groups = 'drop'),
            by = "site_id") %>%
  mutate(ref_pos = floor(ref_pos * args$genome_size),
         m_type = "m1") %>%
  rename(gen_0 = maf)
Info




# list.dirs(args$sim_dir, recursive = FALSE, full.names = T)[1] %>%
#   map(function(pop_dir){
#         
#   })

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}


list.dirs(args$sim_dir,
          recursive = FALSE,
          full.names = TRUE)[1:2] %>%
  set_names() %>%
  map(process_popdir,
      n_genomes = args$n_genomes,
      genome_size = args$genome_size,
      filter = FALSE) %>%
  map(function(Pop, Gen0, filter = FALSE){
    # Merge with starting variation
    Pop <- Pop %>%
      full_join(Gen0 %>%
                  select(ref_id, ref_pos, m_type, gen_0),
                by = c("ref_id", "ref_pos", "m_type")) %>%
      select(site_id = snp_id, ref_id, ref_pos, s_coef, m_type, gen_0, everything())
    
    # Remove undetected
    if(filter)
      Pop <- remove_undetected(Pop)
    
    Pop
  }, Gen0 = Info, filter = TRUE) %>%
  imap(function(Pop, pop_dir, base_dir){
    pop_id <- basename(pop_dir)
    outdir <- file.path(base_dir, pop_id)
    dir.create(outdir)
    
    filename <- file.path(outdir, "snps_info.txt")
    write_tsv(Pop %>% select(-starts_with("gen_")),
              filename)
    
    filename <- file.path(outdir, "snps_freq.txt")
    write_tsv(Pop %>% select(site_id,
                             starts_with("gen_")),
              filename)
    
  }, base_dir = args$outdir)









