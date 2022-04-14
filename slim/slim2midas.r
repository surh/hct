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


#' Title
#'
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
read_slimout <- function(file){
  # file <- "sim_x/pop_1/gen_1.slimout"
  
  slim <- readr::read_lines(file)
  
  mut_flag <- FALSE
  geno_flag <- FALSE
  info <- NULL
  geno <- NULL
  for(line in slim){
    if(str_detect(line, "^#"))
      next
    
    if(line == "Mutations:"){
      mut_flag <- TRUE
      geno_flag <- FALSE
      cat("\t>Detected beginning of mutations...\n")
      next
    }
    if(line == "Genomes:"){
      mut_flag <- FALSE
      geno_flag <- TRUE
      cat("\t>Detected beginning of genomes...\n")
      next
    }
    
    
    if(!mut_flag && !geno_flag){
      print(line)
      stop("ERROR: unexpected file region", call. = TRUE)
    }
    
    # line <- slim[6]
    line <- str_split(line, pattern = " ")[[1]]
    # line
    
    if(mut_flag && !geno_flag){
      line <- matrix(line, nrow = 1)
      colnames(line) <- paste0("V", 1:ncol(line))
      info <- info %>%
        bind_rows(as_tibble(line))
    }else if(!mut_flag && geno_flag){
      
      geno <- geno %>%
        bind_rows(tibble(genome_id = line[1],
                         mut_id = line[-(1:2)]))
      
    }
    
  }
  
  colnames(info) <- c("mut_id", "site_id", "m_type",
                      "ref_pos", "s_coef", "d_coef",
                      "orig", "generation", "prevalence")
  # info
  
  geno <- geno %>%
    mutate(genome_id = str_remove(genome_id, "^p[*][:]")) %>%
    mutate(genome_id = paste0("i", genome_id)) %>%
    mutate(count = 1) %>%
    pivot_wider(names_from = "genome_id",
                values_from = "count",
                values_fill = 0)
  # geno
  
  mut2site <- set_names(info$site_id, nm = info$mut_id)
  geno <- geno %>%
    mutate(mut_id = mut2site[mut_id] %>% as.character()) %>%
    rename(site_id = mut_id)
  # geno
  
  info <- info %>%
    transmute(site_id,
              ref_id = "chr1",
              ref_pos = as.numeric(ref_pos),
              m_type,
              s_coef = as.numeric(s_coef),
              generation = as.numeric(generation),
              prevalence = as.numeric(generation))
  
  
  return(list(info = info, freq = geno))
}

#' Process simulation of one population
#'
#' @param pop_dir 
#' @param n_genomes 
#'
#' @return
#' @export
#'
#' @examples
process_popdir <- function(pop_dir,
                           n_genomes = 10){
  # pop_dir <- file.path(args$sim_dir, "pop_10")
  # n_genomes <- 10
  
  pop_id <- basename(pop_dir)
  cat("Processing population ", pop_id, "\n")
  
  # Process population directory filenames
  pop_files <- list.files(pop_dir, recursive = FALSE, full.names = T)
  # ms_files <- pop_files[ basename(pop_files) %>% str_detect("[.]ms$") ]
  slimout_files <- pop_files[ basename(pop_files) %>% str_detect("[.]slimout$") ]
  generations <- basename(slimout_files) %>% str_remove("^gen_") %>%
    str_remove("[.]slimout$") %>% as.numeric()
  
  
  Pop <- generations %>%
    # set_names() %>%
    map_dfr(function(g, pop_dir, n_genomes = 10){
      # g <- 1
      # Reading ms file for current generation sample
      cat("  >Generation ", g, "\n")
      # dat <- read_ms( file.path(pop_dir, paste0("gen_", g, ".ms")) )
      dat <- read_slimout( file.path(pop_dir, paste0("gen_", g, ".slimout")) )
      
      
      # Calculate derived allele frequency from sample and merge with
      # info
      gen_label <- paste0("gen_", g)
      dat$info %>%
        left_join(dat$freq[1:(n_genomes + 1)] %>%
                    pivot_longer(-site_id, names_to = "genome",
                                 values_to = "genotype") %>%
                    group_by(site_id) %>%
                    summarise(maf = mean(genotype),
                              .groups = 'drop'),
                  by = "site_id") %>%
        mutate(gen = gen_label)
      # dat$info
      
    }, pop_dir = pop_dir,
    n_genomes = n_genomes)
  # Pop
  
  # Rearrange into table
  cat("  Rearranging allele frequenciess..\n")
  Pop <- Pop %>% 
    select(-prevalence) %>%
    pivot_wider(values_from = "maf", names_from = "gen",
                values_fill = 0)
  return(Pop)
}

# remove_undetected <- function(Pop){
#   # Discard absent sites
#   cat("Discarding undetected...\n")
#   undetected_ii <- Pop %>%
#     select(starts_with("gen_")) %>%
#     is.na %>%
#     apply(1,all)
#   table(undetected_ii)
#   Pop %>%
#     filter(!undetected_ii)
# }

args <- list(sim_dir = "sim_x/",
             n_genomes = 10,
             outdir = "output/",
             seed = 2308123)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

Pops <- list.dirs(args$sim_dir,
          recursive = FALSE,
          full.names = TRUE) %>%
  set_names() %>%
  map(process_popdir,
      n_genomes = args$n_genomes) %>%
  imap(function(Pop, pop_dir, base_dir){
    # Pop <- Changes[[1]]
    # pop_dir <- names(Changes)[1]
    # base_dir <- args$outdir
    
    pop_id <- basename(pop_dir)
    outdir <- file.path(base_dir, pop_id)
    dir.create(outdir)
    
    # Make sure de-novo mutations have unique ids
    denovo_ii <- Pop$m_type != "m1"
    Pop <- Pop %>%
      mutate(new_id = replace(site_id,
                              denovo_ii,
                              paste0(pop_id,
                                     ".",
                                     site_id[denovo_ii]))) %>%
      mutate(site_id = new_id) %>%
      select(-new_id)
      
    filename <- file.path(outdir, "snps_info.txt")
    write_tsv(Pop %>% select(-starts_with("gen_")),
              filename)
    
    filename <- file.path(outdir, "snps_freq.txt")
    write_tsv(Pop %>% select(site_id,
                             starts_with("gen_")),
              filename)
    
    
    Pop
  }, base_dir = args$outdir)

# Writing overall snp info file
Info <- Pops %>%
  imap_dfr(function(Pop, pop_dir){
    pop_id <- basename(pop_dir)

    Pop %>%
      select(-starts_with("gen_")) %>%
      mutate(pop_src = "standing") %>%
      mutate(pop_src = replace(pop_src, m_type != "m1", pop_id))
  }) 
Info <- Info [ !duplicated(Info), ] %>%
  mutate(selected = s_coef != 0)
filename <- file.path(args$outdir, "snps_info.txt")
write_tsv(Info, filename)

# Calculate changes
Changes <- Pops %>%
  imap_dfr(function(Pop, pop_dir){
    
    pop_id <- basename(pop_dir)
    
    # Calculate maf changes used for sites counts (for bern)
    generations <- colnames(Pop)
    generations <- generations[ str_detect(generations, "^gen_") ]
    generations <- str_remove(generations, "^gen_") %>% as.numeric()
    
    t_0 <- min(generations)
    t_n <- max(generations)
    samples <- paste0("gen_", c(t_0, t_n))
    
    # Select start and end and calculate maf change
    Pop <- Pop %>%
      select(site_id, all_of(samples))
    colnames(Pop) <- c("site_id", "t_0", "t_n")
    Pop %>%
      mutate(maf_change = t_n - t_0,
             pop_id = pop_id)
    
  })
# Write maf changes file, useful for s_coef method
Changes
filename <- file.path(args$outdir, "maf_changes.tsv")
write_tsv(Changes, filename)

# Calculate sites file for bern model
Sites <- Changes %>%
  group_by(site_id) %>%
  summarise(n_decrease = sum(maf_change < 0),
            n_equal = sum(maf_change == 0),
            n_increase = sum(maf_change > 0)) %>%
  mutate(n_patients = n_decrease + n_equal + n_increase)
Sites
filename <- file.path(args$outdir, "sites.tsv")
write_tsv(Sites, filename)



