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


args <- list(midas_dir = "~/micropopgen/exp/2021/2021-08-20.dnds/Fplautii_merged/",
             refdb = "~/micropopgen/exp/2021/2021-08-20.dnds/refdb/Fplautii.rda",
             mode = "all_dummy",
             maf_thres = 0.8,
             max_coding_muts_per_sample = Inf,
             max_muts_per_gene_per_sample = Inf,
             genetic_code = 1,
             seed = NULL,
             outdir = "output")
print(args)

library(tidyverse)
library(HMVAR)
library(dndscv)

cat("Read midas data...\n")
Dat <- HMVAR::read_midas_data(args$midas_dir, cds_only = TRUE)

if(args$mode == "random_top_thres"){
  cat(paste("random_top: Choosing one random sample per SNV among those that",
            "have the greatest minor allele frequency above some threshold.\n"))
  
  info <- Dat$info %>%
    select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                     depth_thres = 1)
  
  cat("\tThresholding & random selection...\n")
  set.seed(args$seed)
  Dat <- Dat %>%
    filter(!is.na(ref)) %>%
    filter(freq >= args$maf_thres) %>%
    split(.$site_id) %>%
    map_dfr(function(d){
      ii <- which(d$freq == max(d$freq))
      if(length(ii) > 1){
        ii <- sample(ii, size = 1)
      }
      
      d[ii, ]
    })
  # Dat
  
  # Preparing mutations object
  cat("\tPreparing mutations table...\n")
  Dat <- Dat %>%
    select(sampleID = sample, chr, pos, ref, mut)
}else if(args$mode == "all_thres"){
  cat(paste("all_thres: Include all samples with minor allele frequency above",
            "some threshold\n"))
  
  info <- Dat$info %>%
    select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                     depth_thres = 1)
  
  cat("\tThresholding...\n")
  Dat <- Dat %>%
    filter(freq >= args$maf_thres)
  
  # Preparing mutations object
  cat("\tPreparing mutations table...\n")
  Dat <- Dat %>%
    select(sampleID = sample, chr, pos, ref, mut)
}else if(args$mode == "all"){
  cat(paste("all: Include all samples with minor allele at any frequency\n"))
  
  info <- Dat$info %>%
    select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                     depth_thres = 1)
  
  cat("\tPreparing mutations table...\n")
  Dat <- Dat %>%
    select(sampleID = sample, chr, pos, ref, mut)
}else if(args$mode %in% c("true_singletons", "dummy_singletons")){
  cat(paste("true_singletons: All variants detected in only one sample",
            "at any frequency\n"))
  
  info <- Dat$info %>%
    select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                     depth_thres = 1)
  
  cat("\tFinding true singletons...\n")
  Dat <- Dat %>%
    # filter(depth > 0) %>%
    split(.$site_id) %>%
    map_dfr(function(d){
      d <- d %>%
        filter(freq > 0)
      
      if(nrow(d) == 1){
        return(d)
      }else{
        return(NULL)
      }
    })
  
  cat("\tPreparing mutations table...\n")
  if(args$mode == "true_singletons"){
    Dat <- Dat %>%
      select(sampleID = sample, chr, pos, ref, mut)
  }else if(args$mode == "dummy_singletons"){
    cat(paste("\tdummy_singletons: pretending all singletons came", 
              "from one sample.\n"))
    
    Dat <- Dat %>%
      transmute(sampleID = "dummySample", chr = chr, pos = pos,
                ref = ref, mut = mut)
  }else{
    stop("ERROR: unknown mode", call. = TRUE)
  }
}else if(args$mode == "thres_singletons"){
  cat(paste("thres_singletons: All variants detected above some threshold in",
            "only one sample.\n"))
    
  info <- Dat$info %>%
    select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                     depth_thres = 1)
  
  cat("\tFinding thresholded singletons...\n")
  Dat <- Dat %>%
    split(.$site_id) %>%
    map_dfr(function(d){
      d <- d %>%
        filter(freq >= args$maf_thres)
      
      if(nrow(d) == 1){
        return(d)
      }
      
      return(NULL)
    })
  
  cat("\tPreparing mutations table...\n")
  Dat <- Dat %>%
    select(sampleID = sample, chr, pos, ref, mut)
}else if(args$mode == "all_dummy"){
  cat(paste("all_dummy: All variants detected, pretend they came from one",
            "sample\n"))
  
  cat("\tPreparing mutations table...\n")
  Dat <- Dat$info %>%
    transmute(sampleID = "dummySample", chr = ref_id, pos = ref_pos,
              ref = major_allele, mut = minor_allele)
}else{
  stop("ERROR: mode not recognized")
}


cat("Running dndscv...\n")
res <- dndscv(mutations = Dat, refdb = args$refdb,
              max_coding_muts_per_sample = args$max_coding_muts_per_sample,
              max_muts_per_gene_per_sample = args$max_muts_per_gene_per_sample,
              numcode = args$genetic_code)

# Prepare output dir
cat("Creating output directory...\n")
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

cat("Writing local dnds...\n")
filename <- file.path(args$outdir, "dnds_loc.tsv")
res$sel_loc %>% as_tibble() %>%
  select(-n_spl, -wspl_loc) %>%
  # filter(qall_loc < 0.1) %>%
  arrange(desc(wmis_loc)) %>%
  write_tsv(filename)

cat("Writing cv dnds...\n")
filename <- file.path(args$outdir, "dnds_cv.tsv")
res$sel_cv %>%
  as_tibble() %>%
  select(-n_spl, -wspl_cv) %>%
  # filter(qallsubs_cv < 0.1) %>%
  arrange(desc(wmis_cv)) %>%
  write_tsv(filename)

