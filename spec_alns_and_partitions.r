#!/usr/bin/env Rscri
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

library(tidyverse)
library(HMVAR)

# Based on earlier version that simply created alignment file.
# This version now creates a partitions file in RAxML style

args <- list(dir = "./",
             map = "meta.txt",
             outdir = "output",
             min_maf = 0.5,
             min_samples = 5)

# Read metadata
meta <- read_tsv(args$map)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

dirs <- list.dirs(args$dir, full.names = TRUE, recursive = FALSE)
for(d in dirs){
  spec <- basename(d)
  cat(spec, "\n")
  
  # Read alns
  seqs <- read_tsv(file.path(d, "snps_alleles.txt"), 
                   col_types = cols(site_id = col_character()))
  # seqs
  sites <- seqs$site_id
  seqs <- seqs %>%
    select(-site_id)
  # seqs
  
  if(ncol(seqs) < args$min_samples){
    warn("SKIPPING. Not enough samples")
    next
  }
  
  # Change colnames (useful for beast)
  pts <- colnames(seqs)
  new_names <- set_names(paste0(meta$sample, "_", meta$pt, "_", meta$date), meta$sample)
  colnames(seqs) <- as.character(new_names[colnames(seqs)])
  
  # Remove sites missing half or more of the samples
  ii <- rowSums(!is.na(seqs)) > ( ncol(seqs) / 2 )
  seqs <- seqs[ ii, ]
  sites <- sites[ii]
  # seqs
  
  # summary(rowSums(is.na(seqs)) / ncol(seqs), decreasing = T)
  # Filter sites by number of Ns and MAF
  ii <- apply(seqs, 1, function(x, maf = 0.2, min_samples = 5){
    x <- x[!is.na(x)]
    tab <- table(x)
    
    if( length(tab) != 2 ){
      return(FALSE)
    }else if( (min(tab) >= min_samples) && ((min(tab) / length(x)) >= maf) ){
      return(TRUE)
    }else{
      return(FALSE)
    }
    
  }, maf = args$min_maf, min_samples = args$min_samples)
  seqs <- seqs[ii, ]
  sites <- sites[ii]
  
  # Preparing partition table
  info <- read_tsv(file.paths(d, "snps_info.txt"),
                   col_types = cols(site_id = col_character(),
                                    gene_id = col_character()))
  
  # info %>% print(n = 100)
  # Select sites in and rearrange
  info <- info %>%
    select(site_id, ref_id, ref_pos, locus_type, gene_id) %>%
    filter(site_id %in% sites) %>%
    arrange(ref_id, ref_pos) 
  # info %>% print(n = 100)
  
  # Rearrange alns
  ii <- match(info$site_id, sites)
  seqs <- seqs[ ii, ]
  sites <- sites[ii]
  
  if(nrow(seqs) < 1){
    warn("SKIPPING. No sites left after filtering")
    next
  }
  if(!all(info$site_id == sites)){
    stop("ERROR. Site id's don't match")
  }
  
  
  # Get aln index
  info$n <- 1:nrow(info)
  parts <- info %>%
    group_by(gene_id) %>%
    summarise(aln_start = min(n),
              aln_end = max(n),
              .groups = 'drop')
  
  # # Print first 20 genes
  # parts %>%
  #   head(20) %>%
  #   transmute(part = paste("DNA,", gene_id, "=", paste0(aln_start,"-",aln_end))) %>%
  #   write_tsv(file = "Fplautii.20genes.partitions", col_names = FALSE)
  # 
  # ii <- 1:parts$aln_end[20]
  # seqinr::write.fasta(sequences = seqs[ii, ] %>%
  #                       as.list() %>%
  #                       map(function(x){
  #                         x[ is.na(x) ] <- "N"
  #                         x
  #                       }),
  #                     names = colnames(seqs), file.out = "Fplautii.20genes.fasta")
  
  
  # Write full aln
  filename <- file.path(args$outdir, paste0(spec,".partitions"))
  parts %>%
    transmute(part = paste("DNA,", gene_id, "=", paste0(aln_start,"-",aln_end))) %>%
    write_tsv(file = filename, col_names = FALSE)
  filename <- file.path(args$outdir, paste0(spec,".fasta"))
  seqinr::write.fasta(sequences = seqs %>%
                        as.list() %>%
                        map(function(x){
                          x[ is.na(x) ] <- "N"
                          x
                        }),
                      names = colnames(seqs), file.out = filename)
}








# set.seed(83)
# seqinr::write.fasta(sequences = seqs[sample(1:nrow(seqs), size = 300, replace = FALSE), ] %>%
#                       as.list() %>%
#                       map(function(x){
#                         x[ is.na(x) ] <- "N"
#                         x
#                       }),
#                     names = colnames(seqs), file.out = "Fplautii.rand300.fasta")
# 
# 
# seqinr::write.fasta(sequences = seqs[sample(1:nrow(seqs), size = 1000, replace = FALSE), ] %>%
#                       as.list() %>%
#                       map(function(x){
#                         x[ is.na(x) ] <- "N"
#                         x
#                       }),
#                     names = colnames(seqs), file.out = "Fplautii.rand1000.fasta")
# 
# 
# 
# 
# 
# 
# 
