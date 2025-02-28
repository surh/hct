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
  p <- arg_parser(paste("Script that takes an UHGG genome, annotations, and",
                        "MIDAS SNPs and creates a dndscv reference."))

  # Positional arguments
  p <- add_argument(p, "gff",
                    help = paste("GFF3 format with for extracting CDSs.",
                                 "Third column must be 'CDS' for coding genes.",
                                 "Ninth column must begin with ID=<gene_id>;"),
                    type = "character")
  p <- add_argument(p, "enogg",
                    help = "eggNOG annotation for the genome",
                    type = "character")
  p <- add_argument(p, "genome_fna",
                    help = "FASTA of genome assembly",
                    type = "character")
  p <- add_argument(p, "midas_dir",
                    help = "Directory with output of midas_merge.py",
                    type = "character")

  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--enogg_format",
                     help = paste("IF UHGG annotations used the  'uhgg'.",
                                  "Alternative is 'emapper1'"),
                     type = "character",
                     default = "uhgg")

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)

  # Process arguments

  return(args)
}

args <- process_arguments()
# args <- list(gff,
#              enogg,
#              genome_fna,
#              midas_dir,
#              enogg_format = "uhgg",
#              outdir)
# args <- list(gff = "genes.gff",
#              enogg = "eggnog.tsv",
#              genome_fna = "genome.fna",
#              midas_dir = "midas_dir/",
#              enogg_format = 'uhgg',
#              outdir = "out")
print(args)

library(tidyverse)
library(dndscv)
library(seqinr)

# Prepare object to save general numbers
statistics <- list()

# Build custom reference
cat("Preparing reference genome...\n")
# Read midas data
cat("\tReading MIDAS data...\n")
Dat <- HMVAR::read_midas_data(args$midas_dir, cds_only = TRUE)

# Identify positions to change in reference
cat("\tIdentifying sites to change...\n")
to_change <- Dat$info %>%
  filter(ref_allele != major_allele) %>%
  select(ref_id, ref_pos, ref_allele, major_allele, minor_allele) %>%
  mutate(major_allele = tolower(major_allele),
         minor_allele = tolower(minor_allele),
         ref_allele = tolower(ref_allele))

# Read fasta
cat("\tReading original genome fasta...\n")
seq <- read.fasta(args$genome_fna)

# Change positions
cat("\tChanging sites...\n")
for(i in 1:nrow(to_change)){
  seq[[to_change$ref_id[i]]][to_change$ref_pos[i]] <- to_change$major_allele[i]
}

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Write custom ref for dnds
cat("\tWriting new reference...\n")
ref_fna <- file.path(args$outdir, "reference.fna")
write.fasta(seq, names = names(seq), file.out = ref_fna)

# Basic genome stats
cat("Collecting basic genome stats\n")
statistics$n_contigs <- length(seq)
statistics$genome_size <- seq %>% map_int(length) %>% sum
statistics$ambiguous_bases <- seq %>% map_int(~ sum(!(.x %in% c("a", "c", "g", "t")))) %>% sum
statistics$changed_ref_bases <- nrow(to_change)
# statistics

# Read gene defintions
cat("Reading gene defintions from GFF3 file...\n")
gff <- read_tsv(args$gff,
                col_names = FALSE,
                col_types = cols(X1 = col_character(),
                                 X3 = col_character(),
                                 X4 = col_number(),
                                 X5 = col_number(),
                                 X7 = col_character(),
                                 X9 = col_character()))
cat("\tExtracting gene IDs...\n")
gff$gene.id <- gff$X9 %>%
  map_chr(function(x){
    id <- str_split(x, pattern = ";")[[1]][1]
    id <- str_remove(id, "^ID=")
    id
  })
# gff

statistics$n_feats <- nrow(gff)
statistics$n_cds <- sum(gff$X3 == "CDS")
# statistics

cat("\tDeleting features with ambiguous bases...\n")
feats_to_remove <- NULL
for(i in 1:nrow(gff)){
  gene_seq <- seq[[gff$X1[i]]][ gff$X4[i]:gff$X5[i] ]
  if(any(!(gene_seq %in% c("a", "c", "g", "t")))){
    cat("\t\tRemoving feature", gff$gene.id[i], "\n")
    feats_to_remove <- c(feats_to_remove, gff$gene.id[i])
  }
}
gff <- gff %>%
  filter(!(gene.id %in% feats_to_remove))


statistics$n_feats_noamb <- nrow(gff)
statistics$n_cds_noamb <- sum(gff$X3 == "CDS")
# statistics

# Annotations provide gene names
cat("Reading eggNOG annotations...\n")
annot <- HMVAR::read_eggnog(args$enogg,
                            format = args$enogg_format)


statistics$n_annotated <- nrow(annot)
statistics$n_named <- sum(!is.na(annot$predicted_gene_name))

# Create cds object for cds file
cat("Creating CDS file...\n")
cds <- gff %>%
  filter(X3 == "CDS") %>%
  transmute(gene.id,
            cds.id = gene.id,
            chr = X1,
            chr.coding.start = X4,
            chr.coding.end = X5,
            strand = X7) %>%
  mutate(length = (chr.coding.end - chr.coding.start + 1)) %>%
  mutate(cds.start = 1,
         cds.end = length,
         strand_num = 1) %>%
  mutate(strand_num = replace(strand_num, strand == '-', -1)) %>%
  left_join(annot %>%
              select(gene.id = query_name, gene.name = predicted_gene_name),
            by = "gene.id") %>%
  select(gene.id, gene.name, cds.id, chr, chr.coding.start, chr.coding.end,
         cds.start, cds.end, length, strand = strand_num) %>%
  mutate(gene_num = gene.id %>%
           str_split(pattern = "_") %>%
           map_chr(~ .x[3])) %>%
  mutate(gene.name = replace(gene.name, !is.na(gene.name),
                             paste(gene.name[!is.na(gene.name)],
                                   gene_num[!is.na(gene.name)], sep = ":"))) %>%
  mutate(gene.name = replace(gene.name, is.na(gene.name),
                             gene.id[ is.na(gene.name) ])) %>%
  select(-gene_num)
cds_tsv <- file.path(args$outdir, "cdsfile.tsv")
write_tsv(cds, cds_tsv)

# Statsitics about the processing
cat("Writing statistics...\n")
filename <- file.path(args$outdir, "statistics.tsv")
write_tsv(tibble(statistic = names(statistics),
                 value = unlist(statistics)), filename)

# Build reference
cat("Building reference...\n")
buildref(cdsfile = cds_tsv,
         genomefile = ref_fna, 
         outfile = file.path(args$outdir, "reference.rda"),
         numcode = 1, excludechrs = NULL,
         onlychrs = NULL, useids = FALSE)
