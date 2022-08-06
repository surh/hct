#!/usr/bin/env Rscript



# setwd("/cashew/users/sur/exp/fraserv/2022/today")

library(tidyverse)
library(ape)

args <- list(genome_meta = "genomes-all_metadata.tsv",
             pdir = "preHCT_posHCT/",
             tree = "uhgg_phylogenies/bac120_iqtree.nwk",
             outdir = "specs_misc")


Genome_meta <- read_tsv(args$genome_meta,
                        col_types = cols(.default = col_character()))

uhgg_specs <- list.files(args$pdir) %>%
  str_remove("[.]tsv[.]gz$")


#' We confirm that Species_rep & Mgnify accession contain the same
#' information and they correspond to the "species" that we
#' mapped via MIDAS
Genome_meta$Genome %>% unique %>% length()
Genome_meta$Species_rep %>% unique %>% length()
Genome_meta$MGnify_accession %>% unique %>% length()

flag <- FALSE
flag <- all(duplicated(Genome_meta$Species_rep) == duplicated(Genome_meta$MGnify_accession))
flag

if(!flag)
  stop("ERROR")


#' Simplify
Genome_meta <- Genome_meta %>%
  select(Genome, Species_rep, MGnify_accession,
         Lineage)
Genome_meta

#' Get representative genome & lineage per species
Specs <- Genome_meta %>%
  group_by(spec = MGnify_accession) %>%
  summarise(rep_genome = unique(Species_rep),
            Lineage = unique(Lineage),
            .groups = 'drop')
Specs



#' Now read tree
tre <- read.tree(args$tree)

genomes2specs <- Specs %>%
  filter(spec %in% uhgg_specs)
genomes2spec <- set_names(genomes2specs$spec, nm = genomes2specs$rep_genome)
genomes2spec

#' Filter tree and change labels to "species" ids
tre <- keep.tip(tre, names(genomes2spec))
tre$tip.label <- as.character(genomes2spec[ tre$tip.label ] )
plot(tre)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}


filename <- file.path(args$outdir, "species_tree.nwk")
write.tree(tre, file = filename)
filename <- file.path(args$outdir, "species_meta.tsv")
write_tsv(Specs, file = filename)






