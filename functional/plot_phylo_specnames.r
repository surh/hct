library(tidyverse)
library(ape)
library(ggtree)
# setwd("/cashew/users/sur/exp/fraserv/2022/today")

args <- list(summaries = "pdir_hit_counts.tsv",
             tre = "specs_misc/species_tree.nwk",
             spec_meta = "specs_misc/species_meta.tsv",
             pdir_info = "pdir_info/",
             outdir = "phylo")

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}



tre <- read.tree(args$tre)
tre

meta <- read_tsv(args$spec_meta,
                 col_types = cols(.default = col_character()))
meta

meta <- meta %>%
  separate(Lineage, into = c("domain",
                             "phylum",
                             "class",
                             "order",
                             "family",
                             "genus",
                             "species"),
           sep = ";") %>%
  filter(spec %in% tre$tip.label) %>%
  mutate(genus = str_remove(genus, "^g__"),
         species = str_remove(species, "^s__"))

#' No duplicated species
any(duplicated(paste(meta$species)))


spec_lut <- set_names(meta$species, meta$spec)
tre$tip.label <- as.character(spec_lut[ tre$tip.label ])




#' Basic tree plot
p1 <- ggtree(tre) + geom_tiplab(align = TRUE) +
  # hexpand(0.3) +
  # theme(plot.margin = unit(c(0,0.3,0,0), "npc"))
  xlim(c(0,3))
p1
filename <- file.path(args$outdir, "phylo_specnames.png")
ggsave(filename, p1, width = 7, height = 5)
filename <- file.path(args$outdir, "phylo_specnames.svg")
ggsave(filename, p1, width = 7, height = 5)
