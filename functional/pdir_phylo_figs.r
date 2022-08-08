library(tidyverse)
library(ape)
library(ggtree)
# setwd("/cashew/users/sur/exp/fraserv/2022/today")

args <- list(summaries = "pdir_hit_counts.tsv",
             tre = "specs_misc/species_tree.nwk",
             spec_meta = "specs_misc/species_meta.tsv",
             pdir_info = "pdir_info/",
             outdir = "phylo_figs")

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}


Dat <- read_tsv(args$summaries,
                col_types = cols(spec = col_character()))
Dat


tre <- read.tree(args$tre)
tre


Dat$spec == tre$tip.label

#' Basic tree plot
p1 <- ggtree(tre) + geom_tiplab(align = TRUE) +
  # hexpand(0.3) +
  # theme(plot.margin = unit(c(0,0.3,0,0), "npc"))
  xlim(c(0,2.7))
p1
# dev.off()

#' Plot number of populations
p2 <- Dat %>%
  select(spec, n_pops)  %>%
  rename(label = spec) %>%
  ggplot(aes(y = n_pops, x = label)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank())
p2

#' Plot of synonymous to non synonymous counts
p3 <- Dat %>%
  select(spec, n_hits, n_hits_ns)  %>%
  mutate(n_hits_s = n_hits - n_hits_ns) %>%
  select(-n_hits) %>%
  pivot_longer(-spec, names_to = "type",
               values_to = "n_hits") %>%
  rename(label = spec) %>%
  ggplot(aes(y = n_hits, x = label)) +
  geom_bar(aes(fill = type), stat = "identity") +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank())
p3


p4 <- Dat %>%
  select(spec, n_hits, n_hits_ns)  %>%
  mutate(n_hits_s = n_hits - n_hits_ns) %>%
  select(-n_hits) %>%
  pivot_longer(-spec, names_to = "type",
               values_to = "n_hits") %>%
  rename(label = spec) %>%
  ggplot(aes(y = n_hits, x = label)) +
  geom_bar(aes(fill = type), stat = "identity", position = "fill") +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank())
p4

#' Combine everything with aplot
pp <- (p2 +
  theme(axis.text.y = element_blank())) %>% 
  aplot::insert_left(p1 + 
                       xlim(c(0,5))) %>%
  aplot::insert_right(p3 +
                        theme(axis.text.y = element_blank())) %>%
  aplot::insert_right(p4 +
                        theme(axis.text.y = element_blank()))
pp  
filename <- file.path(args$outdir, "phylo_n_hits.png")
ggsave(filename, pp, width = 12, height = 5)
filename <- file.path(args$outdir, "phylo_n_hits.svg")
ggsave(filename, pp, width = 12, height = 5)


#' Get ns vs s for all tested snps.
#' Ignoring non-coding snps
bg_type <- list.files(args$pdir_info,
           full.names = TRUE,
           recursive = FALSE) %>%
  map_dfr(function(f){
    # f <- list.files(args$pdir_info,
    #                 full.names = TRUE,
    #                 recursive = FALSE)[2]

    pdir <- read_tsv(f,
                     col_types = cols(site_id = col_character(),
                                      ref_id = col_character(),
                                      spec = col_character(),
                                      site_type = col_character()))
    
    pdir %>%
      filter(!is.na(p_directional)) %>%
      filter(!is.na(snp_effect)) %>%
      group_by(spec) %>%
      summarise(n_tested_ns = sum(snp_effect == "non-synonymous"),
                n_tested_s = sum(snp_effect == "synonymous"),
                .groups = 'drop')
  })
bg_type

p5 <- bg_type %>%
  select(spec, n_tested_s, n_tested_ns)  %>%
  pivot_longer(-spec, names_to = "type",
               values_to = "n_tested") %>%
  rename(label = spec) %>%
  ggplot(aes(y = n_tested, x = label)) +
  geom_bar(aes(fill = type), stat = "identity", position = "fill") +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank())
p5


Ns_enrich <- bg_type %>%
  full_join(Dat %>%
              select(spec, n_hits, n_hits_ns)  %>%
              mutate(n_hits_s = n_hits - n_hits_ns) %>%
              select(-n_hits),
            by = "spec") %>%
  pmap_dfr(function(spec, n_tested_ns, n_tested_s, n_hits_ns, n_hits_s){
    tibble(spec = spec, n_tested_ns, n_tested_s, n_hits_ns, n_hits_s) %>%
      bind_cols(fisher.test(matrix(c(n_hits_ns, n_hits_s, 
                                     n_tested_ns, n_tested_s),
                                   ncol = 2)) %>%
                  broom::tidy() %>%
                  select(-method) %>%
                  rename(OR = estimate))
    
  })
Ns_enrich
filename <- file.path(args$outdir, "ns_hits_enrich_test.tsv")
write_tsv(Ns_enrich, filename) %>%
  



p6 <- Ns_enrich %>%
  select(spec, OR, p.value)  %>%
  rename(label = spec) %>%
  ggplot(aes(y = OR, x = label)) +
  geom_bar(aes(fill = -log10(p.value)), stat = "identity") +
  geom_hline(yintercept = 1, col = "red", size = 2) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y = element_blank())
p6


#' Overall enrichment when there are hits
Ns_enrich %>%
  filter(n_hits_ns + n_hits_s > 0) %>%
  summarise(n_tested_ns = sum(n_tested_ns),
            n_tested_s = sum(n_tested_s),
            n_hits_ns = sum(n_hits_ns),
            n_hits_s = sum(n_hits_s)) %>%
  mutate(group = "with_tests") %>%
  bind_rows(Ns_enrich %>%
              summarise(n_tested_ns = sum(n_tested_ns),
                        n_tested_s = sum(n_tested_s),
                        n_hits_ns = sum(n_hits_ns),
                        n_hits_s = sum(n_hits_s)) %>%
              mutate(group = "all")) %>%
  pmap_dfr(function(n_tested_ns, n_tested_s, n_hits_ns, n_hits_s, group){
  tibble(group, n_tested_ns, n_tested_s, n_hits_ns, n_hits_s) %>%
    bind_cols(fisher.test(matrix(c(n_hits_ns, n_hits_s, 
                                   n_tested_ns, n_tested_s),
                                 ncol = 2)) %>%
                broom::tidy() %>%
                select(-method) %>%
                rename(OR = estimate))
  
})


pchisq( -2*sum(log(Ns_enrich$p.value)), 2*length(Ns_enrich$p.value), lower.tail=FALSE)
pchisq( -2*sum(log(Ns_enrich$p.value[(Ns_enrich$n_hits_ns + Ns_enrich$n_hits_s) > 0])), 2*sum((Ns_enrich$n_hits_ns + Ns_enrich$n_hits_s) > 0), lower.tail=FALSE)

pp <- (p4 +
         geom_hline(yintercept = 0.75) +
         theme(axis.text.y = element_blank())) %>% 
  aplot::insert_left(p1 + 
                       xlim(c(0,5))) %>%
  aplot::insert_right(p5 +
                        geom_hline(yintercept = 0.75) +
                        theme(axis.text.y = element_blank())) %>%
  aplot::insert_right(p6 +
                        theme(axis.text.y = element_blank()))
pp  

filename <- file.path(args$outdir, "phylo_ns_enrich_hits.png")
ggsave(filename, pp, width = 12, height = 5)
filename <- file.path(args$outdir, "phylo_ns_enrich_hits.svg")
ggsave(filename, pp, width = 12, height = 5)


