library(tidyverse)
# setwd("/cashew/users/sur/exp/fraserv/2022/today")

args <- list(spec_meta = "specs_misc/species_meta.tsv",
             pdir_info = "pdir_info/",
             uhh_catalogue = "/cashew/shared_data/mgnify/v1.0/uhgg_catalogue",
             pdir_thres = 0.8,
             outdir = "gene_hits")

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

Genes <- list.files(args$pdir_info,
           full.names = TRUE,
           recursive = FALSE) %>%
  map_dfr(function(f){
    # f <- list.files(args$pdir_info,
    #                 full.names = TRUE,
    #                 recursive = FALSE)[2]
    # f
    
    id <- basename(f) %>% str_remove("[.]tsv[.]gz$")
    cat(id, "\n")
    
    enog_file <- file.path(args$uhh_catalogue,
              str_sub(id, 1, 13),
              id, "genome",
              paste0(id, "_eggNOG.tsv"))
    
    # ipro_file <- file.path(args$uhh_catalogue,
    #                        str_sub(id, 1, 13),
    #                        id, "genome",
    #                        paste0(id, "_InterProScan.tsv"))
    
    pdir <- read_tsv(f,
                     col_types = cols(site_id = col_character(),
                                      ref_id = col_character(),
                                      spec = col_character(),
                                      site_type = col_character()))
    
    
    eggnog <- HMVAR::read_eggnog(enog_file, format = "uhgg")
    
    Res <- pdir %>%
      filter(p_directional >= args$pdir_thres)
    
    if(nrow(Res) == 0){
      return(NULL)
    }
    
    Res <- Res %>%
      group_by(gene_id) %>%
      summarise(n_hits = length(site_id),
                n_hits_ns = sum(snp_effect == "non-synonymous"),
                n_hits_s = sum(snp_effect == "synonymous"),
                max_pdir = max(p_directional),
                spec = unique(spec),
                .groups = 'drop') %>%
      arrange(desc(n_hits_ns),
              desc(n_hits),
              desc(max_pdir)) %>%
      left_join(eggnog %>%
                  rename(gene_id = query_name) %>%
                  select(gene_id, predicted_gene_name,eggNOG_annot, GO_terms,
                         KEGG_KOs, OGs, COG_cat, EC_terms,  KEGG_pathways,
                         KEGG_modules,
                         KEGG_reaction, KEGG_rclass, BRITE,
                         KEGG_TC, CAZy, BiGG_reactions ),
                by = "gene_id") 
    filename <- file.path(args$outdir, paste0(id, ".tsv"))
    write_tsv(Res, filename)
    
    Res
    })
Genes




Genes %>%
  filter(n_hits_ns > 0) %>%
  arrange(desc(n_hits_ns),
          desc(n_hits),
          desc(max_pdir)) %>%
  print(n = 20) %>%
  write_tsv("genes_ns_hits_with_annot.tsv")


gene_ogs <- Genes %>%
  filter(n_hits_ns > 0) %>%
  select(gene_id, terms = OGs) %>%
  # transmute(terms = is.na(term)) %>% table
  pmap(HMVAR:::expand_annot) %>%
  map_dfr(function(d){
    d <- d %>%
      mutate(term = str_remove(term, "@\\w+$"))
    
    d[!duplicated(d),]
  })
gene_ogs  %>%
  left_join(Genes %>%
              select(gene_id, spec,
                     n_hits, n_hits_ns, n_hits_s, max_pdir),
            by = "gene_id")
sort(table(gene_ogs$term), decreasing = TRUE) [1:10]

  