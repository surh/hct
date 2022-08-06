library(tidyverse)
# setwd("/cashew/users/sur/exp/fraserv/2022/today")


args <- list(pdir = "preHCT_posHCT/",
             midas_dir = "snps_merged/",
             pdir_thres = 0.8,
             outdir = "pdir_info")


# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

Res <- list.files(args$pdir,
           full.names = TRUE,
           recursive = FALSE) %>%
  map(function(f){
    # f <- list.files(args$pdir,
    #                 full.names = TRUE,
    #                 recursive = FALSE)[1]
    spec <- f %>% basename %>% str_remove("[.]tsv[.]gz$")
    cat(spec, "\n")
    
    pdir <- read_tsv(f, col_types = cols(site_id = col_character()))
    pdir
  
    
    info <- HMVAR::read_midas_info(file.path(args$midas_dir, 
                                             spec, 
                                             "snps_info.txt"))
    info <- HMVAR::determine_snp_effect(info)
    info
    
    
    
    pdir <- pdir %>%
      left_join(info %>%
                  select(site_id, ref_id, ref_pos, 
                         locus_type, gene_id, site_type, 
                         snp_effect),
                by = "site_id") %>%
      mutate(spec = spec)
    
    

    
    filename <- file.path(args$outdir, paste0(spec, ".tsv.gz"))
    write_tsv(pdir, filename)
    

    
    n_pops <- max(pdir$n_patients)
    
    n_hits <- sum(pdir$p_directional >= args$pdir_thres)
    n_genes_hits <- pdir %>%
      filter(!is.na(gene_id)) %>%
      group_by(gene_id) %>%
      summarise(n_hits = sum(p_directional >= args$pdir_thres),
                .groups = 'drop') %>%
      filter(n_hits > 0) %>%
      nrow
    
    n_hits_pos <- sum(pdir$p_pos >= args$pdir_thres)
    n_genes_hits_pos <- pdir %>%
      filter(!is.na(gene_id)) %>%
      group_by(gene_id) %>%
      summarise(n_hits = sum(p_pos >= args$pdir_thres),
                .groups = 'drop') %>%
      filter(n_hits > 0) %>%
      nrow
    
    
    n_hits_ns <- pdir %>%
      filter(snp_effect == "non-synonymous") %>%
      filter(p_directional >= args$pdir_thres) %>%
      nrow
    n_genes_hits_ns <- pdir %>%
      filter(snp_effect == "non-synonymous") %>%
      group_by(gene_id) %>%
      summarise(n_hits = sum(p_directional >= args$pdir_thres),
                .groups = 'drop') %>%
      filter(n_hits > 0) %>%
      nrow
      
    
    
    tibble(spec, n_pops,
           n_hits, n_hits_ns, n_hits_pos,
           n_genes_hits, n_genes_hits_ns, n_genes_hits_pos)
    
  })


write_tsv(pdir, "pdir_hit_counts.tsv")