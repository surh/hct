#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
# This file is part of hct.
# 
# hctis free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hctis distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hct.  If not, see <https://www.gnu.org/licenses/>.
library(argparser)

run_fgsea <- function(dat, score_type = "pos", min_size = 5){
  score_type <- "pos"
  paths <- dat %>%
    select(gene_id, terms) %>%
    pmap_dfr(HMVAR:::expand_annot) %>%
    split(.$term) %>%
    map(~.x$gene_id)
  ranks <- set_names(x = dat$score, nm = dat$gene_id)
  
  res <- fgsea(paths, ranks, scoreType = score_type, minSize = min_size)
  
  return(res)
}

process_arguments <- function(){
  p <- arg_parser(paste("Script that takes output from bernoulli mix model",
                        "and performs GSEA (via fgsea) and/or manhattan",
                        "plotting of the results"))
  
  # Positional arguments
  p <- add_argument(p, "",
                    help = paste(""),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--",
                     help = paste(""),
                     type = "character",
                     default = "")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
args <- list(bern = "test_data/p_directional/MGYG-HGUT-00099.tsv.gz",
             info = "test_data/snps_info/MGYG-HGUT-00099.txt",
             annot = "test_data/core_genes/MGYG-HGUT-00099.txt",
             contig_sizes = "test_data/seq_lengths/MGYG-HGUT-00099.tsv",
             boolean = TRUE,
             bool_vals = c("accessory", "core"),
             min_size = 5,
             annot_test = TRUE,
             locus_test = FALSE,
             ns_test = FALSE,
             plot_probs = FALSE,
             OR_trans = TRUE,
             outdir = "outs/MGYG-HGUT-00099/core_genes")



library(tidyverse)
library(HMVAR)
library(fgsea)






# Read data
bern <- read_tsv(args$bern,
                 col_types = cols(site_id = col_character(),
                                  id = col_character(),
                                  p_directional = col_number()))
info <- read_tsv(args$info,
                 col_types = cols(site_id = col_character(),
                                  ref_id = col_character(),
                                  ref_pos = col_number()))

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

if(args$plot_probs){
  contig_sizes <- read_tsv(args$contig_sizes,
                           col_names = FALSE,
                           col_types = cols(X1 = col_character(),
                                            X2 = col_number()))
  
  contig_sizes <- contig_sizes %>%
    arrange(desc(X2))
  
  point_cols <- rep(c("darkblue", "skyblue"), length.out = nrow(contig_sizes))
  
  p1 <- info %>%
    select(site_id, ref_id, ref_pos) %>%
    right_join(bern %>%
                 select(site_id, p_directional),
               by = "site_id") %>%
    mutate(ref_id = factor(ref_id, levels = contig_sizes$X1)) %>%
    arrange(ref_id, ref_pos) %>%
    split(.$ref_id) %>%
    map_dfr(function(d, contig_sizes){
      if(nrow(d) == 0){
        return(NULL)
      }
      id <- unique(d$ref_id)
      ii <- which(contig_sizes$X1 == id)
      # print(ii)
      if(ii == 1){
        size <- 0
      }else{
        size <- sum(contig_sizes$X2[1:(ii - 1)])
      }
      
      d %>%
        mutate(cum_pos = ref_pos + size)
    }, contig_sizes = contig_sizes) %>%

    ggplot(aes(x = cum_pos, y = p_directional / (1 - p_directional))) +
    geom_point(aes(col = ref_id)) +
    scale_color_manual(values = point_cols) +
    xlim(c(1, sum(contig_sizes$X2)))  +
    guides(col = FALSE) +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          panel.grid = element_blank(),
          axis.title.x = element_blank())
    
  # p1 <- info %>%
  #   select(site_id, ref_id, ref_pos) %>%
  #   right_join(bern %>%
  #                select(site_id, p_directional),
  #              by = "site_id") %>%
  #   bind_rows(contig_sizes %>%
  #               group_by(ref_id = X1) %>%
  #               summarise(site_id = "dummy",
  #                         ref_pos = c(1, X2),
  #                         p_directional = 0,
  #                         .groups = 'drop')) %>%
  #   mutate(ref_id = factor(ref_id, levels = contig_sizes$X1)) %>%
  #   ggplot(aes(x = ref_pos, y = p_directional)) +
  #   facet_grid(~ ref_id, space = "free_x", scales = "free_x") +
  #   geom_point(aes(col = ref_id)) +
  #   scale_x_continuous(limits = function(lims){lims}) +
  #   guides(col = FALSE) +
  #   theme(panel.background = element_blank(),
  #         panel.border = element_blank(),
  #         strip.placement = "none",
  #         strip.background = element_blank(),
  #         strip.text = element_blank(),
  #         axis.line.x = element_line(),
  #         axis.line.y = element_line(),
  #         panel.spacing.x = unit(0, 'lines'),
  #         panel.grid = element_blank())
  # p1
  filename <- file.path(args$outdir, "bern_manhattan.png")
  ggsave(filename, p1, width = 10, height = 4, dpi = 150)
}


if(args$locus_test){
  dat <- bern %>%
    select(site_id, p_directional) %>%
    left_join(info %>%
                select(site_id, locus_type),
              by = "site_id")

  # dat$locus_type[ sample(nrow(dat), size = floor(nrow(dat)/2)) ] <- "IGR"
  if(length(unique(dat$locus_type)) > 1 && sum(table(dat$locus_type) >= args$min_size) > 1){
    if(args$OR_trans){
      res <- dat %>%
        # head(1000) %>%
        transmute(gene_id = site_id,
                  terms = locus_type,
                  score = p_directional / (1-p_directional)) %>%
        run_fgsea(score_type = "pos", min_size = args$min_size)
    }else{
      res <- dat %>%
        transmute(gene_id = site_id,
                  terms = locus_type,
                  score = p_directional) %>%
        run_fgsea(score_type = "pos", min_size = args$min_size)
    }
    filename <- file.path(args$outdir, "locus_type_test.tsv")
    res %>%
      as_tibble() %>%
      mutate(leadingEdge = leadingEdge %>%
               map_chr(~paste(.x, collapse = ","))) %>%
      arrange(pval) %>%
      write_tsv(filename)
  }else{
    warning("Not enough locus types", call. = TRUE)
  }
}

if(args$ns_test){
  dat <- info %>%
    select(site_id, major_allele, minor_allele, amino_acids) %>%
    determine_snp_effect() %>%
    select(site_id, snp_effect) %>%
    filter(!is.na(snp_effect)) %>%
    inner_join(bern %>%
                 select(site_id, p_directional),
               by = 'site_id')
  
  if(args$OR_trans){
    res <- dat %>%
      # head(1000) %>%
      transmute(gene_id = site_id,
                terms = snp_effect,
                score = p_directional / (1 - p_directional)) %>%
      run_fgsea(score_type = "pos", min_size = args$min_size)
  }else{
    res <- dat %>%
      transmute(gene_id = site_id,
                terms = snp_effect,
                score = p_directional) %>%
      run_fgsea(score_type = "pos", min_size = args$min_size)
  }
    
  filename <- file.path(args$outdir, "ns_test.tsv")
  res %>%
    as_tibble() %>%
    mutate(leadingEdge = leadingEdge %>%
             map_chr(~paste(.x, collapse = ","))) %>%
    arrange(pval) %>%
    write_tsv(filename)
}

if(args$annot_test){
  annot <- read_tsv(args$annot,
                    col_types = cols(gene_id = col_character()))
  # annot
  
  if(args$boolean){
    # annot <- read_tsv("test_data/core_genes/MGYG-HGUT-00044.txt")
    # annot
    cat("Converting boolean annotation...\n")
    annot$terms <- args$bool_vals[ as.vector(unlist(annot[,2])) + 1 ]
    annot <- annot %>%
      select(gene_id, terms)
  }
  
  dat <- bern %>%
    select(site_id, p_directional) %>%
    inner_join(info %>%
                 select(site_id, gene_id),
               by = "site_id") %>%
    group_by(gene_id) %>%
    summarise(score = max(p_directional),
              .groups = 'drop') %>%
    inner_join(annot, by = "gene_id")
  
  if(args$OR_trans){
    dat <- dat %>%
      mutate(score = score / (1 - score))
  }
    
  res <- run_fgsea(dat = dat, score_type = "pos", min_size = args$min_size) %>%
    as_tibble() %>%
    mutate(leadingEdge = leadingEdge %>%
             map_chr(~paste(.x, collapse = ","))) %>%
    arrange(pval)
  # res
  filename <- file.path(args$outdir, "annot_test.tsv")
  write_tsv(res, filename)
}








