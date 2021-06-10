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
  p <- add_argument(p, "bern",
                    help = paste("TSV file with results from bernoulli mix",
                                 "model. It should have site_id and",
                                 "p_directional columns"),
                    type = "character")
  p <- add_argument(p, "info",
                    help = paste("snps_info.txt in MIDAS format (from",
                                 "midas_merge.py. It should have an entry for",
                                 "every SNP tested"))
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--manhattan",
                     help = paste("Flag indicating whether to make a",
                                  "manhattan plot of p_directional,"),
                     flag = TRUE)
  p <- add_argument(p, "--contig_sizes",
                    help = paste("File with contig sizes for manhattan plot."),
                    type = "character",
                    default = NULL)
  p <- add_argument(p, "--locus_test",
                    help = paste("Flag indicating whether to test if there is",
                                 "an enrichment for a particular locus type.",
                                 "If passed, there must be a locus_type column",
                                 "in the info file."),
                    flag = TRUE)
  p <- add_argument(p, "--ns_test",
                    help = paste("Flag indicating whether to test if there is",
                                 "an enrichment non-synonynmous or synonymous.",
                                 "SNPs amon p_directional values. Only coding",
                                 "sites are included in the test."),
                    flag = TRUE)
  p <- add_argument(p, "--annot_test",
                    help = paste("Flag indicating whether to test if there are",
                                 "enrichments of gene-level annotations for",
                                 "p_directional values. For each gene the max",
                                 "p_directional value within it is used. If",
                                 "passed, there must be a gene_id column in",
                                 "the info file, and a valid annotation file",
                                 "must be provided (see --annot)."),
                    flag = TRUE)
  p <- add_argumentS(p, "--annot",
                     help = paste("File with gene level annotations. It must",
                                  "be a TSV file with a gene_id column. It",
                                  "should also have a terms column with",
                                  "comma-separated annotation terms for each",
                                  "gene, unless it is a boolean annotation.",
                                  "If it is a boolean annotation, then the ",
                                  "second column must contain the annotation",
                                  "value (TRUE/FALSE)."))
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$manhattan && !file.exists(args$contig_sizes))
    stop("ERROR: If --manhattan is passed, then a valid --contig_sizes must be provided", call. = TRUE)
  if(args$annot_test && !file.exists(args$annot))
    stop("ERROR: If --annot_test is passed, then a valid --anot must be provided", call. = TRUE)
  
  return(args)
}

# args <- process_arguments()
args <- list(bern = "/home/sur/micropopgen/exp/2021/2021-06-10.gsea_hct_test/test_data/p_directional/MGYG-HGUT-00099.tsv.gz",
             info = "/home/sur/micropopgen/exp/2021/2021-06-10.gsea_hct_test/test_data/snps_info/MGYG-HGUT-00099.txt",
             annot = "/home/sur/micropopgen/exp/2021/2021-06-10.gsea_hct_test/test_data/core_genes/MGYG-HGUT-00099.txt",
             contig_sizes = "/home/sur/micropopgen/exp/2021/2021-06-10.gsea_hct_test/test_data/seq_lengths/MGYG-HGUT-00099.tsv",
             annot_test = TRUE,
             locus_test = FALSE,
             ns_test = FALSE,
             manhattan = FALSE,
             outdir = "/home/sur/micropopgen/exp/2021/2021-06-10.gsea_hct_testouts/MGYG-HGUT-00099/core_genes",
             
             bool_labs = c("accessory", "core"),
             min_size = 5,
             OR_trans = TRUE)



library(tidyverse)
library(HMVAR)
library(fgsea)

# Read data
cat("Reading data...\n")
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

if(args$manhattan){
  cat("Readign contig sizes...\n")
  contig_sizes <- read_tsv(args$contig_sizes,
                           col_names = FALSE,
                           col_types = cols(X1 = col_character(),
                                            X2 = col_number()))
  
  contig_sizes <- contig_sizes %>%
    arrange(desc(X2))
  
  point_cols <- rep(c("darkblue", "skyblue"), length.out = nrow(contig_sizes))
  
  cat("Creating manhattan plot...\n")
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
    }, contig_sizes = contig_sizes)
  
  if(args$OR_trans){
    p1 <- p1 %>% ggplot(aes(x = cum_pos, y = p_directional / (1 - p_directional)))
  }else{
    p1 <- p1 %>% ggplot(aes(x = cum_pos, y = p_directional))
  }

  p1 <- p1 +
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
  cat("Matching p_directional with locus type...\n")
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
  
  if( !all(c("gene_id", "terms") %in%  colnames(annot)) ){
    cat("Missing 'terms' column, checking if boolean...\n")
    
    if( "gene_id" %in% colnames(annot) && is.logical(unlist(annot[,2])) ){
      cat("Found boolean annotation in second column, proceeding...\n")
      
      if( length(args$bool_labs) != 2 || any(is.na(args$bool_labs)) ){
        cat("No valid boolean labels provided, using default...\n")
        args$bool_labs <- c("false", "true")
      }
      
      cat("Converting boolean annotation...\n")
      annot$terms <- args$bool_labs[ as.vector(unlist(annot[,2])) + 1 ]
      annot <- annot %>%
        select(gene_id, terms)
    }else{
      stop("ERROR: annotation file must either have a 'terms' column or a boolean second column", call. = TRUE)
    }
  }
  
  cat("Matching sites_to genes  and annotations...\n")
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








