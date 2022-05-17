#!/usr/bin/env Rscript

# (C) Copyright 2022 Sur Herrera Paredes
# This file is part of This software.
# 
# This software is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with This software.  If not, see <https://www.gnu.org/licenses/>.

#+ read arguments, echo = FALSE
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Script to compare P(directional) detection of",
                        "selection with methods from Feder et al.",
                        "Script ready for knitr::spin()."))
  
  # Optional arguments
  p <- add_argument(p, "--s_coef",
                     help = paste("File with selection coefficient estimates."),
                     type = "character",
                     default = "")
  p <- add_argument(p, "--FIT",
                    help = paste("File with estimates from FIT."),
                    type = "character",
                    default = "")
  p <- add_argument(p, "--pdir",
                    help = paste("File with P(directional) estimates."),
                    type = "character",
                    default = "")
  p <- add_argument(p, "--info",
                    help = paste("File with ground truth for SNPs."),
                    type = "character",
                    default = "")
  p <- add_argument(p, "--outdir",
                    help = paste(""),
                    type = "character",
                    default = "output")
  
  
  ## ADD MAF
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  return(args)
}

# args <- process_arguments()
args <- list(s_coef = "/home/sur/micropopgen/exp/2022/today4/sim_results/sim_1/s_coef.tsv",
             FIT = "/home/sur/micropopgen/exp/2022/today4/sim_results/sim_1/FIT.tsv",
             pdir = "/home/sur/micropopgen/exp/2022/today4/sim_results/sim_1/p_directional.tsv.gz",
             maf = "/home/sur/micropopgen/exp/2022/today4/sim_results/sim_1/maf_changes.tsv.gz",
             info = "/home/sur/micropopgen/exp/2022/today4/sim_results/sim_1/snps_info.txt.gz",
             outdir = "comp_test/")
library(tidyverse)

#+ print args
print(args)

#+ functions, echo=FALSE
#' Title
#'
#' @param d tibble or data frame. Must have columns `truth` & `score`. The
#' former must be interpretable as logical & the latter must be interpretable
#' as numeric with higher values for stronger evidence of bein a true positive.
#' @return
#' @export
#'
#' @examples
roc <- function(d){
  d <- d %>%
    arrange(desc(score))
  # d
  
  
  roc <- NULL
  tp <- 0
  fp <- 0
  
  n_tp <- sum(d$truth)
  n_tf <- nrow(d) - n_tp
  for(i in 1:nrow(d)){
    if(d$truth[i]){
      tp <- tp + (1 / n_tp)
    }else{
      fp <- fp + (1 / n_tf )
    }
    
    roc <- roc %>%
      bind_rows(tibble(rank = i,
                       tpr = tp,
                       fpr = fp))
    
  }
  
  return(roc)
}

#' We read the datA
#+ read data
cat("Reading data...\n")
s_coef <- read_tsv(args$s_coef,
                   col_types = cols(site_id = col_character()))
FIT <- read_tsv(args$FIT,
                col_types = cols(site_id = col_character()))
pdir <- read_tsv(args$pdir,
                 col_types = cols(site_id = col_character()))
maf <- read_tsv(args$maf,
                col_types = cols(site_id = col_character()))
info <- read_tsv(args$info,
                 col_types = cols(site_id = col_character()))


#' maf in this case is per site per pop. We need to aggregate per site
#' to compare with other methods
#+ process maf
maf <- maf %>%
  group_by(site_id) %>%
  summarise(freq_change = mean(maf_change),
            freq_sd = sd(maf_change),
            n_pops = length(maf_change),
            .groups = "drop") %>%
  mutate(t.value = freq_change / freq_sd) %>%
  mutate(pval = 2*(1 - pt(abs(t.value), df = n_pops - 1)))

#' We create a directory to store output
#+ outdir
cat("Prepare output dir...\n")
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
  
  args$pvaldir <- file.path(args$outdir, "pvals")
  dir.create(args$pvaldir)
  
  args$statplot <- file.path(args$outdir, "statplot")
  dir.create(args$statplot)
  
  args$confdir <- file.path(args$outdir, "confmat")
  dir.create(args$confdir)
  
  args$rocdir <- file.path(args$outdir, "roc")
  dir.create(args$rocdir)
  
  args$ppvdir <- file.path(args$outdir, "ppv")
  dir.create(args$ppv)
  
}else{
  stop("ERROR: output directory already exists", call. = TRUE)
}

#' # Method comparison
#' ## Comparison of *p-value* distributions
#' 
#' This is just a visual comparison and only for the methods in Feder et. al.
#' 
#' The p-value distributions should be uniform-like with a spike at low p-values
#' (left hand side of the plot). Deviations from this indicate that the
#' tests are not performing as well as expected.

#+ pvalue distributions
p1 <- maf %>%
  ggplot(aes(x = pval)) +
  geom_histogram(bins = 20) +
  ggtitle(label = "Allele frequency change") +
  AMOR::theme_blackbox()
filename <- file.path(args$pvaldir, "maf_pvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- s_coef %>%
  ggplot(aes(x = pval)) +
  geom_histogram(bins = 20) +
  ggtitle(label = "Selection coefficient (s)") +
  AMOR::theme_blackbox()
filename <- file.path(args$pvaldir, "s_coef_pvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- FIT %>%
  ggplot(aes(x = pval)) +
  geom_histogram(bins = 20) +
  ggtitle(label = "Frequency Increment Test (s)") +
  AMOR::theme_blackbox()
filename <- file.path(args$pvaldir, "FIT_pvals.png")
ggsave(filename, p1, width = 5, height = 4)


#' ## Compare main statistics with ground truth
#' 
#' The Feder *et. al.* methods use p-values as their main statistic. We expect smaller
#' p-values for sites truly under selection in the simulation.
#' 
#' The P(directional) and its partitions (P(directional,+), P(directional, -))
#' produce odds ratios of the posterior probability that a site is under directional
#' selection, and so sites under true selection in the simulation must have higher
#' ORs.

#+ statistic plots, fig.cap = "Comparison of main statistics between different selection tests. For *p-value* based methods the red line is at nominal $\\alpha = 0.05$ significance & for P(directional) methods the line is at $OR = 1$"
cat("Plotting main statistic of each test vs ground truth...\n")
p1 <- maf %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = pval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  scale_y_log10() +
  geom_hline(yintercept = 0.05, col = "red") +
  ggtitle(label = "Allele frequency change") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "maf_change_pvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- maf %>%
  left_join(info, by = "site_id") %>%
  mutate(pval = replace(pval, pval == 0, min(pval[ pval > 0 ]))) %>%
  mutate(qval = p.adjust(pval, method = 'fdr')) %>%
  ggplot(aes(x = selected, y = qval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  scale_y_log10() +
  geom_hline(yintercept = 0.05, col = "red") +
  ggtitle(label = "Allele frequency change") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "maf_change_qvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- s_coef %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = pval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 0.05, col = "red") +
  scale_y_log10() +
  ggtitle(label = "Selection coefficient (s)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "s_coef_pvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- s_coef %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr')) %>%
  ggplot(aes(x = selected, y = qval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 0.05, col = "red") +
  scale_y_log10() +
  ggtitle(label = "Selection coefficient (s)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "s_coef_qvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- FIT %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = pval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 0.05, col = "red") +
  scale_y_log10() +
  ggtitle(label = "Frequency Increment Test (s)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "FIT_pvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- FIT %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr')) %>%
  ggplot(aes(x = selected, y = qval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 0.05, col = "red") +
  scale_y_log10() +
  ggtitle(label = "Frequency Increment Test (s)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "FIT_qvals.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- pdir %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = p_directional / (1 - p_directional))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 1, col = "red") +
  scale_y_log10() +
  ggtitle(label = "P(directional)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "pdir_or.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- pdir %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = p_pos / (1 - p_pos))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 1, col = "red") +
  scale_y_log10() +
  ggtitle(label = "P(directional,+)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "ppos_or.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- pdir %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = p_neg / (1 - p_neg))) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 1, col = "red") +
  scale_y_log10() +
  ggtitle(label = "P(directional,-)") +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$statplot, "pneg_or.png")
ggsave(filename, p1, width = 5, height = 4)

#' ## Confusion matrices
#' 
#' Next we print the confusion matrices of each method. For these
#' matrices we use nominal $\alpha=0.5$ or $OR = 1$ as threshold values.
#' A limitation of confusion matrices is that they only show us the performance
#' at one threshold.

#+ confusion tables
Conf <- NULL
cat("Getting confussion tables...\n")

#
cat("\tallele frequency change pval...\n")
tab <- maf %>%
  left_join(info, by = "site_id")
tab <- caret::confusionMatrix(data = factor(tab$pval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "Selection coefficient (s) p-values") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'),
        plot.title = element_text(hjust = 1.2))
filename <- file.path(args$confdir, "maf_pval_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "maf_pval"))

#
cat("\tallele frequency change qval...\n")
tab <- maf %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr'))
tab <- caret::confusionMatrix(data = factor(tab$qval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "Selection coefficient (s) q-values") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'),
        plot.title = element_text(hjust = 1.2))
filename <- file.path(args$confdir, "maf_qval_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "maf_qval"))

#
cat("\ts_coef pval...\n")
tab <- s_coef %>%
  left_join(info, by = "site_id")
tab <- caret::confusionMatrix(data = factor(tab$pval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "Selection coefficient (s) p-values") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'),
        plot.title = element_text(hjust = 1.2))
filename <- file.path(args$confdir, "s_coef_pval_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "s_coef_pval"))

#
cat("\ts_coef qval...\n")
tab <- s_coef %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr'))
tab <- caret::confusionMatrix(data = factor(tab$qval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "Selection coefficient (s) q-values") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'),
        plot.title = element_text(hjust = 1.2))
filename <- file.path(args$confdir, "s_coef_qval_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "s_coef_qval"))

#
cat("\tFIT pval...\n")
tab <- FIT %>%
  left_join(info, by = "site_id")
tab <- caret::confusionMatrix(data = factor(tab$pval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "FIT p-values") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'))
filename <- file.path(args$confdir, "FIT_pval_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "FIT_pval"))

#
cat("\tFIT qval...\n")
tab <- FIT %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr')) 
tab <- caret::confusionMatrix(data = factor(tab$qval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "FIT q-values") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'))
filename <- file.path(args$confdir, "FIT_qval_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "FIT_qval"))

#
cat("\tP(directional)...\n")
tab <- pdir %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = p_directional > 0.5, selected)
tab <- caret::confusionMatrix(data = factor(tab$significant),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "P(directional)") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'))
filename <- file.path(args$confdir, "pdir_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "pdir"))

#
cat("\tP(directional,+)...\n")
tab <- pdir %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = p_pos > 0.5, selected)
tab <- caret::confusionMatrix(data = factor(tab$significant),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "P(directional,+)") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'))
filename <- file.path(args$confdir, "ppos_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "ppos"))

#
cat("\tP(directional,-)...\n")
tab <- pdir %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = p_neg > 0.5, selected)
tab <- caret::confusionMatrix(data = factor(tab$significant),
                              reference = factor(tab$selected),
                              positive = "TRUE")
p1 <- tab$table %>%
  as_tibble() %>%
  ggplot(aes(x = Reference, y = Prediction)) +
  geom_point(aes(size = n), shape = 19, col = 'steelblue') +
  geom_text(aes(label = n)) +
  geom_vline(xintercept = 1.5) +
  geom_hline(yintercept = 1.5) +
  ggtitle(label = "P(directional,-)") +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(color = "black", face = 'bold'))
filename <- file.path(args$confdir, "pneg_confmat.png")
ggsave(filename, p1, width = 3, height = 3)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "pneg"))

cat("Writing confussion table params...\n")
filename <- file.path(args$confdir, "confusion_table_params.tsv")
write_tsv(Conf, filename)

#+ conf comparison
# Conf %>%
#   print(n = 14)
# Conf %>% filter(term == "accuracy")
p1 <- Conf %>%
  filter(term == "pos_pred_value") %>%
  ggplot(aes(x = test, y = estimate)) +
  geom_bar(stat = "identity") +
  ggtitle(label = "Positive Predictive Value") +
  AMOR::theme_blackbox()
filename <- file.path(args$confdir, "ppv.png")
ggsave(filename, p1, width = 5, height = 4)

p1 <- Conf %>%
  filter(term == "balanced_accuracy") %>%
  ggplot(aes(x = test, y = estimate)) +
  geom_bar(stat = "identity") +
  ggtitle(label = "(Sensitivity + Specificity) / 2") +
  AMOR::theme_blackbox()
filename <- file.path(args$confdir, "balanced_accuracy.png")
ggsave(filename, p1, width = 5, height = 4)

#' ## ROC curves
#' 
#' ROC curves allow us to compare the performance of classifier methods
#' at all thresholds.
#' 
#' An important caveat is that the number of tests is not identical
#' across methods since the selection coefficient & FIT tests discard some
#' observations and sites that don't fit their inclusion criteria, while
#' the P(directional) statistics always include all sites and populations
#' simulated.

#+ roc curves
ROC <- bind_rows(maf %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = -log10(pval)) %>%
                   roc() %>%
                   mutate(test = "maf"),
                 
                 
                 maf %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = -log10(p.adjust(pval, method = 'fdr'))) %>%
                   roc() %>%
                   mutate(test = "maf_FDR"),
                 
                 s_coef %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = -log10(pval)) %>%
                   roc() %>%
                   mutate(test = "s_coef"),
                 
                 s_coef %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = -log10(p.adjust(pval, method = 'fdr'))) %>%
                   roc() %>%
                   mutate(test = "s_coef_FDR"),
                 
                 FIT %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = -log10(pval)) %>%
                   roc() %>%
                   mutate(test = "FIT"),
                 
                 FIT %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = -log10(p.adjust(pval, method = 'fdr'))) %>%
                   roc() %>%
                   mutate(test = "FIT_FDR"),
                 
                 pdir %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = p_directional / (1-p_directional)) %>%
                   roc() %>%
                   mutate(test = "P(directional)"),
                 
                 pdir %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = p_pos / (1-p_pos)) %>%
                   roc() %>%
                   mutate(test = "P(directional,+)"),
                 
                 pdir %>%
                   left_join(info, by = "site_id") %>%
                   transmute(truth = selected,
                             score = p_neg / (1-p_neg)) %>%
                   roc() %>%
                   mutate(test = "P(directional,-)"))

#' Write all ROCs
filename <- file.path(args$rocdir, "roc_curves.tsv")
write_tsv(ROC, filename)


#' Plot all ROCs
p1 <- ROC %>%
  ggplot(aes(x = fpr, y = tpr, group = test))+
  geom_line(aes(col = test)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$rocdir, "roc_curves.png")
ggsave(filename, p1, width = 6, height = 4)

#' We plot a subset of ROCs that are likely to be included in final paper
selected_tests <- c("maf_FRD", "s_coef_FDR", "FIT_FDR", "P(directional)")
ROC %>%
  filter(test %in% selected_tests)



#' When lines of different ROC curves intersect it can be hard to make a decision
#' regarding a method. A standard way to summarise the information of ROC
#' curves is by calculating the Area Under the Curve (AUC). EWe use the trapezoid
#' rule.

#+ r ROC AUC
AUC <- ROC %>%
  split(.$test) %>%
  map_dfr(function(d){
    tibble(AUC = integrate(f = approxfun(x = d$fpr, y = d$tpr,
                                         ties = min, yleft = 0, yright = 1),
              lower = 0, upper = 1,
              abs.tol = 0.01,
              subdivisions = 1000)$value)
  }, .id = "test")
filename <- file.path(args$rocdir, "roc_auc.tsv")
write_tsv(AUC, filename)

#+ ROC AUC table, results='asis'
knitr::kable(AUC, caption = "AUC of ROC curves")
