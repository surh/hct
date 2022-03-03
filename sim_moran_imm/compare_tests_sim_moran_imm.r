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



args <- list(s_coef = "/home/sur/micropopgen/exp/2022/today2/comparison_results/s_coef/sim_1.tsv",
             FIT = "/home/sur/micropopgen/exp/2022/today2/comparison_results/FIT/sim_1.tsv",
             pdir = "/home/sur/micropopgen/exp/2022/today2/comparison_results/p_directional/sim_1.tsv.gz",
             info = "/home/sur/micropopgen/exp/2022/today2/sims/sim_1/snps_info.txt",
             outdir = "comp_test/")
library(tidyverse)

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
info <- read_tsv(args$info,
                 col_types = cols(site_id = col_character()))


#' We create a directory to store output
#+ outdir
cat("Prepare output dir...\n")
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
  
  args$pvaldir <- file.path(args$outdir, "pvals")
  dir.create(args$pvaldir)
  
  args$statplot <- file.path(args$outdir, "statplot")
  dir.create(args$statplot)
  
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
p1 <- s_coef %>%
  left_join(info, by = "site_id") %>%
  ggplot(aes(x = selected, y = pval)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2) +
  geom_point(size = 0.1, position = position_jitter(width = 0.2)) +
  geom_hline(yintercept = 0.05, col = "red") +
  scale_y_log10() +
  ggtitle(label = "Selection coefficient (s)") +
  AMOR::theme_blackbox()
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
filename <- file.path(args$statplot, "pneg_or.png")
ggsave(filename, p1, width = 5, height = 4)

#' ## Confusion matrices
#' 
#' Next we print the confusion matrices of each method. For these
#' matrices we use nominal $\alpha=0.5$ or $OR = 1$ as threshold values.
#' A limitation of confusion matrices is that they only show us the performance
#' at one threshold.

#+ confusion tables
cat("Getting confussion tables...\n")
# s_coef %>%
#   left_join(info, by = "site_id") %>%
#   transmute(significant = pval < 0.05, selected) %>%
#   ftable 

tab <- s_coef %>%
  left_join(info, by = "site_id")
tab <- caret::confusionMatrix(data = factor(tab$pval < 0.05),
                              reference = factor(tab$selected),
                              positive = "TRUE")
print(tab$table)
Conf <- Conf %>% 
  bind_rows(tab %>%
              broom::tidy() %>%
              select(-class) %>%
              mutate(test = "s_coef_pval"))

s_coef %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr')) %>%
  transmute(significant = qval < 0.05, selected) %>%
  ftable

FIT %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = pval < 0.05, selected) %>%
  ftable

FIT %>%
  left_join(info, by = "site_id") %>%
  mutate(qval = p.adjust(pval, method = 'fdr')) %>%
  transmute(significant = qval < 0.05, selected) %>%
  ftable

pdir %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = p_directional > 0.5, selected) %>%
  ftable

pdir %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = p_pos > 0.5, selected) %>%
  ftable

pdir %>%
  left_join(info, by = "site_id") %>%
  transmute(significant = p_neg > 0.5, selected) %>%
  ftable



## ROC curves

ROC curves allow us to compare the performance of classifier methods
at all thresholds.

An important caveat is that the number of tests is not identical
across methods since the selection coefficient & FIT tests discard some
observations and sites that don't fit their inclusion criteria, while
the P(directional) statistics always include all sites and populations
simulated.

```{r roc curves}
ROC <- bind_rows(s_coef %>%
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

p1 <- ROC %>%
  ggplot(aes(x = fpr, y = tpr, group = test))+
  geom_line(aes(col = test)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "roc_curves.png")
ggsave(filename, p1, width = 6, height = 4)

filename <- file.path(args$outdir, "roc_curves.tsv")
write_tsv(ROC, filename)
```

When lines of different ROC curves intersect it can be hard to make a decision
regarding a method. A standard way to summarise the information of ROC
curves is by calculating the Area Under the Curve (AUC). EWe use the trapezoid
rule.

```{r ROC AUC}
AUC <- ROC %>%
  split(.$test) %>%
  map_dfr(function(d){
    tibble(AUC = integrate(f = approxfun(x = d$fpr, y = d$tpr, ties = min),
              lower = 0, upper = 1,
              abs.tol = 0.01,
              subdivisions = 1000)$value)
  }, .id = "test")
filename <- file.path(args$outdir, "roc_auc.tsv")
write_tsv(AUC, filename)
```

```{r ROC AUC table, results='asis'}
knitr::kable(AUC, caption = "AUC of ROC curves")


