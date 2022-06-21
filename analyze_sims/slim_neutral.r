#!/usr/bin/env Rscript

# (C) Copyright 2022 Sur Herrera Paredes
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
# along with This program.  If not, see <https://www.gnu.org/licenses/>.

#' # Setup
library(tidyverse)

args <- list(comps_dir = "2022-06-16.comparisons/comparisons/",
             estimates_dir = "2022-06-03.selection_methods/",
             meta = "pars_per_simulation.tsv",
             nsites = "sites_per_sim.tsv",
             outdir = "neutral_output/")

#' We choose the methods and colors to display
# selected_tests <- c("maf_FDR", "s_coef_FDR", "FIT_FDR",
#                     "P(directional)", "P(directional,-)")
# test_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")
selected_tests <- c("maf_FDR", "s_coef_FDR", "FIT_FDR",
                    "P(directional)")
test_colors <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")

#' Read metadata of slimulations
Meta <- read_tsv(args$meta,
                 col_types = cols(sim_id = col_character(),
                                  run_id_short = col_character(),
                                  .default = col_number()))
Meta

Nsites <-read_tsv(args$nsites,
                  col_types = cols(sim_id = col_character(),
                                   
                                   
                                   .default = col_number()))
Nsites
Meta <- Meta %>%
  left_join(Nsites %>%
              select(sim_id, n_sites, n_standing, n_selected),
            by = "sim_id")
Meta

#' Select slimulations without sites under selection
Meta <- Meta %>%
  filter(n_selected == 0)
Meta


# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

#' The code below calculates number of false positive as desired thresholds
#' It is for demonstration purposes, need to be incorporated at compare script:w
FPRs <- list.dirs(args$comps_dir, recursive = FALSE, full.names = TRUE) %>%
  map_dfr(function(slim_dir){
    sim_id <- basename(slim_dir)
    
    subdir <- file.path(slim_dir, "summaries")
    filename <- file.path(subdir, "n_fps_neutral.tsv")
    if(file.exists(filename)){
      dat <- read_tsv(filename, col_types = cols(test = col_character(),
                                                 .default = col_number())) %>%
        mutate(sim_id = sim_id)
    }else{
      dat <- NULL
    }
    
    return(dat)
  })
FPRs


#' Whe check that we get the right number of fpr values per threshold and sim
FPRs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>% ftable(test ~ thres, .)

#' After reading the false positive rates we plot the overall
#' results

p1 <- FPRs %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  ggplot(aes(x = thres, y = fpr)) +
  facet_wrap(~ test, scales = "free_x") +
  # ggbeeswarm::geom_beeswarm(aes(fill = thres), size = 2, col = "black", shape = 21, cex = 0.8) +
  ggbeeswarm::geom_beeswarm(aes(col = thres), size = 2, cex = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("#bae4bc", "#7bccc4", "#43a2ca", "#0868ac",
                                "#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  geom_hline(yintercept = 0.05, col = "black") +
  ylim(c(0, 1)) +
  AMOR::theme_blackbox() 
p1
filename <- file.path(args$outdir, "fpr_neutral.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "fpr_neutral.svg")
ggsave(filename, p1, width = 6, height = 4)

#' The FPR for maf & FIT is really hight. For s_coef it is very low, but
#' this comes at cost of extreme lost of power according to what the AUC
#' of the ROC suggests.
#' 
#' For P(directional) FDR is around the classic 0.05 threshold at an OR
#' of 2 and it is not much higher for OR1, while clearly below this
#' threshold for OR 3 and 4.
#' 
#' IMPORTANT: note that the p-value-like thresholds cannot be compared to
#' the odds rations, they are measuring different things and thus they
#' belong in different categories.


#' We quickly focus on P(directional) to see if specific slimulation
#' attributes have any effect

p1 <- FPRs %>%
  filter(test %in% "P(directional)") %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(1, 2, 3, 4))) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen),
            by = "sim_id") %>%
  ggplot(aes(x = factor(2 * Ne * Mu), y = fpr)) +
  # facet_wrap(~ thres) +
  ggbeeswarm::geom_beeswarm(aes(col = thres),dodge.width = 1, cex = 2) +
  scale_color_manual(values = c("#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "fpr_neutral_pdir_theta.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "fpr_neutral_pdir_theta.svg")
ggsave(filename, p1, width = 6, height = 4)

#' the overall amount of diversity doesn't affect the FPR for P(directional).
#' Note that $Ne$ is constant so $\Mu$ defines theta.

p1 <- FPRs %>%
  filter(test %in% "P(directional)") %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(1, 2, 3, 4))) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen),
            by = "sim_id") %>%
  ggplot(aes(x = factor(Rho), y = fpr)) +
  # facet_wrap(~ thres) +
  ggbeeswarm::geom_beeswarm(aes(col = thres),dodge.width = 1, cex = 1.5) +
  scale_color_manual(values = c("#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "fpr_neutral_pdir_rho.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "fpr_neutral_pdir_rho.svg")
ggsave(filename, p1, width = 6, height = 4)

p1 <- FPRs %>%
  filter(test %in% "P(directional)") %>%
  filter(sim_id %in% Meta$sim_id) %>%
  mutate(thres = factor(thres, levels = c(1, 2, 3, 4))) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen),
            by = "sim_id") %>%
  ggplot(aes(x = factor(tractlen * Rho), y = fpr)) +
  # facet_wrap(~ thres) +
  ggbeeswarm::geom_beeswarm(aes(col = thres),dodge.width = 1, cex = 1.5) +
  scale_color_manual(values = c("#fdcc8a", "#fc8d59", "#e34a33", "#b30000")) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "fpr_neutral_pdir_hgt.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "fpr_neutral_pdir_hgt.svg")
ggsave(filename, p1, width = 6, height = 4)

#' There might be a small effect of Rho, but not tractlen on the FPR. Higher,
#' Rho may lead to higher FPR
#' 
#' We do a model to test this
dat <- FPRs %>%
  filter(test %in% "P(directional)") %>%
  filter(sim_id %in% Meta$sim_id) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen),
            by = "sim_id") %>%
  mutate(Mu = log10(Mu),
         Rho = log10(Rho),
         tractlen = log10(tractlen))
dat

m1 <- lmerTest::lmer(fpr ~ Mu + Rho + tractlen + (1|sim_id), data = dat)
summary(m1)
summary(m1)$coefficients

#' random effect is higly significant
lmerTest::ranova(m1)


#' to deal with singluarity use brms
m1.brms <- brms::brm(fpr ~ Mu + Rho + tractlen + (1 | sim_id),
          data = dat, chains = 4, cores = 4)
summary(m1.brms)$fixed


m2.brms <- brms::brm(fpr ~ Rho + (1 | sim_id),
                data = dat, chains = 4, cores = 4)
summary(m2.brms)$fixed


#' Next step is to count the bias in negative p_directional. 
#' Make sure it is a general phenomenon

FPs_sign <- Meta$sim_id %>%
  map_dfr(function(sim_id, sel_dir, comp_dir){
    # sim_id <- "slim_1"
    # sel_dir <- args$estimates_dir
    # comp_dir <- args$comps_dir
    
    alpha_thres <- c(0.1, 0.05, 0.01, 0.001)
    or_thres <- c(1, 2, 3, 4)
    
    
    # MAF
    MAF <- read_tsv(file.path(comp_dir, sim_id, "MAF_change", "maf_change.tsv"),
                       col_types = cols(site_id = col_character())) %>%
      filter(!is.na(pval))
    alpha_mat <- matrix(rep(alpha_thres, each = nrow(MAF)),
                        ncol = length(alpha_thres))
    m_np <- colSums(p.adjust(MAF$pval, 'fdr') < alpha_mat & MAF$t.value > 0)
    m_nn <- colSums(p.adjust(MAF$pval, 'fdr') < alpha_mat & MAF$t.value < 0)
    
    # s_coef count pos & neg
    s_coef <- read_tsv(file.path(sel_dir, "s_coef", paste0(sim_id, ".tsv")),
                       col_types = cols(site_id = col_character())) %>%
      filter(!is.na(pval))
    alpha_mat <- matrix(rep(alpha_thres, each = nrow(s_coef)),
                        ncol = length(alpha_thres))
    s_np <- colSums(p.adjust(s_coef$pval, 'fdr') < alpha_mat & s_coef$t.s > 0)
    s_nn <- colSums(p.adjust(s_coef$pval, 'fdr') < alpha_mat & s_coef$t.s < 0)
    
    # FIT
    FIT <- read_tsv(file.path(sel_dir, "FIT", paste0(sim_id, ".tsv")),
                       col_types = cols(site_id = col_character())) %>%
      filter(!is.na(pval))
    alpha_mat <- matrix(rep(alpha_thres, each = nrow(FIT)),
                        ncol = length(alpha_thres))
    f_np <- colSums(p.adjust(FIT$pval, 'fdr') < alpha_mat & FIT$z.score > 0)
    f_nn <- colSums(p.adjust(FIT$pval, 'fdr') < alpha_mat & FIT$z.score < 0)
    
    # P(directional)
    pdir <- read_tsv(file.path(sel_dir, "p_directional", paste0(sim_id, ".tsv.gz")),
                    col_types = cols(site_id = col_character())) %>%
      filter(!is.na(p_directional))
    or_mat <- matrix(rep(or_thres, each = nrow(pdir)),
                        ncol = length(or_thres))
    p_np <- colSums((pdir$p_directional / (1- pdir$p_directional)) > or_mat & (pdir$p_pos > pdir$p_neg))
    p_nn <- colSums((pdir$p_directional / (1- pdir$p_directional)) > or_mat & (pdir$p_pos < pdir$p_neg))
    
    tibble(test = rep(selected_tests, each = 4),
           thres = c(alpha_thres, alpha_thres, alpha_thres, or_thres),
           n_pos = c(m_np, s_np, f_np, p_np),
           n_neg = c(m_nn, s_nn, f_nn, p_nn),
           sim_id = sim_id)
    
  }, sel_dir = args$estimates_dir, comp_dir = args$comps_dir)
FPs_sign

#' The steps above need to be incorporated in the compare methods script of
#' the general pipeline


p1 <- FPs_sign %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  pivot_longer(-all_of(c("test", "thres", "sim_id")), names_to = "sign", values_to = "n_fp") %>%
  filter(test != "P(directional)") %>%
  ggplot(aes(x = sim_id, y = n_fp)) +
  facet_grid(thres ~ test, scales = "free_x", space = "free_x", ) +
  geom_bar(aes(fill = sign), stat = "identity", position = "fill", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_blank())
p1


p2 <- FPs_sign %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  mutate(thres = factor(thres, levels = c(0.1, 0.05, 0.01, 0.001, 1, 2, 3, 4))) %>%
  pivot_longer(-all_of(c("test", "thres", "sim_id")), names_to = "sign", values_to = "n_fp") %>%
  filter(test == "P(directional)") %>%
  ggplot(aes(x = sim_id, y = n_fp)) +
  facet_grid(thres ~ test, scales = "free_x", space = "free_x", ) +
  geom_bar(aes(fill = sign), stat = "identity", position = "fill", width = 1) +
  theme_classic() +
  theme(axis.text.x = element_blank())
p2

p_legend <- cowplot::get_legend(p2)
p_legend

p_graph <- cowplot::plot_grid(p1 + theme(legend.position = "none"),
                         p2 + theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    legend.position = "none"),
                         nrow = 1,
                         rel_widths = c(0.74, 0.26))
p_graph

pp <- cowplot::plot_grid(p_legend, p_graph,
                        nrow = 2,
                        rel_heights = c(0.15, 0.85))
pp

filename <- file.path(args$outdir, "fpr_neutral_sign.png")
ggsave(filename, pp, width = 6, height = 4)
filename <- file.path(args$outdir, "fpr_neutral_sign.svg")
ggsave(filename, pp, width = 6, height = 4)

#' We see that ofr P(directional) there is a big bias in the
#' fase positives towards negative directionality. This
#' probably means drift is detected at some level as selection.
#' One possibility is to focus only on alleles that start as minor alleles
#' in the real data

#' # Session Info
sessionInfo()

