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

# setwd("~/micropopgen/exp/2022/today/")
library(tidyverse)
source(file.path(this.path::this.dir(), "functions.r"))
args <- list(comps_dir = "old_moran_sim_selmethods/2022-05-17.comparisons/",
             meta = "pars_per_simulation_moran.tsv",
             # nsites = "sites_per_sim.tsv",
             outdir = "auc_output/")

#' Read metadata of slimulations
Meta <- read_tsv(args$meta,
                 col_types = cols(sim_id = col_character(),
                                  .default = col_number()))
#' Number of sites under selection per simulation
Meta <- Meta %>%
  mutate(n_selected = nsites * prop_selected) 
Meta

#' Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

#' # AUC summaries
#' Read all AUC files for when there was selection
AUC <- list.dirs(args$comps_dir, recursive = FALSE, full.names = TRUE) %>%
  map_dfr(function(slim_dir){
    sim_id <- basename(slim_dir)
    cat(sim_id, "\n")
    
    roc_dir <- file.path(slim_dir, "roc")
    if(dir.exists(roc_dir)){
      auc_file <- file.path(roc_dir, "roc_auc.tsv")
      AUC <- read_tsv(auc_file, col_types = cols(test = col_character(),
                                                 AUC = col_number())) %>%
        mutate(sim_id = sim_id)
    }else{
      AUC <- NULL
    }
    
    return(AUC)
  })
AUC

#' Now we plot the AUCs
#' # NOTE: This uses a slightly different graphical representation than
#' the SLiMulation. That is because the beeswarm looks really ugly with so many
#' overlapping points near the top
table(AUC$test) %>%
  as.data.frame() %>% as_tibble() %>%
  rename(test = Var1) %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) 

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  ggplot(aes(x = test, y = AUC)) +
  # geom_boxplot(outlier.colour = NA) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.5)) +
  # ggbeeswarm::geom_beeswarm(priority = "none") +
  geom_point(position = position_jitter(width = 0.3), alpha = 0.2) +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_overall_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_overall_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' #STILL NEED TO FIX ROC CALCULATIONS 
#' Value below should be zero
1 - max(AUC$AUC)
warning("FIX ROC")

#' We see that P(directional) is on par with s_coef and both are better than
#' other methods. Though in some cases, all perform well.



#' ## Stats overall AUC
#' Perform simple test to determine significance
#' Need to account for paired structure
dat <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) 
m1 <- lmerTest::lmer(AUC ~ test + (1 | sim_id), data = dat)
m2 <- lm(AUC ~ test, data = dat)

#' We use the Likelihood ratio test to confirm that the random effect is
#' significant. Two tests agree
anova(m1, m2, "LRT")
summary(m1)

lmerTest::ranova(m1)


#' Finally we make all pairwise comparisons
lmerTest::difflsmeans(m1)
p.adjust(lmerTest::difflsmeans(m1)$`Pr(>|t|)`, 'fdr')

#' No significant difference between s_cief & P(directional).
#' Other methods are worse with MAF in last place

#' # AUC over params
Meta
#' Compare vs npops

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta,
            by = "sim_id") %>%
  ggplot(aes(x = test, y = AUC)) +
  facet_wrap(~ factor(npops), nrow = 1) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  # ggbeeswarm::geom_beeswarm(aes(col = test), size = 1, cex = 2) +
  geom_point(aes(col = test), size = 1.5, cex = 2, position = position_jitter(width = 0.4),
             alpha = 0.4) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(label = "npops") +
  AMOR::theme_blackbox() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
p1
filename <- file.path(args$outdir, "auc_by_npops_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_npops_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' We see that as the number of populations increase, also the
#' ROC AUC increases. 


#' Now we compare by popsize
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta,
            by = "sim_id") %>%
  ggplot(aes(x = test, y = AUC)) +
  facet_wrap(~ factor(popsize), nrow = 1) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  # ggbeeswarm::geom_beeswarm(aes(col = test), size = 1, cex = 2) +
  geom_point(aes(col = test), size = 1.5,
             cex = 2,
             position = position_jitter(width = 0.4),
             alpha = 0.4) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(label = "popsize") +
  AMOR::theme_blackbox() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
p1
filename <- file.path(args$outdir, "auc_by_popsize_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_popsize_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' Interestingly population size (i.e. number of genomes simulated
#' within each population) has no apparent effect.

#' Now we plot by number of sites
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta,
            by = "sim_id") %>%
  ggplot(aes(x = test, y = AUC)) +
  facet_wrap(~ factor(nsites), nrow = 1) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  # ggbeeswarm::geom_beeswarm(aes(col = test), size = 1, cex = 2) +
  geom_point(aes(col = test), cex = 2, 
             size = 1.5,
             position = position_jitter(width = 0.4),
             alpha = 0.4) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(label = "nsites") +
  AMOR::theme_blackbox() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
p1
filename <- file.path(args$outdir, "auc_by_nsites_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_nsites_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' Limited effect of number of sites. Probably not unexpected given that
#' each site is simulated independently by the Moran process.

#' Plot by elative strength of selection (g)
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta,
            by = "sim_id") %>%
  ggplot(aes(x = test, y = AUC)) +
  facet_wrap(~ factor(g), nrow = 1) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  # ggbeeswarm::geom_beeswarm(aes(col = test), size = 1, cex = 2) +
  geom_point(aes(col = test), cex = 2, 
             size = 1.5,
             position = position_jitter(width = 0.4),
             alpha = 0.4) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(label = "Relative selection (g)") +
  AMOR::theme_blackbox() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
p1
filename <- file.path(args$outdir, "auc_by_sel_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_sel_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' This has the largest impact, as expected. Although
#' FIT behaves in unexpected manner.


#' Plot by normalized time ny po size.
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta,
            by = "sim_id") %>%
  ggplot(aes(x = test, y = AUC)) +
  facet_wrap(~ factor(T / popsize), nrow = 1) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  # ggbeeswarm::geom_beeswarm(aes(col = test), size = 1, cex = 2) +
  geom_point(aes(col = test), cex = 2, 
             size = 1.5,
             position = position_jitter(width = 0.4),
             alpha = 0.4) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(label = "T / popsize") +
  AMOR::theme_blackbox() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
p1
filename <- file.path(args$outdir, "auc_by_time_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_time_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' We see that the more "Time" has passed, more effective most methods
#' Not perfect but at least improvement seems monotonic. Some variation
#' is normal. It reflects more time for selection to act.
#' 
#' T is the number of individual replacements attempted in the simulation. So
#' T/popsize is the number of replacements as a function of the population
#' size. Note that it is **NOT** the fraction of population replaced, nor
#' the expected value of that fraction, because individuals can survive
#' replacement attempt, and can be chosen again, can be replaced even if
#' they are not founded, and selection affects the probs of replacement.



#' Finally plot by number of selected sites
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta,
            by = "sim_id") %>%
  ggplot(aes(x = test, y = AUC)) +
  facet_wrap(~ factor(n_selected), nrow = 1) +
  geom_point(aes(col = test), cex = 2, 
             size = 1.5,
             position = position_jitter(width = 0.4),
             alpha = 0.4) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  # scale_x_log10() +
  ggtitle(label = "n_selected") +
  AMOR::theme_blackbox() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
p1
filename <- file.path(args$outdir, "auc_by_nselected_moran.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_nselected_moran.svg")
ggsave(filename, p1, width = 6, height = 4)

#' No differences under this variable. Expected given that all
#' sites are simulated independently
