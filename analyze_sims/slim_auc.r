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
             meta = "pars_per_simulation.tsv",
             nsites = "sites_per_sim.tsv",
             outdir = "auc_output/")

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

                  
# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

#' # AUC summaries
#' Read all AUC files for when there was selection
AUC <- list.dirs(args$comps_dir, recursive = FALSE, full.names = TRUE) %>%
  map_dfr(function(slim_dir){
    sim_id <- basename(slim_dir)
    
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
  ggbeeswarm::geom_beeswarm() +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_overall_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_overall_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' We see that P(directional clearly outperforms other methods). FIT
#' is not as bad here, and it is on par with s_coef while looking at the
#' derived allele frequency change (maf) is clearly not great


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

#' Pairwise difference are significant, indicating that P(directional)
#' outperforms all other methods, after we correct for multiple testing.

#' # AUC over params

#' Compare vs theta

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection),
            by = "sim_id") %>%
  ggplot(aes(x = factor(2 *Ne * Mu), y = AUC)) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_by_theta_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_theta_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' We see that P(directional) always performs always perform well, while for
#' s_coef, there are few instances with low genetic variation where it is on
#' par with P(directional) but not in the majoriity. For High diversity
#' scoef always overperforms

# p1 <- AUC %>%
#   filter(test %in% selected_tests) %>%
#   mutate(test = factor(test, levels = selected_tests)) %>%
#   left_join(Meta %>%
#               select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection,
#                      n_sites, n_standing),
#             by = "sim_id") %>%
#   ggplot(aes(x = cut(n_sites, breaks = c(0, 6000,7000,50000, 59000, 64000)), y = AUC)) +
#   # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
#   # scale_fill_manual(values = test_colors) +
#   ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
#   scale_color_manual(values = test_colors) +
#   # geom_point(aes(col = test), size = 2) +
#   # scale_color_manual(values = test_colors) +
#   scale_y_continuous(limits = c(0, 1)) +
#   # scale_x_log10() +
#   AMOR::theme_blackbox()
# p1

#' Below we show that the final number of sites matches (roughly) the
#' preduction from the watterson estimator (Theta)
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection,
                     n_sites, n_standing),
            by = "sim_id") %>%
  ggplot(aes(x = factor(2 * Ne * Mu), y = n_sites)) +
  # geom_point(aes(col = test)) +
  # scale_color_manual(values = test_colors) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  scale_color_manual(values = test_colors) +
  scale_y_log10() +
  # scale_x_log10() +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "nsites_vs_theta_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "nsites_vs_theta_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' Now we plot by recombination rate
p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection),
            by = "sim_id") %>%
  ggplot(aes(x = factor(Rho * tractlen), y = AUC)) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_by_hgt_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_hgt_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' No obvious effect of recombination rate. Maube a bit underperformance
#' At very high recombination foor P(directional), is it in the other
#' direction for the rest?

#' What about selection strength

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection),
            by = "sim_id") %>%
  ggplot(aes(x = factor(scoef), y = AUC)) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_by_scoef_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_scoef_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' When the selection coefficient is high enough s_coef can sometimes
#' match P(directional), but in general it underperforms. P(directional)
#' itself improves when the selection coefficient is bigger, which 
#' is what we expect

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection),
            by = "sim_id") %>%
  ggplot(aes(x = factor(prop_selection), y = AUC)) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_by_prop_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_prop_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' The s_coef methods seems to have more extreme values when prop_selection is low,
#' but sometimes it matches P(directional). P(directional) is consistently
#' good, and it differentiates itself form the other methods better
#' when there is a lot of selection.
#' 
#' What about absolute number of selected sites

p1 <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Ne, Mu, Rho, tractlen, scoef, prop_selection,
                     n_selected),
            by = "sim_id") %>%
  ggplot(aes(x = n_selected, y = AUC)) +
  # ggbeeswarm::geom_beeswarm(aes(fill = test), col = "black", shape = 21, size = 2, cex = 2) +
  # scale_fill_manual(values = test_colors) +
  # ggbeeswarm::geom_beeswarm(aes(col = test), size = 2, cex = 2) +
  # scale_color_manual(values = test_colors) +
  geom_point(aes(col = test), size = 2) +
  scale_color_manual(values = test_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_log10(breaks = c(1, 10, 100, 500)) +
  AMOR::theme_blackbox()
p1
filename <- file.path(args$outdir, "auc_by_nselected_slim.png")
ggsave(filename, p1, width = 6, height = 4)
filename <- file.path(args$outdir, "auc_by_nselected_slim.svg")
ggsave(filename, p1, width = 6, height = 4)

#' When there is only one site under selection, s_coef can sometimes catch
#' it, probably when that one variant manages to increase, but as things
#' get more complicated, the advantage of p_directional increases.
#' It is not clear if P(directional is affected by n_selected)

#' ## Stats across params

#' We now fit a full model that tries to capture
#' how all the different variables influence
#' ROC AUC performance
dat <- AUC %>%
  filter(test %in% selected_tests) %>%
  mutate(test = factor(test, levels = selected_tests)) %>%
  left_join(Meta %>%
              select(sim_id, Mu, Rho, tractlen, scoef, prop_selection,
                     n_selected, n_sites),
            by = "sim_id") %>%
  mutate(scoef = log10(scoef),
         tractlen = log10(tractlen),
         Rho = log10(Rho),
         Mu = log10(Mu),
         n_selected = log10(n_selected))
dat
# m1 <- lmerTest::lmer(AUC ~ test + n_selected + log10(scoef) + log10(tractlen) + log10(Rho) + log10(Mu) + (1 | sim_id), data = dat)
m1 <- lmerTest::lmer(AUC ~ test + n_selected + scoef + tractlen + Rho + Mu + (1 | sim_id), data = dat)


#' The random effect controls for *pairedness* and all the other
#' variables are continuous, we log transform to keep all variables
#' in a similar scale

#' There is perhaps a small tendency to over fit towards the
#' P(directional) results, but it is mild
plot(m1)
boxplot(residuals(m1, "pearson") ~ dat$test)

#' We confirm that the random effect is significant
lmerTest::ranova(m1)

#' Then we look at the overall results
summary(m1)

#' We see that each method is better than MAF, but importabtly, only
#' the selection coefficient impacts performance overall
lmerTest::difflsmeans(m1, test.effs = "test")

#' We try a reduced model 
drop1(m1)
m2 <- lmerTest::lmer(AUC ~ test + n_selected + scoef + Rho + Mu + (1 | sim_id), data = dat)
summary(m2)
drop1(m2)
m2 <- lmerTest::lmer(AUC ~ test + n_selected + scoef + Mu + (1 | sim_id), data = dat)
summary(m2)
drop1(m2)
m2 <- lmerTest::lmer(AUC ~ test + scoef + Mu + (1 | sim_id), data = dat)
summary(m2)
drop1(m2)
m2 <- lmerTest::lmer(AUC ~ test + scoef + (1 | sim_id), data = dat)
summary(m2)
lmerTest::ranova(m2)
lmerTest::difflsmeans(m2, test.effs = "test")

#' Overall, it is clear that P(directional) performs better than the
#' other methods

#' # Session Info
sessionInfo()
