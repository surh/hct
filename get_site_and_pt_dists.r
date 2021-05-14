#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
#
# This file is part of HMVAR.
#
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Read site level data and produce number of",
                        "maf increase and decrease from beginning and",
                        "end per site and patient"))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Site-level data, each row must be a site",
                                 "in a specific patient. Already QC'd"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--prop_thres",
                     help = paste("Threshold of proportion of patients",
                                  "with site, and sites per patient."),
                     type = "numeric",
                     default = 0.8)
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
# args <- list(input = "site_data.tsv.gz",
#              outdir = "output",
#              prop_thres = 0.8)

# Other dependencies
library(tidyverse)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Read data
cat("Reading site data...\n")
Dat <- read_tsv(args$input,
                col_types = cols(site_id = col_character(),
                                 sample = col_character(),
                                 depth = col_number(),
                                 freq = col_number(),
                                 pt = col_character(),
                                 date = col_date(format = "%Y-%m-%d"),
                                 day = col_number(),
                                 ref_id = col_character(),
                                 ref_pos = col_number(),
                                 gene_id = col_character(),
                                 recode = col_logical()))

# For every site-patient combination, calculate maf change
cat("Calculating MAF change per site-patient combination...\n")
Dat <- Dat %>%
  # head(10000) %>%
  group_by(site_id, pt) %>%
  summarise(change = freq[ day == max(day) ] - freq[ day == min(day) ],
            depth_start = depth[ day == min(day) ],
            depth_end = depth[ day == max(day) ],
            n_days = max(day),
            # recode = unique(recode),
            .groups = 'drop')

cat("Plotting time vs maf change...\n")
p1 <- Dat %>%
  ggplot(aes(x = n_days, y = change)) +
  geom_hex(aes(fill = ..count.. + 1)) +
  scale_fill_gradientn(colors = c("#fff7ec", "#fdd49e", "#fc8d59",
                                  "#d7301f", "#7f0000"),
                       trans = "log10") +
  theme_classic()
filename <- file.path(args$outdir, "hexbin_days_vs_change.jpeg")
ggsave(filename, p1, width = 6, height = 4, dpi = 150)

# Get distributions per site and patient & filter by prop patientes
cat("Calculating sites distribution...\n")
Sites <- Dat %>%
  group_by(site_id) %>%
  summarise(mean_depth_start = mean(depth_start),
            mean_depth_end = mean(depth_end),
            n_decrease = sum(change < 0),
            n_equal = sum(change == 0),
            n_increase = sum(change > 0),
            n_patients = length(pt),
            .groups = 'drop') %>%
  filter(n_patients >= args$prop_thres * max(n_patients))

# Make binomial test on sites
cat("\tBinomial test...\n")
f_increase <- Sites$n_increase / Sites$n_patients
f_decrease <- Sites$n_decrease / Sites$n_patients
f_change <- (Sites$n_increase + Sites$n_decrease) / Sites$n_patients
Sites <- Sites %>%
  mutate(pval_inc = 1 - (pbinom(q = n_increase - 1, size = n_patients, prob = mean(f_increase))),
         pval_dec = 1 - (pbinom(q = n_decrease - 1, size = n_patients, prob = mean(f_decrease))))
cat("Writing site distribution...\n")
filename <- file.path(args$outdir, "sites_dist.tsv.gz")
write_tsv(Sites, filename)

# Some plots for sites
cat("Comparing p-values of increase vs decrease in MAF...\n")
p1 <- Sites %>%
  ggplot(aes(x = -log10(pval_inc), y = -log10(pval_dec))) +
  # geom_point() + 
  geom_hex(aes(fill = log10(..count..))) +
  # geom_smooth(method = "lm") +
  theme_classic()
filename <- file.path(args$outdir, "hexbin_binomial_inc_vs_dec.jpeg")
ggsave(filename, p1, width = 6, height = 4, dpi = 150)

cat("Relationship between sequencing depth and change enrichment...\n")
p1 <- Sites %>%
  transmute(depth = (mean_depth_end - mean_depth_start ), pval_inc, pval_dec) %>%
  pivot_longer(-depth, names_to = "test", values_to = "pval") %>%
  ggplot(aes(x = depth, y = -log10(pval))) +
  facet_wrap(~ test, scales = "free") +
  geom_point() + 
  geom_smooth(method = "lm") +
  xlab(label = "Change in mean depth over time") +
  theme_classic()
filename <- file.path(args$outdir, "depthchange_vs_binompval.jpeg")
ggsave(filename, p1, width = 6, height = 4, dpi = 150)


# Sites %>%
#   filter(pval_inc < 0.05 & pval_dec < 0.05) %>%
#   arrange(desc(n_increase + n_decrease))

# Sites %>%
#   arrange(desc(prank_increase)) %>%
#   print(n = 200)
# Sites %>%
#   arrange(desc(prank_decrease)) %>%
#   print(n = 200)

# Patient level distribution
cat("Calculating patient distributions...\n")
Pts <- Dat %>%
  group_by(pt) %>%
  summarise(n_decrease = sum(change < 0),
            n_equal = sum(change == 0),
            n_increase = sum(change > 0),
            n_sites = length(site_id),
            .groups = 'drop') %>%
  filter(n_sites >= args$prop_thres * max(n_sites)) %>%
  mutate(p_increase = n_increase / n_sites,
         p_decrease = n_decrease / n_sites,
         p_change = (n_increase + n_decrease) / n_sites)
# Pts %>%
#    arrange(desc(prank_change)) %>%
#    print(n = 50)
# boxplot(Pts$p_change)
# qqnorm(scale(Pts$p_change))
# abline(a = 0, b = 1)
cat("Writing patient level distribution...\n")
filename <- file.path(args$outdir, "pts_dist.tsv.gz")
write_tsv(Pts, filename)
