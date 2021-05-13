#!/usr/bin/env Rscript
library(tidyverse)

args <- list(input = "site_data.tsv.gz",
             outdir = "output",
             prop_thres = 0.8)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Read data
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
Dat <- Dat %>%
  # head(10000) %>%
  group_by(site_id, pt) %>%
  summarise(change = freq[ day == max(day) ] - freq[ day == min(day) ],
            depth_start = depth[ day == min(day) ],
            depth_end = depth[ day == max(day) ],
            n_days = max(day),
            # recode = unique(recode),
            .groups = 'drop')

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
f_increase <- Sites$n_increase / Sites$n_patients
f_decrease <- Sites$n_decrease / Sites$n_patients
f_change <- (Sites$n_increase + Sites$n_decrease) / Sites$n_patients
Sites <- Sites %>%
  mutate(pval_inc = 1 - (pbinom(q = n_increase - 1, size = n_patients, prob = mean(f_increase))),
         pval_dec = 1 - (pbinom(q = n_decrease - 1, size = n_patients, prob = mean(f_decrease))))
filename <- file.path(args$outdir, "sites_dist.tsv.gz")
write_tsv(Sites, filename)

# Some plots for sites
p1 <- Sites %>%
  ggplot(aes(x = -log10(pval_inc), y = -log10(pval_dec))) +
  # geom_point() + 
  geom_hex(aes(fill = log10(..count..))) +
  # geom_smooth(method = "lm") +
  theme_classic()
filename <- file.path(args$outdir, "hexbin_binomial_inc_vs_dec.jpeg")
ggsave(filename, p1, width = 6, height = 4, dpi = 150)

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
filename <- file.path(args$outdir, "pts_dist.tsv.gz")
write_tsv(Pts, filename)
