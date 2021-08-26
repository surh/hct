#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


library(tidyverse)
library(HMVAR)

# args <- list(dndscv = "",
#              cds = "~/micropopgen/exp/2021/2021-08-20.dnds/refdb/Fplautii.cdsfile.tsv",
#              pdir = "~/micropopgen/exp/2021/2021-08-20.dnds/p_directional/MGYG-HGUT-00099.tsv.gz",
#              info = "~/micropopgen/exp/2021/2021-08-20.dnds/Fplautii_merged/snps_info.txt",
#              outdir = "output")
args <- list(dndscv = "all_dummy/dnds_cv.tsv",
             cds = "test_refdb/cdsfile.tsv",
             pdir = "test_in/MGYG-HGUT-00099.tsv.gz",
             info = "test_in/MGYG-HGUT-00099/snps_info.txt",
             outdir = "test_plot")

cat("Reading map of sites to genes...\n")
site2gene <- HMVAR::read_midas_info(args$info) %>%
  select(site_id, gene_id)

cat("Reading map of gene names to gene IDs...\n")
name2.id <- read_tsv(args$cds,
         col_types = cols(gene.id = col_character(),
                          gene.name = col_character())) %>%
  select(gene_id = gene.id, gene_name = gene.name)

cat("Read p_directional...\n")
pdir <- read_tsv(args$pdir, 
                 col_types = cols(site_id = col_character()))

cat("Get gene-level p_directional...\n")
pdir <- pdir %>%
  select(site_id, p_directional) %>%
  left_join(site2gene, by = "site_id") %>%
  group_by(gene_id) %>%
  summarise(p_directional = max(p_directional),
            .groups = 'drop')
  
cat("Reading dndscv results...\n")
Res <- read_tsv(args$dndscv,
                col_types = cols(gene_name = col_character(),
                                 .default = col_number))

cat("Matching dndscv and p_directional...\n")
Res <- Res %>%
  select(gene_name, n_syn, n_mis, n_non,
         wmis_cv, wnon_cv,
         pmis_cv, ptrunc_cv, pallsubs_cv,
         qmis_cv, qtrunc_cv, qallsubs_cv) %>%
  left_join(name2id, by = "gene_name") %>%
  inner_join(pdir, by = "gene_id")

# Prepare output dir
cat("Creating output directory...\n")
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

cat("Plotting...\n")
p1 <- Res %>%
  mutate(pmis_cv = -log10(pmis_cv),
         ptrunc_cv = -log10(ptrunc_cv),
         pallsubs_cv = -log10(pallsubs_cv),
         qmis_cv = -log10(qmis_cv),
         qtrunc_cv = -log10(qtrunc_cv),
         qallsubs_cv = -log10(qallsubs_cv)) %>%
  mutate(pmis_cv = replace(pmis_cv, pmis_cv == Inf, min(pmis_cv[pmis_cv != Inf])),
         ptrunc_cv = replace(ptrunc_cv, ptrunc_cv == Inf, min(ptrunc_cv[ptrunc_cv != Inf])),
         pallsubs_cv = replace(pallsubs_cv, pallsubs_cv == Inf, min(pallsubs_cv[pallsubs_cv != Inf])),
         qmis_cv = replace(qmis_cv, qmis_cv == Inf, min(qmis_cv[qmis_cv != Inf])),
         qtrunc_cv = replace(qtrunc_cv, qtrunc_cv == Inf, min(qtrunc_cv[qtrunc_cv != Inf])),
         qallsubs_cv = replace(qallsubs_cv, qallsubs_cv == Inf, min(qallsubs_cv[qallsubs_cv != Inf]))) %>%
  mutate(wmis_cv = log2(wmis_cv),
         wnon_cv = log2(wnon_cv)) %>%
  pivot_longer(cols = c(-gene_name, -gene_id, -p_directional),
               names_to = "statistic", values_to = "value") %>%
  filter(!(statistic == "wmis_cv" & value == -Inf)) %>%
  filter(!(statistic == "wnon_cv" & value == -Inf)) %>%
  
  ggplot(aes(x = log2(p_directional / (1 - p_directional)), y = value)) +
  facet_wrap(~ statistic, scales = "free_y", ncol = 3) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()
filename <- file.path(args$outdir, "dndscv_vs_pdirectional.loglog.png")
ggsave(filename, p1,
       width = 6, height = 6, dpi = 150)

# Res %>%
#   ggplot(aes(x = p_directional >= 0.8, y = wmis_cv)) +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_point(position = position_jitter(width = 0.1)) +
#   scale_y_log10() +
#   AMOR::theme_blackbox()
# 
# Res %>%
#   filter(wnon_cv > 0) %>%
#   ggplot(aes(x = p_directional >= 0.8, y = wnon_cv)) +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_point(position = position_jitter(width = 0.1)) +
#   scale_y_log10() +
#   AMOR::theme_blackbox()
