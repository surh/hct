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
library(phybreak)

args <- list(burnin = 2e3,
             iter = 2e3,
             thin = 5,
             alleles = "filtered_alleles/snps_alleles.txt",
             map = "meta.txt")

alleles <- read_tsv(args$alleles,
                    col_types = cols(site_id = col_character(),
                                     .default = col_character()))
alleles <- alleles %>%
  select(-site_id)
alleles

meta <- read_tsv(args$map,
                 col_types = cols(sample = col_character(),
                                  pt = col_character(),
                                  date = col_date(format = "%Y-%m-%d")))
meta <- meta %>%
  mutate(newname = paste0(pt, ".", sample))
meta

# Create data object
alleles <- tolower(t(alleles))

# Change names to include pt
newnames <- set_names(meta$newname, meta$sample)
row.names(alleles) <- newnames[ row.names(alleles) ]

sample_names <- row.names(alleles)

host_names <- set_names(meta$pt, meta$newname)
host_names <- host_names[ sample_names ]

sample_times <- set_names(meta$date, meta$newname)
sample_times <- sample_times[ sample_names ]

dat <- phybreakdata(sequences = alleles,
                    sample.times = sample_times,
                    host.names = host_names,
                    sample.names = sample_names)

Sys.time()
m1 <- phybreak(dataset = dat)
Sys.time()
m1 <- burnin_phybreak(x = m1, ncycles = args$burnin)
Sys.time()
m1 <- sample_phybreak(x = m1, nsample = args$iter, thin = args$thin)
Sys.time()

save(m1, file = "model.rdat")


# 
# ESS(m1)
# transtree(m1)
# infectorsets(m1)
# phylotree(m1)
# 
# 
# plotTrans(m1)
# plotPhylo(m1)
# plotPhyloTrans(m1)
# 
# 
# 
# str(m1$s)
# dim(m1$s$inftimes)
# # acf(t(m1$s$inftimes))
# # pairs(t(m1$s$inftimes))
# 
# plotgg_phybreak_trace(x, params = "main"){
#   x <- m1
# 
#   if(length(params) == 1 && params == "main"){
#     dat <- tibble(iter = 1:length(x$s$logLik),
#                   mu = x$s$mu,
#                   mG = x$s$mG,
#                   mS = x$s$mS,
#                   wh.s = x$s$wh.s,
#                   wh.e = x$s$wh.0,
#                   dist.e = x$s$dist.e,
#                   dist.s = x$s$dist.s,
#                   dist.m = x$s$dist.m)
#     dat %>%
#       pivot_longer(-iter, values_to = "estimate", names_to = "parameter") %>%
#       ggplot(aes(x = iter, y = estimate)) +
#       facet_wrap(~ parameter, scales = "free_y") +
#       geom_line() +
#       theme_classic()
#   }else{
#     stop("ERROR: not implemented", call. = TRUE)
#   }
#   
# }
# 
# acf(x$s$mu[seq(from = 1, to = 100, by = 5)])







