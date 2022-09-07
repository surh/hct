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

# setwd("/cashew/users/sur/exp/fraserv/2022/today3")


library(tidyverse)
source(file.path(this.path::this.dir(), "sim_functions.r"))
args <- list(simdir = "sims/")





d <- list.dirs(args$simdir, recursive = FALSE)[1]

maf_changes <- list.dirs(file.path(d, "freqs"), recursive = FALSE) %>%
  map_dfr(function(popdir){
    cat(popdir, "\n")
    pop <- basename(popdir)
    freq <- HMVAR::read_midas_abun(file.path(popdir, "snp_freqs.txt"))
    tpoints <- setdiff(colnames(freq), "site_id")
    tpoints <- str_remove(tpoints, "^t") %>% as.numeric()
    get_site_maf_changes(freq = freq, t_0 = min(tpoints), t_n = max(tpoints)) %>%
      mutate(pop_id = as.character(pop))
  })
HMVAR::s_coefficient(Dat = maf_changes, aggregate = TRUE, time = 500)













