#!/usr/bin/env Rscript

# (C) Copyright 2022 Sur Herrera Paredes
# 
# This file is part of HCT.
# 
# HCT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HCT is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with HCT. If not, see <http://www.gnu.org/licenses/>.

# Script to check model results from from stan model

# setwd("/cashew/users/sur/exp/fraserv/2022/today2/output")
args <- list(stan_model = "stan_models/MGYG-HGUT-03439.stan.rdat",
             output = "checks")

library(tidyverse)
library(rstan)
library(bayesplot)


vars <- c("P", "Q")
load(args$stan_model)


# Show key parameters and trace
print(m1.stan, pars = vars)
bayesplot::mcmc_trace(m1.stan, pars = vars)


sampler_params <- get_sampler_params(m1.stan, inc_warmup = TRUE)
sampler_params %>% bind_rows


bayesplot::mcmc_pairs(m1.stan, pars = vars)

bayesplot::mcmc_pairs(m1.stan, pars = c(vars, "lp__"))

m1.stan






# setwd("/cashew/users/sur/exp/fraserv/2022/today3/preHCT_posHCT_2/output/model_summaries/")
# setwd("/cashew/users/sur/exp/fraserv/2022/today3/preFOS_posFOS_control_2/output/model_summaries/")
setwd("/cashew/users/sur/exp/fraserv/2022/today2/output/model_summaries/")


library(tidyverse)


# Dat <- read_tsv("MGYG-HGUT-02438.tsv.gz")
# Dat <- read_tsv("MGYG-HGUT-01029.tsv.gz")
# Dat <- read_tsv("MGYG-HGUT-00144.tsv.gz")
Dat <- read_tsv("MGYG-HGUT-03439.tsv.gz")
Dat

Dat %>%
  mutate(var_type = str_remove(var, "[\\[][0-9]+[\\]]")) %>%
  filter(var_type == "p") %>%
  ggplot(aes(x = `50%`)) +
  geom_histogram(bins = 5)


Dat %>%
  mutate(var_type = str_remove(var, "[\\[][0-9]+[\\]]")) %>%
  filter(var_type == "q") %>%
  ggplot(aes(x = `50%`)) +
  geom_histogram(bins = 10)
Dat

Dat %>%
  mutate(var_type = str_remove(var, "[\\[][0-9]+[\\]]")) %>%
  group_by(var_type) %>%
  summarise(avg = mean(mean),
            min = min(mean),
            max = max(mean))
Dat


n <- 9218

m <- 0.0158
v <- 7
hist(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v), breaks = 20)
range(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v))

m <- 0.476
v <- 30
hist(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v), breaks = 20)
range(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v))

