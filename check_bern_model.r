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
             output = "checks",
             max_treedepth = 10)

library(tidyverse)
library(rstan)
library(bayesplot)


vars <- c("P", "Q")
load(args$stan_model)


# Show key parameters and trace
print(m1.stan, pars = vars)
bayesplot::mcmc_trace(m1.stan, pars = vars)


sampler_params <- get_sampler_params(m1.stan, inc_warmup = FALSE)
dat <- do.call(rbind, sampler_params) %>%
  as_tibble() %>% as.list() %>%
  map_dfr(my_summary, .id = "var") %>%
  print(digits = 3)

# Check for exceeded treedepth & divergent transitions
if((dat %>% filter(var == "treedepth__"))["max"] > args$max_treedepth)
  warning("some iterations exceeded max_treedepth")

if((dat %>% filter(var == "divergent__"))["mean"] > 0)
  warning("There were divergent transitions")


# c(args[c("output", "stan_model")], list(test = 1))
# 
# sampler_params %>%
#   map(my_summary, digits = 3)
# 

par <- "energy__"
dat <- c(rstan::extract(m1.stan, pars = c(vars, "lp__")),
  list(sampler_params %>%
         map(~ .x[ , par]) %>%
         unlist
  ))


# Prepare pairs plot with lp & energy
dat <- as.array(m1.stan, pars = c(vars, "lp__")) %>%
  apply(2, function(mat){
    mat
  }, simplify = FALSE) 
 
# str(sampler_params)
# str(dat)

if(length(dat) != length(sampler_params))
  stop("ERROR")

for(i in 1:length(dat)){
  dat[[i]] <- cbind(dat[[i]], energy__ = sampler_params[[i]][,"energy__"])
}

bayesplot::mcmc_pairs(dat)
# There should be a negative relationship between lp__ and energy__.
# Primitive parameters that are correlated with the energy__ margin in the pairs
# plot are a good place to start thinking about reparameterizations.

# bayesplot::mcmc_pairs(sampler_params, pars = c("accept_stat__", "energy__"))

# Switch to checking summaries
# Dat <- read_tsv("MGYG-HGUT-03439.tsv.gz")
Dat <- summary(m1.stan)$summary %>%
  as.data.frame %>%
  rownames_to_column(var = "var") %>%
  as_tibble()
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

