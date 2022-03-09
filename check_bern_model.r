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

my_summary <- function(x){
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75))
  tibble(min = min(x),
         Q1 =qs[1],
         median = qs[2],
         mean = mean(x),
         Q3 = qs[3],
         max = max(x))
}



# setwd("/cashew/users/sur/exp/fraserv/2022/today2/output")
# args <- list(stan_model = "/cashew/users/sur/exp/fraserv/2022/today3/preFOS_posFOS_control_3/output/stan_models/MGYG-HGUT-00144.stan.rdat",
#              output = "checks",
#              max_treedepth = 10)

knitr::opts_chunk$set(fig.width = 10)
args <- list(stan_model = opts[1],
             max_treedepth = as.numeric(opts[2]))

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


# ShinyStan
# shinystan::launch_shinystan(m1.stan)



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

Dat %>%
  mutate(var_type = str_remove(var, "[\\[][0-9]+[\\]]")) %>%
  group_by(var_type) %>%
  summarise(avg = mean(mean),
            min = min(mean),
            max = max(mean))

################ ps ######################
dat <- Dat %>%
  mutate(var_type = str_remove(var, "[\\[][0-9]+[\\]]")) %>%
  filter(var_type == "p") %>%
  select(`50%`) %>%
  unlist %>% as.vector
beta.fit <- MASS::fitdistr(x = dat, densfun = "beta",
                           start = list(shape1 = mean(dat) * 10,
                                        shape2 = (1 - mean(dat)) * 10))
print(beta.fit)

# hist(rbeta(n = 1000, shape1 = mean(dat) * 10, shape2 = 1 - mean(dat) * 10))
# hist(rbeta(n = 1000, shape1 = beta.fit$estimate['shape1'], shape2 = beta.fit$estimate['shape2']))

# m = shape1 / v
# v = shape2 / (1 - m)
# 
# m = shape1 / ( shape2 / (1 -m) )
# m = shape1 * (1 - m) / shape2
# m * shape2 = shape1 - shape1 * m
# m * (shape1 + shape2) = shape1
# 
# m = shape1 / (shape1 + shape2)
# v = shape2 / (1 - m)
m <- beta.fit$estimate['shape1'] / (beta.fit$estimate['shape1'] + beta.fit$estimate['shape2'])
v <- beta.fit$estimate['shape2'] / (1 - m)
print(m)
print(v)

print(mean(dat))

tibble(var = dat,
       dbeta =  dbeta(x = dat,
                      shape1 = m * v,
                      shape2 = (1 - m)*v)) %>%
  ggplot(aes(x = var)) +
  geom_density(col = "red", fill = "red", alpha = 0.3) +
  geom_line(aes(y = dbeta)) +
  # ylim(c(0, 0.03)) +
  theme_classic()


rm(dat, m, v, beta.fit)
################ qs ######################
dat <- Dat %>%
  mutate(var_type = str_remove(var, "[\\[][0-9]+[\\]]")) %>%
  filter(var_type == "q") %>%
  select(`50%`) %>%
  unlist %>% as.vector
beta.fit <- MASS::fitdistr(x = dat, densfun = "beta",
                           start = list(shape1 = mean(dat) * 10,
                                        shape2 = (1 - mean(dat)) * 10))
print(beta.fit)

# hist(rbeta(n = 1000, shape1 = mean(dat) * 10, shape2 = 1 - mean(dat) * 10))
# hist(rbeta(n = 1000, shape1 = beta.fit$estimate['shape1'], shape2 = beta.fit$estimate['shape2']))

# m = shape1 / v
# v = shape2 / (1 - m)
# 
# m = shape1 / ( shape2 / (1 -m) )
# m = shape1 * (1 - m) / shape2
# m * shape2 = shape1 - shape1 * m
# m * (shape1 + shape2) = shape1
# 
# m = shape1 / (shape1 + shape2)
# v = shape2 / (1 - m)
m <- beta.fit$estimate['shape1'] / (beta.fit$estimate['shape1'] + beta.fit$estimate['shape2'])
v <- beta.fit$estimate['shape2'] / (1 - m)
print(m)
print(v)

print(mean(dat))

tibble(var = dat,
       dbeta =  dbeta(x = dat,
                      shape1 = m * v,
                      shape2 = (1 - m)*v)) %>%
  ggplot(aes(x = var)) +
  geom_density(col = "red", fill = "red", alpha = 0.3) +
  geom_line(aes(y = dbeta)) +
  # ylim(c(0, 0.03)) +
  theme_classic()




# n <- 1000
# m <- 0.49
# v <- 10
# hist(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v), breaks = 20)
# range(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v))
# 
# m <- 0.476
# v <- 30
# hist(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v), breaks = 20)
# range(rbeta(n = n, shape1 = m * v, shape2 = (1 - m) * v))

