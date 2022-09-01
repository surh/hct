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



#' Transition probabilities for modified Moran process
#' 
#' Given a system state = (i.e. pop size & number of individuals of
#' genotype A), & fitness of both genotypes, it returns a vector of
#' of probabilities that the number of individuals of genotype A decreases,
#' remains the same, or increases by 1.
#' 
#' When no genotypes has reached fixation, the probailities are those of
#' the Standard Moran Process with selection
#' (https://en.wikipedia.org/wiki/Moran_process#Selection).
#' 
#' When one of the genotypes has reached fixation. There is a probability
#' (p_imm) that the other genotype will be reseeded via immigration.
#'
#' @param x Number of individuals with genotype A.
#' @param N Poputlation size
#' @param f Fitness of genotype A.
#' @param g Fitness of genotype B.
#' @param p_imm Probability of reseeding of extinct genotype via immigration
#' (See details).
#'
#' @return A numeric vector of probabilities that the number of individuals
#' of genotype A (x) decreases by 1, remains equal, or increases by 1.
#' @export
#' @references 
#' https://en.wikipedia.org/wiki/Moran_process
transition_probs <- function(x, N, f, g, p_imm){
  if(x == 0){
    p_dec <- 0
    p_inc <- p_imm
  }else if (x == N){
    p_dec <- p_imm
    p_inc <- 0
  }else{
    pop_fitness <- (f * x) + (g * (N - x))
    p_dec <- (g * (N - x) / pop_fitness) * (x / N)
    p_inc <- (f * x / pop_fitness) * ((N - x) / N)
  }
  
  p_eq <- 1 - p_dec - p_inc 
  
  return(c(p_dec, p_eq, p_inc))
}

#' Simulate site under modified Moran process.
#'
#' @param x_0 Starting number of individuals with genotype A
#' @param N Population size
#' @param T Number of Moran Process time steps
#' @param f Fitness of genotype A
#' @param g Fitness of genotype B
#' @param p_imm Probability of re-seeding via immigration after extinction
#' of a genotype.
#'
#' @return A vector with the simulated system state (i.e. number of
#' individuals with genotype A) over time.
#' @export
#' @author Sur Herrera Paredes
sim_site <- function(x_0, N, T, f, g, p_imm){
  # set.seed(12903242)
  # Res <- NULL
  X <- x_0
  for(i in 1:T){
    P <- transition_probs(x = X[i], N = N, f = f, g = g, p_imm = p_imm)
    X[i + 1] <- X[i] + sample(c(-1, 0, 1), size = 1, prob = P)
  }
  
  return(X)
}

#' Simulate population
#' 
#' Simulates the evolution of multiple sites in a population with a 
#' modified Moran Process
#'
#' @param nsites Number of sites (i.e. variants) in the population.
#' @param pop_size Population site
#' @param x_0 Starting number of individuals with genotype A. NOTE: All sites
#' start with the same genotypic state in the population.
#' @param p_imm Probability of re-seeding given extinction of a genotype in
#' a given site in a population.
#' @param f Fitness for genotype A in selected sites. Also used as fitness
#' for all genotypes in non selected sites.
#' @param g fitnesss of genotype B in selected sites
#' @param prop_selected Proportion of sites under selection
#' @param T Number of steps of the Moran process.
#'
#' @return
#' @export
sim_population <- function(nsites, pop_size = 100,
                           x_0 = 99,
                           p_imm = 0.05, f = 1,
                           g = 1.05,
                           prop_selected = 0.01, T = 100){
  n_selected <- round(nsites * prop_selected)
  
  freqs <- NULL
  
  if(n_selected > 0){
    for(i in 1:n_selected){
      freqs <- rbind(freqs, sim_site(x_0 = x_0, N = pop_size, T = T,
                                     f = f, g = g,
                                     p_imm = p_imm) / pop_size)
    }
  }
  
  if(n_selected < nsites){
    for(i in (n_selected + 1):nsites){
      freqs <- rbind(freqs, sim_site(x_0 = x_0, N = pop_size, T = T,
                                     f = f, g = f,
                                     p_imm = p_imm) / pop_size)
    }
    
  }
  
  colnames(freqs) <- paste0("t", 0:T)
  freqs <- as_tibble(freqs)
  freqs$site_id <- paste0("s", 1:nsites)
  freqs <- freqs %>%
    select(site_id, everything())
  
  info <- tibble(site_id = freqs$site_id,
                 selected = rep(c(TRUE,FALSE),
                                times = c(n_selected, nsites-n_selected)))  
  
  return(list(info = info, freq = freqs))
  
}

get_site_maf_changes <- function(freq, t_0 = 0, t_n = 100){
  t_0 <- paste0("t", t_0)
  t_n <- paste0("t", t_n)
  freq %>%
    select(site_id, v0 = all_of(t_0), v1 = all_of(t_n))
    # mutate(maf_change = tn - t0) %>%
    # select(site_id, maf_change)
}


