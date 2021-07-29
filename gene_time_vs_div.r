#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
# This file is part of .
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

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Get number of SNPs between pairs of samples per",
                        "gene, and produce table with number of differences",
                        "from oldest sample"))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Input directory"),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("Mapping file, must have columns sample,",
                                 "pt, and date"),
                    type = "character")
  
  
  # Optional arguments
  p <- add_argument(p, "--type",
                     help = paste("Either a single MIDAS dir or a directory",
                                  "of multiple (multi) MIDAS directories"),
                     type = "character",
                     default = "single")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--snp_type",
                    help = paste("Following inStrain definitions either",
                                 "popSNPs or conSNPs"),
                    type = "character",
                    default = "popSNPs")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$snp_type %in% c("popSNPs", "conSNPs")){
    args$distfun <- args$snp_type
  }else{
    stop("ERROR: bad --snp_type", call. = TRUE)
  }
  
  if(args$type == "multi"){
    cat("Processing multiple species directories...\n")
    args$input <- list.dirs(args$input, recursive = FALSE)
  }else if(args$type == "single"){
    cat("Processing single species directory...\n")
  }else{
    stop("ERROR: invalid --type", call. = TRUE)
  }
  
  return(args)
}

#' Count inStrain popSNPs
#'
#' @param x 
#'
#' @return
#' @export
popSNPs <- function(x){
  x[ !(x == 1 | x == 0) ] <- NA # Keep fixed differences
  
  apply(x,1,function(vec,x.t){
    colSums(vec != x.t, na.rm = TRUE)
  }, x.t = t(x)) %>% as.dist()
}

#' Count inStrain conSNPs
#'
#' @param x 
#'
#' @return
#' @export
conSNPs <- function(x){
  x[ x == 0.5 ] <- NA # ties are not conSNPs
  x <- round(x, digits = 0) # get consensus
  
  apply(x,1,function(vec,x.t){
    colSums(vec != x.t, na.rm = TRUE)
  }, x.t = t(x)) %>% as.dist()
}

# Read arguments
args <- process_arguments()
# args <- list(input = "MGYG-HGUT-00099/",
#              map = "hct_quickmap.txt",
#              distfun = "popSNPs")

# Load more packages
library(tidyverse)
library(HMVAR)

# Read map
meta <- read_tsv(args$map,
                 col_types = cols(sample = col_character(),
                                  pt = col_character(),
                                  date = col_date(format = "%Y-%m-%d")))

# test_genes <- c("GUT_GENOME000518_00001",
#                 "GUT_GENOME000518_00002",
#                 "GUT_GENOME000518_00003",
#                 "GUT_GENOME000518_00004",
#                 "GUT_GENOME000518_00005",
#                 "GUT_GENOME000518_00006",
#                 "GUT_GENOME000518_00007",
#                 "GUT_GENOME000518_00008",
#                 "GUT_GENOME000518_00009",
#                 "GUT_GENOME000518_000010")

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

cat("\n\n=====================\n")
for(d in args$input){
  spec <- basename(d)
  cat("Processing species", spec, "\n")

  cat("\treading data...\n")
  Dat <- read_midas_data(d, cds_only = TRUE, map = meta)
  # print(Dat$info)

  cat("\tGetting SNP counts...\n") 
  Res <- Dat$info %>%
    split(.$gene_id) %>%
    map(function(i, Dat,
                 min_depth = 10,
                 prop = 0.8,
                 min_snps = 20,
                 distfun = "dist"){
      
      distfun <- match.fun(distfun)
      # cat(unique(i$gene_id), "\n")
      if(nrow(i) < min_snps){
        return(NULL)
      }
      
      # Extract data from gene
      d <- Dat$depth %>%
        filter(site_id %in% i$site_id)
      d_sites <- d$site_id
      f <- Dat$freq %>%
        filter(site_id %in% i$site_id)
      f_sites <- f$site_id
      
      if(!all(d_sites %in% f_sites) || !all(f_sites %in% d_sites))
        stop("ERROR: inconsistent sites between freq and depth", call. = TRUE)
      
      # Convert data to matrices
      d <- d %>%
        select(-site_id) %>%
        as.matrix()
      row.names(d) <- d_sites
      f <- f %>%
        select(-site_id) %>%
        as.matrix()
      row.names(f) <- f_sites
      
      # Homogenize matrices
      f <- f[row.names(d), colnames(d)]
      
      # Filter by depth
      f[ d < min_depth ] <- NA
      
      # Filter by missing sites
      f <- f[ , colSums(!is.na(f)) >= nrow(f) * prop, drop = FALSE ]
      
      if(ncol(f) <= 10){
        return(NULL)
      }
      
      # Calculate divergence
      dis <- distfun(t(f))
      # return(dis)
      
      # Third attempt counting distance from first sample
      first_sample <- (meta %>%
                         filter(sample %in% attr(dis, "Labels")) %>%
                         arrange(date) %>%
                         select(sample) %>%
                         unlist)[1] %>%
        as.character
      divs <- as.matrix(dis)[,first_sample]
      divs <- tibble(sample = names(divs), divergence = as.numeric(divs))
      
      return(list(dis = dis, divs = divs))
      
    }, Dat = Dat, distfun = args$distfun) %>%
    compact()
  
  if(length(Res) == 0){
      cat("\t>No genes met conditions. Skipping\n")
      next
  }

  filename <- file.path(args$outdir, paste0(spec,".rdat"))
  save(Res, file = filename)
  
  cat("\tGetting table of 'divergences'...\n")
  Divs <- Res %>%
    map_dfr(function(l, meta){
      l$divs %>%
        left_join(meta, by = "sample")
    }, meta = meta, .id = "gene_id") %>%
    split(.$gene_id) %>%
    map_dfr(function(d){
      d$days <- as.numeric(d$date - min(d$date))
      d
    })
  filename <- file.path(args$outdir, paste0(spec,".tsv"))
  write_tsv(Divs, filename)
 
  cat("\tPlotting overall 'divergences'...\n")
  p1 <- Divs %>%
    ggplot(aes(x = days, y = divergence)) +
    geom_point(aes(col = pt), show.legend = FALSE) +
    geom_smooth(method = "lm", formula = y ~ x) +
    # scale_y_log10() +
    # scale_x_log10() +
    ylab(label = "# SNPs") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90))
  filename <- file.path(args$outdir, paste0(spec,".png"))
  ggsave(filename, width = 6, height = 4, dpi = 150)
  
}
cat("=====================\n")







# Divs %>%
#   group_by(sample, pt, date, days) %>%
#   summarise(divergence = sum(divergence),
#             .groups = 'drop') %>%
#   arrange(pt, days) %>%
#   # print(n = 100)
#   ggplot(aes(x = date, y = log10(divergence + 1))) +
#   # facet_wrap(~ pt, scales = "free") +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_classic()
# ggsave("date_vs_popSNPs_species.jpeg", width = 6, height = 4, dpi = 150)
# 



