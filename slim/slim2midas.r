library(tidyverse)

#' Read genotype data in ms format
#'
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
read_ms <- function(file){

  ms <- read_lines(file)
  
  # Get first two lines
  # Not actually mandatory in SLiM output
  # ms_cmd <- ms[[1]]
  # ms_seeds <- ms[[2]]

  Res <- NULL
  Geno <- NULL
  i <- 1
  sample_i <- 0
  while(i <= length(ms)){
    # cat(i, "\n")
    if( str_detect(ms[[i]], "^//") ){
      sample_i <- sample_i + 1
      
      if(sample_i > 1){
        stop("ERROR: script only prepared for ms files with 1 sample",
             call. = TRUE)
      }
      
      # Finding sample beginning
      i <- i + 1 
      next
    }
    
    if(sample_i > 0){
      # Only if we are withing the samples region of the file
      if( str_detect(ms[[i]], "^segsites:") ){
        n_sites <- str_remove(ms[[i]], "^segsites:")
        n_sites <- as.numeric(n_sites)
      }else if( str_detect(ms[[i]], "^positions:") ){
        positions <- str_remove(ms[[i]], "^positions:")
        positions <- str_split(trimws(positions, which = "both"), " ")[[1]]
        positions <- as.numeric(positions)
        
        if(length(positions) != n_sites){
          stop("ERROR: number of positions doesn't match number of sites",
               call. = TRUE)
        }
        
        Res[[sample_i]] <- tibble(site_id = as.character(1:n_sites),
                                  ref_id = "chr1",
                                  ref_pos = positions)
        Geno[[sample_i]] <- tibble(site_id = Res[[sample_i]]$site_id)
        
      }else if(trimws(ms[[i]]) == ""){
        cat("Skipping line", i, "\n")
      } else{
        genotype <- as.numeric(str_split(ms[[i]], "")[[1]])
        
        if(length(genotype) != n_sites){
          stop("ERROR: number of genotyped positions doesn't match number of sites",
               call. = TRUE)
        }
        
        Geno[[sample_i]][ paste0("g", ncol(Geno[[sample_i]])) ] <- genotype
      }
    }else{
      cat("Skipping line", i, "\n")
    }

    
    i <- i + 1
  }
  
  return(list(info = Res[[1]], freq = Geno[[1]]))
}



args <- list(start_dir = "standing_variation/",
             sim_dir = "sim_x/",
             n_genomes = 10,
             seed = 2308123)


set.seed(args$seed)
dat <- read_ms("standing_variation/standing_variation.ms")
# genomes <- sample(colnames(dat$freq)[-1], size = args$n_genomes)
genomes <- colnames(dat$freq)[-1]

Info <- dat$info %>%
  left_join(dat$freq %>%
  select(site_id, all_of(genomes)) %>%
  pivot_longer(-site_id) %>%
  group_by(site_id) %>%
  summarise(maf = mean(value),
            .groups = 'drop'),
  by = "site_id") %>%
  mutate(ref_pos = floor(ref_pos * 1e6)) %>%
  rename(g_0 = maf)
Info



pop_dir <- file.path(args$sim_dir, "pop_1")

pop_id <- basename(pop_dir)
  
# Process population directory filenames
pop_files <- list.files(pop_dir, recursive = FALSE, full.names = T)
ms_files <- pop_files[ basename(pop_files) %>% str_detect("[.]ms$") ]
generations <- basename(ms_files) %>% str_remove("^gen_") %>% str_remove("[.]ms$") %>% as.numeric()



Info <- NULL
Freq <- NULL
for(g in generations){
  g <- 10
  # Reading ms file for current generation sample
  dat <- read_ms( file.path(pop_dir, paste0("gen_", g, ".ms")) )
  
  dat$info <- dat$info %>%
    left_join(dat$freq[1:(args$n_genomes + 1)] %>%
                pivot_longer(-site_id, names_to = "genome", values_to = "genotype") %>%
                group_by(site_id) %>%
                summarise(maf = mean(genotype),
                          .groups = 'drop'),
              by = "site_id") %>%
    mutate(ref_pos = floor(1e6 * ref_pos))
  # dat$info
  
  gen_label <- paste0("gen_", g)
  
  # Read info from same generation and homogenize
  info <- HMVAR::read_midas_info(file.path(pop_dir, paste0("gen_", g, "_snp_info.txt")))
  # info
  
  info %>%
    full_join(dat$info %>%
                select(ref_id, ref_pos, maf)) %>%
    rename(`gen_label` = maf)
  
  
  
}  
  








