library(tidyverse)
library(HMVAR)

#' Calculate selection coefficient
#' 
#' Selection coefficient for alleles that change frequency between two
#' time points.
#' 
#' According to https://www.genetics.org/content/196/2/509
#'
#' @param Dat A tibble or data frame. It must have
#'
#' @return
#' @export
s_coefficient <- function(Dat){
  
  Dat %>%
    dplyr::filter(freq > 0 & freq < 1) %>%
    # dplyr::filter(site_id %in% c("100035", "100073", "100087", "69024", "69024", "69027")) %>%
    split(list(.$site_id, .$pt)) %>%
    purrr::map_dfr(function(d){
      
      if(nrow(d) < 2)
        return(NULL)
      
      pt <- unique(d$pt)
      site_id <- unique(d$site_id)
      ref_id <- unique(d$ref_id)
      ref_pos <- unique(d$ref_pos)
      
      d %>%
        # filter(date == max(date) | date == min(date)) %>%
        summarise(v0 = freq[ date == min(date)],
                  v1 = freq[ date == max(date)],
                  delta.t = as.numeric(max(date) - min(date))) %>%
        mutate(s = (1 / delta.t) * log( (v1 / (1 - v1)) * ((1 - v0) / v0) )) %>%
        mutate(pt = pt, site_id = site_id, ref_id = ref_id, ref_pos = ref_pos)
    })
}


meta <- read_tsv("meta.txt")
meta
Dat <- read_midas_data("Fplautii_merged/")
Dat <- match_freq_and_depth(freq = Dat$freq,
                            depth = Dat$depth,
                            info = Dat$info %>%
                              select(site_id, ref_id, ref_pos),
                            map = meta,
                            depth_thres = 1)
Dat

Res <- s_coefficient(Dat) %>%
  arrange(desc(abs(s)))
Res
write_tsv(Res, "output/Fplautii.s_coefficient.tsv")