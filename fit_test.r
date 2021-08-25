library(tidyverse)
library(HMVAR)

#' Frequency Increment Test (FIT)
#' 
#' From https://www.genetics.org/content/196/2/509
#'
#' @param Dat Tibble, must have site_id, freq, pt,
#' date, ref_id & ref_pos columns
#'
#' @return A tibble
#' @export
FIT <- function(Dat){
  Dat %>%
    dplyr::filter(freq != 1 & freq != 0) %>%
    # head(10000) %>%
    # filter(site_id %in% c("100035", "100073", "100087", "69024", "69024", "69027")) %>%
    # arrange(site_id, pt, date) %>%
    # group_by(site_id, pt) %>% 
    # summarise(Y = diff(freq))
    split(list(.$site_id, .$pt)) %>%
    purrr::map_dfr(function(d){
      
      if(nrow(d) < 3)
        return(NULL)
      
      pt <- unique(d$pt)
      site_id <- unique(d$site_id)
      ref_id <- unique(d$ref_id)
      ref_pos <- unique(d$ref_pos)
      
      d %>%
        dplyr::arrange(date) %>%
        dplyr::mutate(vprev = dplyr::lag(freq, n = 1, default = NA),
                      t = as.numeric(date - min(date))) %>%
        dplyr::mutate(tprev = dplyr::lag(t, n = 1, default = NA)) %>%
        dplyr::mutate(Y = (freq - vprev) / sqrt(2 * vprev * (1 - vprev) * (t - tprev))) %>%
        dplyr::summarise(Y.mean = mean(Y, na.rm = TRUE),
                         Y.s2 = var(Y, na.rm = TRUE),
                         L = length(Y) - 1) %>%
        dplyr::mutate(t.fi = Y.mean / sqrt(Y.s2 / L)) %>%
        dplyr::mutate(pval = 2 * (1 - pt(q = abs(t.fi), df = L - 1)),
                      pt = pt,
                      site_id = site_id,
                      ref_id = ref_id,
                      ref_pos = ref_pos)
      
    })
}


# Trying Fit test
# https://www.genetics.org/content/196/2/509


meta <- read_tsv("meta.txt")
meta
Dat <- read_midas_data("Fplautii_merged/")
Dat <- match_freq_and_depth(freq = Dat$freq,
                            depth = Dat$depth,
                            info = Dat$info %>%
                              select(site_id, ref_id, ref_pos),
                            map = meta,
                            depth_thres = 1)

Res <- FIT(Dat) %>%
  arrange(pval)
Res
write_tsv(Res, "output/Fplautii.FIT.tsv")