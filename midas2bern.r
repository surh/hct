library(tidyverse)

# This is taking from orginal pipeline and cleaning & organizing. It is only going
# to take two timepoints. Should simplify future analysis.
# This can probably go in HMVAR, but bern model needs to be its own thing

args <- list(midas_dir = "../../exp/2021/2021-04-19.monotonic_hct/MGYG-HGUT-00099/",
             depth_thres = 5,
             n_thres = 3)


meta <- tibble(pt = letters[1:5],
               start = c("P316", "P102", "P16", "P373", "P219"),
               end = c("P246", "P215", "P228", "P391", "P325" ))
meta




midas2sitesdist <- function(midas_dir, meta, outdir = "./",
                            write_tables = FALSE,
                            depth_thres = 5,
                            npop_thres = 3){
  
  if(missing(midas_dir) || missing(meta)){
    stop("ERROR: You must provide a midas_dir & meta(data) tibble.",
         call. = TRUE)
  }
  if(!is_tibble(meta)){
    stop("ERROR: meta must be a tibble", call. = TRUE)
  }
  if(!all(c("pt", "start", "end") %in% colnames(meta))){
    stop("ERROR: columns pt, start &* end must be present in meta",
         call. = TRUE)
  }
  
  ### Read and process midas data
  midas_map <- meta %>%
    pivot_longer(-pt,values_to = "sample", names_to = "timepoint") %>%
    rename(Group = pt)
  Dat <- HMVAR::read_midas_data(midas_dir = midas_dir,
                                map = midas_map)
  Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq,
                                     depth = Dat$depth,
                                     info = Dat$info %>%
                                       dplyr::select(site_id, ref_id, ref_pos),
                                     map = midas_map %>%
                                       dplyr::rename(pt = Group), 
                                     depth_thres = depth_thres)
  
  cat("Calculating change at every position in every population...\n")
  Dat <- Dat %>%
    split(.$pt) %>%
    map_dfr(function(d){
      d %>%
        pivot_wider(id_cols = c("site_id", "ref_id", "ref_pos", "pt"),
                    names_from = "timepoint",
                    values_from = c("sample", "freq", "depth")) %>%
        filter(!(is.na(freq_start) | is.na(freq_end)))
    }) %>%
    group_by(site_id) %>%
    filter(n() >= npop_thres) %>%
    ungroup() %>%
    mutate(maf_change = freq_end - freq_start)
  
  cat("Calculating per-site distribution\n")
  Sites <- Dat %>%
    group_by(site_id) %>%
    summarise(n_decrease = sum(maf_change < 0),
              n_equal = sum(maf_change == 0),
              n_increase = sum(maf_change > 0),
              n_patients = length(pt),
              .groups = "drop")
  
  cat("Calculating per-population distribution")
  Pts <- Dat %>%
    group_by(pt) %>%
    summarise(n_decrease = sum(maf_change < 0),
              n_equal = sum(maf_change == 0),
              n_increase = sum(maf_change > 0),
              n_sites = length(site_id),
              .groups = 'drop') %>%
    # filter(n_sites >= args$prop_thres * max(n_sites)) %>%
    mutate(p_increase = n_increase / n_sites,
           p_decrease = n_decrease / n_sites,
           p_change = (n_increase + n_decrease) / n_sites)
  
  
  if(write_tables){
    # Prepare output dir
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
    cat("Writing tables..\n")
    write_tsv(Sites, file.path(outdir, "sites.tsv"))
    write_tsv(Pts, file.path(outdir, "pops.tsv"))
  }
  
  return(list(Sites = Sites, Pops = Pts))
}


Res <- midas2sitesdist(midas_dir = args$midas_dir, meta = meta,
                       outdir = "test/",
                       write_tables = FALSE,
                       depth_thres = 5, npop_thres = 3)
Res


