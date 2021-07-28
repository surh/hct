library(tidyverse)
library(HMVAR)

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

args <- list(dir = "MGYG-HGUT-00099/",
             map = "hct_quickmap.txt",
             distfun = "popSNPs")



meta <- read_tsv(args$map)
meta <- meta %>%
  separate(col = "group", into = c("pt","date"), sep = "[.]") %>%
  mutate(date = replace(date, date == "2017-04-11-bp",  "2017-04-11")) %>%
  mutate(date = replace(date, date == "2017-04-11-sp",  "2017-04-11")) %>%
  mutate(date = replace(date, date == "NA",  NA)) %>%
  # select(date) %>% unlist %>% unique
  filter(!is.na(date))

test_genes <- c("GUT_GENOME000518_00001",
                "GUT_GENOME000518_00002",
                "GUT_GENOME000518_00003",
                "GUT_GENOME000518_00004",
                "GUT_GENOME000518_00005",
                "GUT_GENOME000518_00006",
                "GUT_GENOME000518_00007",
                "GUT_GENOME000518_00008",
                "GUT_GENOME000518_00009",
                "GUT_GENOME000518_000010")
test_genes <- NULL
Dat <- read_midas_data(args$dir, cds_only = TRUE, map = meta, genes = test_genes)
Dat


Res <- Dat$info %>%
  split(.$gene_id) %>%
  map(function(i, Dat,
               min_depth = 10,
               prop = 0.8,
               min_snps = 20,
               distfun = "dist"){
    # i <- Dat$info %>% filter(gene_id == "GUT_GENOME000518_00164")
    # min_depth <- 10
    # prop <- 0.8
    # min_snps <- 20
    # distfun <- args$dir
    
    distfun <- match.fun(distfun)
    cat(unique(i$gene_id), "\n")
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
    
    # Second attempt using tree to count divergence
    # tre <- ape::njs(dis)
    # tre <- phytools::midpoint.root(tre)
    # tips <- tre$tip.label
    # divs <- castor::get_all_distances_to_root(tre)
    # 
    # divs <- tibble(sample = tips, divergence = divs[1:length(tips)])
    # 
    # return(list(tre = tre, divs = divs))
    
    # First attempt using matching freq and depth, probably not general and slow
    # match_freq_and_depth(freq = f, depth = d,
    #                      info = i %>%
    #                        select(site_id, gene_id),
    #                      map = meta,
    #                      depth_thres = 1)
    
  }, Dat = Dat, distfun = args$distfun) %>%
  compact()
Res




Res %>%
  map_dfr(function(l, meta){
    l$divs %>%
      left_join(meta, by = "sample")
  }, meta = meta %>%
    mutate(date = parse_date(date, format = "%Y-%m-%d")), .id = "gene_id") %>%
  split(.$gene_id) %>%
  map_dfr(function(d){
    d$days <- as.numeric(d$date - min(d$date))
    d
  }) %>%
  # filter(gene_id == "GUT_GENOME000518_00001") %>%
  ggplot(aes(x = days, y = divergence)) +
  # facet_wrap(~ gene_id, scales = "free") +
  geom_point(aes(col = pt), show.legend = FALSE) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

# Res$GUT_GENOME000518_00001$divs %>%
#   arrange(desc(divergence))
# plot(Res$GUT_GENOME000518_00001$tre)
# 
# Res$GUT_GENOME000518_00001[50:102, ]
# Res$GUT_GENOME000518_00001 %>% as.matrix()
# hclust(Res$GUT_GENOME000518_00001)
# 
# plot(ape::njs(Res$GUT_GENOME000518_00001))
# plot(Res$GUT_GENOME000518_00001)
# plot(Res$GUT_GENOME000518_00002)
# plot(Res$GUT_GENOME000518_00003)
# plot(Res$GUT_GENOME000518_00004)
# plot(Res$GUT_GENOME000518_00005)
# plot(Res$GUT_GENOME000518_00008)




