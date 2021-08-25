library(tidyverse)
library(dndscv)
library(seqinr)
# data(dataset_simbreast)
# res <- dndscv(mutations)

# Read gene defintions
gff <- read_tsv("MGYG-HGUT-00099/MGYG-HGUT-00099.gff",
                col_names = FALSE,
                col_types = cols(X1 = col_character(),
                                 X3 = col_character(),
                                 X4 = col_number(),
                                 X5 = col_number(),
                                 X7 = col_character(),
                                 X9 = col_character()))
gff$gene.id <- gff$X9 %>%
  map_chr(function(x){
    id <- str_split(x, pattern = ";")[[1]][1]
    id <- str_remove(id, "^ID=")
    id
  })
gff

# Annotations provide gene names
annot <- HMVAR::read_eggnog("MGYG-HGUT-00099/MGYG-HGUT-00099_eggNOG.tsv",
                            format = "uhgg")

# Create cds object for cds file
cds <- gff %>%
  filter(X3 == "CDS") %>%
  transmute(gene.id,
            cds.id = gene.id,
            chr = X1,
            chr.coding.start = X4,
            chr.coding.end = X5,
            strand = X7) %>%
  mutate(length = (chr.coding.end - chr.coding.start + 1)) %>%
  mutate(cds.start = 1,
         cds.end = length,
         strand_num = 1) %>%
  mutate(strand_num = replace(strand_num, strand == '-', -1)) %>%
  left_join(annot %>%
              select(gene.id = query_name, gene.name = predicted_gene_name),
            by = "gene.id") %>%
  select(gene.id, gene.name, cds.id, chr, chr.coding.start, chr.coding.end,
         cds.start, cds.end, length, strand = strand_num) %>%
  mutate(gene_num = gene.id %>%
           str_split(pattern = "_") %>%
           map_chr(~ .x[3])) %>%
  mutate(gene.name = replace(gene.name, !is.na(gene.name), paste(gene.name[!is.na(gene.name)], gene_num[!is.na(gene.name)], sep = ":"))) %>%
  mutate(gene.name = replace(gene.name, is.na(gene.name), gene.id[ is.na(gene.name) ])) %>%
  select(-gene_num)
write_tsv(cds, "refdb/Fplautii.cdsfile.tsv")
name2id <- cds %>%
  select(gene_name = gene.name, gene_id = gene.id)

### Build new reference file so I can include all sites
# Read midas data
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
site2gene <- Dat$info %>% select(site_id, gene_id)

# Identify positions to change in reference
to_change <- Dat$info %>%
  filter(ref_allele != major_allele) %>%
  select(ref_id, ref_pos, ref_allele, major_allele, minor_allele) %>%
  mutate(major_allele = tolower(major_allele),
         minor_allele = tolower(minor_allele),
         ref_allele = tolower(ref_allele))

# Read fasta
seq <- read.fasta("MGYG-HGUT-00099/MGYG-HGUT-00099.fna")

# Change positions
for(i in 1:nrow(to_change)){
  seq[[to_change$ref_id[i]]][to_change$ref_pos[i]] <- to_change$major_allele[i]
}

# Write custom ref for dnds
write.fasta(seq, names = names(seq), file.out = "refdb/Fplautii.fasta")


buildref(cdsfile = "refdb/Fplautii.cdsfile.tsv",
         genomefile = "refdb/Fplautii.fasta", 
         outfile = "refdb/Fplautii.rda",
         numcode = 1, excludechrs = NULL,
         onlychrs = NULL, useids = FALSE)


# Read p_directional
pdir <- read_tsv("p_directional/MGYG-HGUT-00099.tsv.gz", 
                 col_types = cols(site_id = col_character()))
pdir <- pdir %>%
  select(site_id, p_directional) %>%
  left_join(site2gene, by = "site_id") %>%
  group_by(gene_id) %>%
  summarise(p_directional = max(p_directional),
            .groups = 'drop')
pdir

#### Analysis
# Pick at random among samples with repeated SNPs

info <- Dat$info %>%
  # filter(major_allele == ref_allele) %>% # not needed later
  select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                            depth_thres = 1)


set.seed(5435)
Dat <- Dat %>%
  filter(!is.na(ref)) %>%
  filter(freq >= 0.8) %>%
  split(.$site_id) %>%
  map_dfr(function(d){
    ii <- which(d$freq == max(d$freq))
    if(length(ii) > 1){
      ii <- sample(ii, size = 1)
    }
    
    d[ii, ]
  })
Dat

Dat <- Dat %>%
  select(sampleID = sample, chr, pos, ref, mut)

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)


res$sel_loc %>% as_tibble() %>%
  filter(qall_loc < 0.1)

res$sel_cv %>%
  as_tibble() %>%
  filter(qallsubs_cv < 0.1) %>%
  arrange(desc(wmis_cv))

res$sel_cv %>%
  as_tibble() %>%
  arrange(desc(wmis_cv)) %>%
  print(n = 50)

Res <- res$sel_cv %>%
  as_tibble() %>%
  select(gene_name, n_syn, n_mis, n_non, wmis_cv, wnon_cv, pmis_cv, ptrunc_cv, pallsubs_cv, qmis_cv, qtrunc_cv, qallsubs_cv) %>%
  # filter(qallsubs_cv < 0.1) %>% arrange(desc(wmis_cv)) %>%
  left_join(name2id, by = "gene_name") %>%
  inner_join(pdir, by = "gene_id")
Res

# Res %>%
#   ggplot(aes(x = p_directional / (1 - p_directional), y = n_syn)) +
#   geom_point(alpha = 0.1, size = 2) +
#   theme_classic()

Res
Res %>%
  mutate(pmis_cv = -log10(pmis_cv),
         ptrunc_cv = -log10(ptrunc_cv),
         pallsubs_cv = -log10(pallsubs_cv),
         qmis_cv = -log10(qmis_cv),
         qtrunc_cv = -log10(qtrunc_cv),
         qallsubs_cv = -log10(qallsubs_cv)) %>%
  mutate(pmis_cv = replace(pmis_cv, pmis_cv == Inf, min(pmis_cv[pmis_cv != Inf])),
         ptrunc_cv = replace(ptrunc_cv, ptrunc_cv == Inf, min(ptrunc_cv[ptrunc_cv != Inf])),
         pallsubs_cv = replace(pallsubs_cv, pallsubs_cv == Inf, min(pallsubs_cv[pallsubs_cv != Inf])),
         qmis_cv = replace(qmis_cv, qmis_cv == Inf, min(qmis_cv[qmis_cv != Inf])),
         qtrunc_cv = replace(qtrunc_cv, qtrunc_cv == Inf, min(qtrunc_cv[qtrunc_cv != Inf])),
         qallsubs_cv = replace(qallsubs_cv, qallsubs_cv == Inf, min(qallsubs_cv[qallsubs_cv != Inf]))) %>%
  mutate(wmis_cv = log2(wmis_cv),
         wnon_cv = log2(wnon_cv)) %>%
  pivot_longer(cols = c(-gene_name, -gene_id, -p_directional),
               names_to = "statistic", values_to = "value") %>%
  filter(!(statistic == "wmis_cv" & value == -Inf)) %>%
  filter(!(statistic == "wnon_cv" & value == -Inf)) %>%
  
  ggplot(aes(x = log2(p_directional / (1 - p_directional)), y = value)) +
  facet_wrap(~ statistic, scales = "free_y", ncol = 3) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic()
  

Res %>%
  ggplot(aes(x = p_directional >= 0.8, y = wmis_cv)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(position = position_jitter(width = 0.1)) +
  scale_y_log10() +
  AMOR::theme_blackbox()

Res %>%
  filter(wnon_cv > 0) %>%
  ggplot(aes(x = p_directional >= 0.8, y = wnon_cv)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point(position = position_jitter(width = 0.1)) +
  scale_y_log10() +
  AMOR::theme_blackbox()

# Include all samples and SNPs with fixed minor allele
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
info <- Dat$info %>%
  select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                   depth_thres = 1)
Dat <- Dat %>%
  filter(freq >= 0.8)

Dat <- Dat %>%
  select(sampleID = sample, chr, pos, ref, mut)
Dat

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)


# All samples and SNPs with present minor alleles
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
info <- Dat$info %>%
  select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                   depth_thres = 1)
Dat <- Dat %>%
  filter(freq > 0)

Dat <- Dat %>%
  select(sampleID = sample, chr, pos, ref, mut)
Dat

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)


# Use only true singletons
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
info <- Dat$info %>%
  select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                   depth_thres = 1)

Dat <- Dat %>%
  filter(depth > 0) %>%
  split(.$site_id) %>%
  map_dfr(function(d){
    d <- d %>%
      filter(freq > 0)
    
    if(nrow(d) == 1){
      return(d)
    }else{
      return(NULL)
    }
  })
Dat


Dat <- Dat %>%
  select(sampleID = sample, chr, pos, ref, mut)
Dat

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)


# Use only fixed singletons
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
info <- Dat$info %>%
  select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                   depth_thres = 1)

Dat <- Dat %>%
  filter(depth > 0) %>%
  split(.$site_id) %>%
  map_dfr(function(d){
    d <- d %>%
      filter(freq > 0)
    
    if(nrow(d) == 1){
      if(d$freq >= 0.8){
        return(d)
      }
    }
    return(NULL)
  })
Dat


Dat <- Dat %>%
  select(sampleID = sample, chr, pos, ref, mut)
Dat

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)

# All sites one dummy sample
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
Dat <- Dat$info %>%
  transmute(sampleID = "dummySample", chr = ref_id, pos = ref_pos,
            ref = major_allele, mut = minor_allele)
Dat

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)

# All true singletons one dummy sample
Dat <- HMVAR::read_midas_data("Fplautii_merged/", cds_only = TRUE)
info <- Dat$info %>%
  select(site_id, chr = ref_id, pos = ref_pos, ref = major_allele, mut = minor_allele)
Dat <- HMVAR::match_freq_and_depth(freq = Dat$freq, depth = Dat$depth, info = info,
                                   depth_thres = 1)


Dat <- Dat %>%
  filter(depth > 0) %>%
  split(.$site_id) %>%
  map_dfr(function(d){
    d <- d %>%
      filter(freq > 0)
    
    if(nrow(d) == 1){
      return(d)
    }else{
      return(NULL)
    }
  })
Dat

sort(table(Dat$sample))

Dat <- Dat %>%
  transmute(sampleID = "dummySample", chr = chr, pos = pos,
            ref = ref, mut = mut)
Dat

res <- dndscv(mutations = Dat, refdb = "refdb/Fplautii.rda",
              max_coding_muts_per_sample = Inf,
              max_muts_per_gene_per_sample = Inf, numcode = 1)
