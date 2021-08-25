
args <- list(gff,
             enogg,
             genome_fna,
             midas_dir,
             enogg_format = "uhgg",
             outdir)

library(tidyverse)
library(dndscv)
library(seqinr)
# data(dataset_simbreast)
# res <- dndscv(mutations)



# Read gene defintions
cat("Reading gene defintions from GFF3 file...\n")
gff <- read_tsv(args$gff,
                col_names = FALSE,
                col_types = cols(X1 = col_character(),
                                 X3 = col_character(),
                                 X4 = col_number(),
                                 X5 = col_number(),
                                 X7 = col_character(),
                                 X9 = col_character()))
cat("\tExtracting gene IDs...\n")
gff$gene.id <- gff$X9 %>%
  map_chr(function(x){
    id <- str_split(x, pattern = ";")[[1]][1]
    id <- str_remove(id, "^ID=")
    id
  })
# gff

# Annotations provide gene names
cat("Reading eggNOG annotations...\n")
annot <- HMVAR::read_eggnog(args$enogg,
                            format = args$enogg_format)

# Create cds object for cds file
cat("Creating CDS file...\n")
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

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

cds_tsv <- file.path(args$outdir, "cdsfile.tsv")
write_tsv(cds, cds_tsv)
# name2id <- cds %>%
#   select(gene_name = gene.name, gene_id = gene.id)

### Build new reference file so I can include all sites
cat("Preparing reference genome...\n")
# Read midas data
cat("\tReading MIDAS data...\n")
Dat <- HMVAR::read_midas_data(args$midas_dir, cds_only = TRUE)
# site2gene <- Dat$info %>% select(site_id, gene_id)

# Identify positions to change in reference
cat("\tIdentifying sites to change...\n")
to_change <- Dat$info %>%
  filter(ref_allele != major_allele) %>%
  select(ref_id, ref_pos, ref_allele, major_allele, minor_allele) %>%
  mutate(major_allele = tolower(major_allele),
         minor_allele = tolower(minor_allele),
         ref_allele = tolower(ref_allele))

# Read fasta
cat("\tReading original genome fasta...\n")
seq <- read.fasta(args$genome_fna)

# Change positions
cat("\tChanging sites...\n")
for(i in 1:nrow(to_change)){
  seq[[to_change$ref_id[i]]][to_change$ref_pos[i]] <- to_change$major_allele[i]
}

# Write custom ref for dnds
cat("\tWriting new reference")
ref_fna <- file.path(args$outdir, "reference.fna")
write.fasta(seq, names = names(seq), file.out = ref_fna)

# Build reference
cat("Build reference...\n")
buildref(cdsfile = cds_tsv,
         genomefile = ref_fna, 
         outfile = file.path(args$outdir, "reference.rda"),
         numcode = 1, excludechrs = NULL,
         onlychrs = NULL, useids = FALSE)
