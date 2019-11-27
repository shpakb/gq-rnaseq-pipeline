args <- commandArgs(T)
gse_table_file <- args[1]
gsm_table_file <- args[2]
gpl_table_file <- args[3]

#filters
organism <- args[3]
min_n_spots <- args[4]
max_n_spots <- args[4]
min_n_gsm <- args[1]
max_n_gsm <- args[1]
srr_df_file <- args[1]

library(tidyverse)

gsm_df <- read.csv(gsm_table_file, sep = "\t", header = T, stringsAsFactors = F)
gsm_df$gse <- gsm_df$gse %>% str_extract("GSE\\d+")
gsm_df[is.na(gsm_df)] <- "NA"
gsm_df <- gsm_df %>% distinct()
gsm_df <- gsm_df %>% select(c(gsm, gse, organism, gpl, gsm_library_selection, gsm_library_strategy))
gsm_df$lab <- organism

gse_df <- read.csv(gse_table_file, sep = "\t", stringsAsFactors = F)
gse_df <- gse_df %>% select(c(gse, super_series))
gse_df$lab <- organism

gpl_df <- read.csv(gpl_table_file, sep="\t")

srr_df <- read.csv(srr_df_file, sep=" ", stringsAsFactors = F, header = F)
srr_df$srr_df <- srr_df$V2 %>% str_extract("GSM\\d+")
srr_df <- srr_df %>% filter(gsm %in% unique(gsm_df$gsm))
srr_df$spots <- srr_df$spots %>% as.integer()
srr_df <- aggregate(spots~gsm, srr_df, sum)

########################################################################################################################
df <-
  gsm_df %>%
  select(gsm, organism) %>%
  distinct()

cDNA <-
  gsm_df %>%
  filter(gsm_library_selection == "cDNA") %>%
  select(gsm) %>%
  unlist %>%
  unique

right_organism <-
  gsm_df %>%
  filter(lab==organism) %>%
  select(gsm) %>%
  unlist %>%
  unique

rna_seq <-
  gsm_df %>%
  filter(gsm_library_strategy == "RNA-Seq") %>%
  select(gsm) %>%
  unlist %>%
  unique

gpl_list <-
  gpl_df %>%
  filter(is_illumina) %>%
  select(gpl) %>%
  unlist %>%
  unique

illumina <-
  gsm_df %>%
  filter(gpl %in% gpl_list) %>%
  select(gsm) %>%
  unlist %>%
  unique

super_series_list <-
  gse_df %>%
  filter(super_series) %>%
  select(gse) %>%
  unlist %>%
  unique

has_srr <- srr_df$V2 %>% unique

spots_filtered <-
  df1 %>%
  filter(spots < max_n_spots) %>%
  filter(spots > min_n_spots) %>%
  select(gsm) %>%
  unlist %>%
  unique

# getting exclusive super super_series gsm
a <- gsm_df %>% filter(gse %in% super_series_list) %>% select(gsm) %>% unlist %>% unique()
b <- gsm_df %>% filter(!(gse %in% super_series_list)) %>% select(gsm) %>% unlist %>% unique()
super_series_gsm <- setdiff(a, b)

df2 <-
  data.frame(
    gsm = df$gsm,
    organism = df$organism,
    not_super_series = !(df$gsm %in% super_series_gsm),
    cDNA = df$gsm %in% cDNA,
    right_organism = df$gsm %in% right_organism,
    rna_seq = df$gsm %in% rna_seq,
    illumina = df$gsm %in% illumina,
    has_srr = df$gsm %in% has_srr,
    spots_cutoff = df$gsm %in% spots_filtered)

write.table(df2, , col.names = T, row.names = F, sep = "\t")


