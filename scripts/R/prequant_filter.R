# TODO output two tables: sra_gsm and gsm_spots
suppressMessages(library(tidyverse))

args <- commandArgs(T)
gse_table_file <- args[1]
gsm_table_file <- args[2]
gpl_table_file <- args[3]
srr_table_file <- args[4]
organism <- args[5] %>% gsub("_", " ", .)
min_n_spots <- as.integer(args[6])
max_n_spots <- as.integer(args[7])
min_n_gsm <- as.integer(args[8])
max_n_gsm <- as.integer(args[9])
gsm_filtering_file <- args[10]
passing_gsm_list_file <- args[11]
srr_gsm_talbe_file <- args[12]
gsm_gse_df_file <- args[13]

cat(sprintf("GSE table file: %s \n", gse_table_file))
cat(sprintf("GSM table file: %s \n", gsm_table_file))
cat(sprintf("GPL table file: %s \n", gpl_table_file))
cat(sprintf("SRR table file: %s \n", srr_table_file))
cat(sprintf("Organism: %s \n", organism))
cat(sprintf("Min number of spots per GSM: %i \n", min_n_spots))
cat(sprintf("Max number of spots per GSM: %i \n", max_n_spots))
cat(sprintf("Min number of GSM per GSE: %i \n", min_n_gsm))
cat(sprintf("Max number of GSM per GSE: %i \n", max_n_gsm))
cat(sprintf("GSM filtering table output file: %s \n", gsm_filtering_file))
cat(sprintf("Passing gsm list: %s \n", passing_gsm_list_file))
cat(sprintf("SRR GSM table file: %s \n", srr_gsm_talbe_file))
cat(sprintf("GSM GSE df file: %s \n", gsm_gse_df_file))


########################################################################################################################

gsm_df <- read.csv(gsm_table_file, sep = "\t", header = T, stringsAsFactors = F)
gsm_df$GSE <- gsm_df$GSE %>% str_extract("GSE\\d+")
gsm_df[is.na(gsm_df)] <- "NA"
gsm_df <- gsm_df %>% distinct()
gsm_df <- gsm_df %>% select(c(GSM, GSE, ORGANISM, GPL, LIBRARY_SELECTION, LIBRARY_STRATEGY))
gsm_df$LAB <- organism

gse_df <- read.csv(gse_table_file, sep = "\t", stringsAsFactors = F)
gse_df <- gse_df %>% select(c(GSE, IS_SUPER_SERIES))

gpl_df <- read.csv(gpl_table_file, sep="\t")

srr_df <- read.csv(srr_table_file, sep="\t", stringsAsFactors = F)
srr_df <- srr_df %>% filter(GSM %in% unique(gsm_df$GSM))
srr_df$SPOTS <- srr_df$SPOTS %>% as.integer()
gsm_spots_df <- aggregate(SPOTS~GSM, srr_df, sum)

########################################################################################################################

df <-
  gsm_df %>%
  select(GSM, GSE, ORGANISM) %>%
  distinct()

cDNA <-
  gsm_df %>%
  filter(LIBRARY_SELECTION == "cDNA") %>%
  select(GSM) %>%
  unlist %>%
  unique

rna_seq <-
  gsm_df %>%
  filter(LIBRARY_STRATEGY == "RNA-Seq") %>%
  select(GSM) %>%
  unlist %>%
  unique

gpl_list <-
  gpl_df %>%
  filter(IS_ILLUMINA) %>%
  select(GPL) %>%
  unlist %>%
  unique

illumina <-
  gsm_df %>%
  filter(GPL %in% gpl_list) %>%
  select(GSM) %>%
  unlist %>%
  unique

super_series_list <-
  gse_df %>%
  filter(IS_SUPER_SERIES)   %>%
  select(GSE) %>%
  unlist %>%
  unique

has_srr <- srr_df$GSM %>% unique

spots_filtered <-
  gsm_spots_df %>%
  filter(SPOTS < max_n_spots) %>%
  filter(SPOTS > min_n_spots) %>%
  select(GSM) %>%
  unlist %>%
  unique

# getting exclusive super super_series gsm
a <- gsm_df %>% filter(GSE %in% super_series_list) %>% select(GSM) %>% unlist %>% unique()
b <- gsm_df %>% filter(!(GSE %in% super_series_list)) %>% select(GSM) %>% unlist %>% unique()
super_series_gsm <- setdiff(a, b)

df1 <-
  data.frame(
    GSM = df$GSM,
    GSE = df$GSE,
    IS_RIGHT_ORGANISM = df$ORGANISM == organism,
    NOT_SUPER_SERIES = !(df$GSM %in% super_series_gsm),
    CDNA = df$GSM %in% cDNA,
    RNA_SEQ = df$GSM %in% rna_seq,
    IS_ILLUMINA = df$GSM %in% illumina,
    HAS_SRR = df$GSM %in% has_srr,
    SPOTS_CUTOFF = df$GSM %in% spots_filtered,
    stringsAsFactors=F)

df2 <-
  df1 %>%
  filter(IS_RIGHT_ORGANISM & NOT_SUPER_SERIES & CDNA & RNA_SEQ & IS_ILLUMINA & HAS_SRR & SPOTS_CUTOFF) %>%
  count(GSE)


passing_gse <-
  df2 %>%
  filter(n > min_n_gsm) %>%
  filter(n < max_n_gsm) %>%
  select(GSE) %>%
  unlist %>%
  unique

df1$GSE_IS_PASSING <-
  df1$GSE %in% passing_gse

write.table(df1, gsm_filtering_file, col.names = T, row.names = F, sep = "\t", quote=F)

passing_gsm <- unique(df1$GSM)

writeLines(passing_gsm, passing_gsm_list_file)

srr_gsm_df <-
  srr_df %>%
  filter(GSM %in% passing_gsm) %>%
  select(SRR, GSM)

write.table(srr_gsm_df, srr_gsm_talbe_file, col.names = T, row.names = F, sep = "\t", quote=F)

gsm_gse_df <-
  gsm_df %>%
  filter(GSE %in% passing_gse) %>%
  filter(GSM %in% passing_gsm)

write.table(gsm_gse_df, gsm_gse_df_file, col.names = T, row.names = F, sep = "\t", quote=F)