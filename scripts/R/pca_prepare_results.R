suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

gsea_results_file <- args[1]
gse_df_file <- args[2]

cat(sprintf("GSEA results file: %s \n", gsea_results_file))
cat(sprintf("GSE df file: %s \n", gse_df_file))

gse_df <-
  read.csv(gse_df_file, sep = "\t", stringsAsFactors=F) %>%
  select(GSE, TITLE)
gse_df$GSE <- gse_df$GSE %>% str_extract("GSE\\d+")
gse_df <- gse_df %>% distinct()

