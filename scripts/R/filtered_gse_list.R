#########################################################################
# Takes gsm stats df and GSM to GSE df and filters for GSM QC
# Writes gsm_gse_filtered.tsv with all GSM passed QC and their mappings
# to GSE
########################################################################
library(tidyverse)

args <- commandArgs(TRUE)

gsm_stats_file <- args[1]
outFile <- args[2]
gsm_gse_file <- args[3]
min_gsm <- as.integer(args[4])
min_exp_genes <- as.integer(args[5])

cat(sprintf("GSM stats file: %s \n", gsm_stats_file))
cat(sprintf("Out file: %s \n", outFile))
cat(sprintf("GSM to GSE file: %s \n", gsm_gse_file))
cat(sprintf("Min GSM: %i \n", min_gsm))
cat(sprintf("Min expressed genes: %s \n", min_exp_genes))

gsm_stats_df <- read.csv(gsm_stats_file, sep="\t")
gsm_gse_df <- read.csv(gsm_gse_file, sep="\t")

good_gsm <-
    gsm_stats_df %>%
    filter(N_GENES_EXP >= min_exp_genes) %>%
    select(GSM) %>%
    unlist() %>%
    as.character()

good_gse <-
    gsm_gse_df %>%
    filter(gsm %in% good_gsm) %>%
    count(gse) %>%
    filter(n >= min_gsm) %>%
    select("gse") %>%
    unlist() %>%
    as.character()

gsm_gse_filtered_df <-
    gsm_gse_df %>%
    filter(gse %in% good_gse) %>%
    filter(gsm %in% good_gsm)

write.table(gsm_gse_filtered_df, outFile, sep="\t", row.names=F, quote=F)


cat("Done. \n")