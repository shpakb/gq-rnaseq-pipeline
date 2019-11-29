#########################################################################
# Takes gsm stats df and GSM to GSE df and filters for GSM QC
# Writes gsm_gse_filtered.tsv with all GSM passed QC and their mappings
# to GSE
########################################################################
suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)
min_n_gsm <- as.integer(args[1])
min_exp_genes <- as.integer(args[2])
prequant_gsm_gse_df_file <- args[3]

# Outs
gsm_stats_df_file <- args[4]
postquant_gsm_gse_df_file <- args[5]
passing_gse_list_file <- args[6]

# Input files
gsm_files <- args[7:length(args)]

cat(sprintf("Aggregating statistics on: %i", length(gsm_files)))
cat(sprintf("Min number of GSM per GSE: %i \n", min_n_gsm))
cat(sprintf("Min number of expressed genes per GSM: %i \n", min_exp_genes))
cat(sprintf("GSM GSE df file: %s \n", prequant_gsm_gse_df_file))
cat(sprintf("GSM stats out file: %s \n", gsm_stats_df_file))
cat(sprintf("GSM GSE df out file: %s \n", postquant_gsm_gse_df_file))
cat(sprintf("Passing gse list file: %s \n", passing_gse_list_file))


prequant_gsm_gse_df <- read.csv(prequant_gsm_gse_df_file, sep="\t")

gsm_stats_df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(gsm_stats_df) <-  c("GSM", "N_GENES_EXP")

for (f in gsm_files) {
    gsm <- read.csv(f,sep="\t")
    gsm_id <- f %>% str_extract("GSM\\d+")
    n_genes_exp <- sum(gsm$est_counts==0)
    gsm_stats_df[nrow(gsm_stats_df)+1,] <- c(gsm_id, n_genes_exp)
}

write.table(gsm_stats_df, gsm_stats_df_file, sep='\t', quote=F, row.names=F)


good_gsm <-
    gsm_stats_df %>%
    filter(N_GENES_EXP >= min_exp_genes) %>%
    select(GSM) %>%
    unlist()

postquant_gsm_gse_df <-
    prequant_gsm_gse_df %>%
    filter(GSE %in% good_gse) %>%
    filter(GSM %in% good_gsm)

write.table(postquant_gsm_gse_df, postquant_gsm_gse_df_file, sep="\t", row.names=F, quote=F)


passing_gse_list <-
    prequant_gsm_gse_df %>%
    filter(GSM %in% good_gsm) %>%
    count(GSE) %>%
    filter(n >= min_n_gsm) %>%
    select("GSE") %>%
    unlist()

writeLines(passing_gse_list, passing_gse_list_file)

cat("Done. \n")