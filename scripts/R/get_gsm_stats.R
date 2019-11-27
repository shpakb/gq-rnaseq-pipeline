#########################################################################
# Takes list of GSM paths
# Returns gsm_stats table in data.frame
#########################################################################
library(tidyverse)

args <- commandArgs(TRUE)

outFile <- args[1]
gsms_f <- args[2:length(args)]

cat(sprintf("Out file: %s \n", outFile))
cat(sprintf("Aggregating statistics on: %i GSM files... \n", length(gsms_f)))

gsm_stats <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(gsm_stats) <-  c("GSM", "N_GENES_EXP")

for (gsm_f in gsms_f) {
    gsm <- read.csv(gsm_f,sep="\t")
    gsm_id <- gsm_f %>% str_extract("GSM\\d+")
    n_genes_exp <- sum(gsm$est_counts==0)
    gsm_stats[nrow(gsm_stats)+1,] <- c(gsm_id, n_genes_exp)
}

print(gsm_stats)

write.table(gsm_stats, outFile, sep='\t', quote=F, row.names=F)

cat("Done. \n")