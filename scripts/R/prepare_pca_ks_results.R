suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

ks_results_df_file <- args[1]
gse_df_file <- args[2]
annotated_output_file <- args[3]

cat(sprintf("KS results file: %s \n", ks_results_df_file))
cat(sprintf("GSE df file: %s \n", gse_df_file))
cat(sprintf("Annotated output file: %s \n", annotated_output_file))

ks_results_df <- read.csv(ks_results_df_file, sep='\t', stringsAsFactors=F)
n_tests <- ks_results_df %>% nrow()

print("Sorting results by NES...")
ks_results_df <- ks_results_df[order(ks_results_df$PVAL, decreasing = F),]

print("Cutting top 1000...")
cuttof <- min(nrow(ks_results_df), 1000)
ks_results_df <- ks_results_df[1:cuttof,]

print("Annotating results with experiment titles...")
gse_df <-
  read.csv(gse_df_file, sep = "\t", stringsAsFactors=F) %>%
  select(GSE, TITLE)

gse_df$GSE <- gse_df$GSE %>% str_extract("GSE\\d+")
gse_df <- gse_df %>% distinct()

ks_results_df$GSE <- ks_results_df$LABEL %>% str_extract("GSE\\d+")

ks_results_df <- merge(ks_results_df, gse_df, all.x=T)

ks_results_df$GSE <- NULL

ks_results_df <- ks_results_df[order(ks_results_df$PVAL, decreasing = F),]

ks_results_df$PADJ <- ks_results_df$PVAL %>% p.adjust(p, method = "bonferroni", n = n_tests)

ks_results_df$LOG_PVAL <- log10(ks_results_df$PVAL)
ks_results_df$LOG_PADJ <- log10(ks_results_df$PADJ)

ks_results_df <- ks_results_df %>% select("LABEL", "LOG_PVAL", "LOG_PADJ", "INTERSECTION", "TITLE")

write.table(ks_results_df, annotated_output_file, col.names = T, row.names = F, sep = "\t", quote=F)

print("Done.")
