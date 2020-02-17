suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

gsea_results_df_file <- args[1]
gse_df_file <- args[2]
annotated_output_file <- args[3]

cat(sprintf("GSEA results file: %s \n", gsea_results_df_file))
cat(sprintf("GSE df file: %s \n", gse_df_file))
cat(sprintf("Annotated output file: %s \n", annotated_output_file))

gsea_results_df <- read.csv(gsea_results_df_file, sep='\t', stringsAsFactors=F)

print("Sorting results by NES...")
gsea_results_df <- gsea_results_df[order(gsea_results_df$PVAL),]

print("Cutting top 1000...")
cuttof <- min(nrow(gsea_results_df), 1000)
gsea_results_df <- gsea_results_df[1:cuttof,]

print("Annotating results with experiment titles...")
gse_df <-
  read.csv(gse_df_file, sep = "\t", stringsAsFactors=F) %>%
  select(GSE, TITLE)

gse_df$GSE <- gse_df$GSE %>% str_extract("GSE\\d+")
gse_df <- gse_df %>% distinct()

gsea_results_df$GSE <- gsea_results_df$LABEL %>% str_extract("GSE\\d+")

gsea_results_df <- merge(gsea_results_df, gse_df, all.x=T)

gsea_results_df$GSE <- NULL

gsea_results_df <- gsea_results_df[order(gsea_results_df$PVAL),]

write.table(gsea_results_df, annotated_output_file, col.names = T, row.names = F, sep = "\t", quote=F)

print("Done.")