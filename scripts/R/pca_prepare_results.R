suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

gsea_results_df_file <- args[1]
gse_df_file <- args[2]
annotated_output_file <- args[3]
gsea_normalization_map_file <- args[4]

cat(sprintf("GSEA results file: %s \n", gsea_results_df_file))
cat(sprintf("GSE df file: %s \n", gse_df_file))
cat(sprintf("Annotated output file: %s \n", annotated_output_file))
cat(sprintf("GSEA ES normalization map file: %s \n", gsea_normalization_map_file))

gsea_normalization_map <- read.csv(gsea_normalization_map_file, sep="\t")
gsea_results_df <- read.csv(gsea_results_df_file, sep='\t', stringsAsFactors=F)
gsea_results_df$NES <- 0

print("Normalizing enrichment score...")
for (i in 1:nrow(gsea_results_df)){
  intersection_size <- gsea_results_df[i,3]
  es <- gsea_results_df[i,2]
  mean_es <-
    gsea_normalization_map %>%
    filter(INTERSECTION_SIZE==intersection_size) %>%
    select(MEAN_ES) %>% unlist %>% unname
  nes <- es/mean_es
  gsea_results_df[i,4] <- nes
}

print("Sorting results by NES...")
gsea_results_df <- gsea_results_df[order(abs(gsea_results_df$NES), decreasing = T),]

print("Cutting top 1000...")
cuttof <- min(nrow(gsea_results_df), 1000)
gsea_results_df <- gsea_results_df[1:cuttof,]

print("Annotating results with experiment titles...")
gse_df <-
  read.csv(gse_df_file, sep = "\t", stringsAsFactors=F) %>%
  select(GSE, TITLE)

gse_df$GSE <- gse_df$GSE %>% str_extract("GSE\\d+")
gse_df <- gse_df %>% distinct()

gsea_results_df$GSE <- gsea_results_df$PC_NAME %>% str_extract("GSE\\d+")

gsea_results_df <- merge(gsea_results_df, gse_df, all.x=T)

gsea_results_df$GSE <- NULL

write.table(gsea_results_df, annotated_output_file, col.names = T, row.names = F, sep = "\t", quote=F)

print("Done.")