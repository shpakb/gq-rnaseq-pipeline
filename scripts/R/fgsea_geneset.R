suppressMessages(library(tidyverse))
suppressMessages(library(fgsea))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))


args <- commandArgs(TRUE)

pc_list_file <- args[1]
geneset_file <- args[2]
output_file <- args[3]

cat(sprintf("PC list file: %s \n", pc_list_file))
cat(sprintf("Geneset file file: %s \n", geneset_file))
cat(sprintf("Output file: %s \n", output_file))

pc_list <- readRDS(pc_list_file)
geneset <- readLines(geneset_file)

cat(sprintf("Number of PC: %i \n", length(pc_list)))
cat(sprintf("Number of genes in geneset: %i \n", (length(geneset) - 2)))

# Emulation of geneset list for fgsea
genesets <- list()
genesets[['geneset']] <- geneset

output_df <-
  data.frame(
    LABEL="blank",
    PVAL=0,
    ES=0,
    NES=0,
    INTERSECTION_SIZE=0,
    stringsAsFactors = FALSE
  )

for (pc_name in names(pc_list)){
  tryCatch({
    print(pc_name)
    pc <- pc_list[[pc_name]]
    fgsea_result <- fgsea::fgsea(pathways = genesets, stats = pc, nperm=500)
    fgsea_result <- fgsea_result %>% select(pval, ES, NES, size)
    fgsea_result$LABEL <- pc_name
    colnames(fgsea_result) <- c("PVAL", "ES", "NES", "INTERSECTION_SIZE", "LABEL")
    output_df <- rbind(output_df, fgsea_result)
 }, error = function(e) {
    print("woops!")
    print(e$message)
  })
}

output_df <- output_df[2:nrow(output_df),]

print("Writing results")

write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)
