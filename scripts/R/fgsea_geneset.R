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

# removing artifact lines from GQ
geneset <- geneset[3:length(geneset)] %>% data.table()


output_df <-
  data.frame(
    PC_NAME="blank",
    GSEA_STAT=0,
    INTERSECTION=0,
    stringsAsFactors = FALSE
  )

for (pc_name in names(pc_list)){
  tryCatch({
    print(pc_name)
    pc <- pc_list[[pc_name]]
    n_genes <- length(pc)
    gene_intersection <- na.omit(match(geneset, names(pc)))
    gsea_stat <- fgsea::calcGseaStat(pc, gene_intersection)
    output_df <- rbind(output_df, c(pc_name, gsea_stat, length(gene_intersection)))
 }, error = function(e) {
    print("woops!")
    print(e$message)
  })
}

print("Writing results")

write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)
