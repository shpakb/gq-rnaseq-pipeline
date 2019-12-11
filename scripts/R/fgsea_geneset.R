suppressMessages(library(tidyverse))
suppressMessages(library(fgsea))

args <- commandArgs(TRUE)

pc_list_file <- args[1]
geneset_file <- args[2]
output_file <- args[3]

cat(sprintf("PC list file: %s \n", pc_list_file))
cat(sprintf("Geneset file file: %s \n", geneset_file))
cat(sprintf("Output file: %s \n", output_file))

pc_list <- readRDS(pc_list_file)
geneset <- readLines(geneset_file)


# removing artifact lines from GQ
geneset <- geneset[3:length(geneset)]

output_df <-
  data.frame(
    PC_NAME="blank",
    GSEA_STAT=0,
    stringsAsFactors = FALSE
  )


for (pc_name in names(pc_list)){
  pc <- pc_list[[pc_name]]
  gsea_stat <- fgsea::calcGseaStat(pc, na.omit(match(geneset, names(pc))))
  print(pc_name)
  print(gsea_stat)
  output_df <- rbind(output_df, c(pc_name, gsea_stat))
}
print(output_df)



write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)