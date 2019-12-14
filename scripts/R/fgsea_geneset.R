suppressMessages(library(tidyverse))
suppressMessages(library(fgsea))
suppressMessages(library(data.table))

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
    PADJ=0,
    NES=0,
    SIZE=0,
    stringsAsFactors = FALSE
  )

for (pc_name in names(pc_list)){
  tryCatch({
    pc <- pc_list[[pc_name]]
    n_genes <- length(pc)
    fgsea_out <- fgsea::fgseaMultilevel(geneset, pc)
    fgsea_out <- fgsea_out %>% select(padj, NES, size) %>% unlist
    output_df <- rbind(output_df, c(pc_name, fgsea_out))
 }, error = function(e) {
    print("woops!")
    print(pc_name)
    print(e$message)
  })
}

output_df <- output_df[2:nrow(output_df),]

print("Writing results")

write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)
