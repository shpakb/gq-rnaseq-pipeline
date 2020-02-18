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
cat(sprintf("Number of genes in geneset: %i \n", (length(geneset) - 2)))

output_df <-
  data.frame(
    LABEL="blank",
    PVAL=0,
    STAT=0,
    INTERSECTION=0
  )

for (pc_name in names(pc_list)){
  tryCatch({
    print(pc_name)
    pc <- pc_list[[pc_name]]
    overlap <- pc[geneset]
    overlap <- overlap[!is.na(names(overlap))]
    ks_result <- ks.test(x = pc, y = overlap, exact = TRUE, alternative = "greater")

    ks_result <- ks_result %>% unlist %>% t %>% as.data.frame %>% select("statistic.D", "p.value")
    ks_result$LABEL <- pc_name
    ks_result$INTERSECTION <- length(overlap)
    colnames(ks_result) <- c("STAT", "PVAL", "LABEL", "INTERSECTION")
    output_df <- rbind(output_df, ks_result)
 }, error = function(e) {
    print("woops!")
    print(e$message)
  })
}

output_df <- output_df[2:nrow(output_df),]

print("Writing results")

write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)
