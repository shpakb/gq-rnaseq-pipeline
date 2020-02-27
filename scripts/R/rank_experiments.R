# Takes results of some statistical test with pvals.
# Outputs tranketed list(1000 entries) of results ranked by pval.

args <- commandArgs(TRUE)

pc_list_file <- args[1]
geneset_file <- args[2]
output_file <- args[3]
test <- args[4]

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
    INTERSECTION=0,
    stringsAsFactors=F
  )

if (test %in% c("ks", "wilcoxGST")){
  suppressMessages(library(tidyverse))
  suppressMessages(library(data.table))
  suppressMessages(library(limma))

  for (pc_name in names(pc_list)){
    tryCatch({
      print(pc_name)
      gst_result <- list()
      pc <- pc_list[[pc_name]]
      indexes <- which(names(pc) %in% geneset)
      pval <- limma::wilcoxGST(index=indexes, statistics=pc, alternative = "either")
      intersection <- length(indexes)
      output_df <- rbind(output_df, c(pc_name, pval, intersection))
   }, error = function(e) {
      print("woops!")
      print(e$message)
    })
  }

  output_df <- output_df[2:nrow(output_df),]
}


if (test %in% c("fgsea", "fgsea_multilevel")){
  suppressMessages(library(tidyverse))
  suppressMessages(library(fgsea))
  suppressMessages(library(data.table))
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))

  print("Parralel execution...")

  cores=detectCores()
  cl <- makeCluster(8)
  registerDoParallel(cl)

  output_df <- foreach(i=1:length(pc_list), .combine=rbind, .packages=c("tidyverse", "fgsea")) %dopar% {
    tryCatch({
      pc_name <- names(pc_list)[i]
      pc <- pc_list[[pc_name]]
      fgsea_result <- fgsea::fgsea(pathways = genesets, stats = pc, nperm=500, gseaParam = gsea_param)
      fgsea_result <- fgsea_result %>% select(pval, ES, NES, size)
      fgsea_result$LABEL <- pc_name
      colnames(fgsea_result) <- c("PVAL", "ES", "NES", "INTERSECTION_SIZE", "LABEL")

      return(fgsea_result)

    }, error = function(e) {})
}


print("Writing results")

write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)

