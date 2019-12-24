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
    PADJ=0,
    NES=0,
    SIZE=0,
    stringsAsFactors = FALSE
  )

# cores <- detectCores()
cl <- makeCluster(10)
registerDoParallel(cl)

# cat(sprintf("Number of cores: %i \n", cores[1]))
output_df <- foreach(pc_name=names(pc_list),
                     .combine=rbind,
                     .export = ls(globalenv()),
                     .packages = c("tidyverse", "fgsea")) %dopar% {
      pc <- pc_list[[pc_name]]
      n_genes <- length(pc)
      fgsea_out <- fgsea::fgseaMultilevel(geneset, pc)
      fgsea_out <- fgsea_out %>% select(padj, NES, size)
      fgsea_out <- cbind(pc_name, fgsea_out)
      colnames(fgsea_out) <- c("PC_NAME", "PADJ", "NES", "INTERSECTION_SIZE")
      # cat(sprintf("%i \n", count),
      #     file="/gscmnt/gc2676/martyomov_lab/shpakb/gq-rnaseq-pipeline/log.txt",
      #     append=TRUE)

      fgsea_out
}

print("FGSEA complete...")
# for(pc_name in names(pc_list)){
#       print(pc_name)
#       pc <- pc_list[[pc_name]]
#
#       n_genes <- length(pc)
#       fgsea_out <- fgsea::fgseaMultilevel(geneset, pc)
#       fgsea_out <- fgsea_out %>% select(padj, NES, size) %>% unlist %>% unname
#
#       fgsea_out <- c(pc_name, fgsea_out)
#
#       output_df <- rbind(output_df, fgsea_out)
# }

stopCluster(cl)

# print(output_df)
# fgsea_out <- as.data.frame(fgsea_out)
print(output_df)

print("Writing results...")

# output_df <- output_df[2:nrow(output_df),]

write.table(output_df, output_file, col.names = T, row.names = F, sep = "\t", quote=F)
