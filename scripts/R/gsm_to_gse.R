suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

gsmFile <- args[1]
gse_id <- gsmFile %>% str_extract("GSE\\d+")

gsm_files <-
    args[4:length(args)] %>%
    read.csv(sep="\t") %>%
    filter(GSE==gse_id) %>%
    select(GSM) %>%
    unlist() %>%
    unique %>%
    paste(gsm_folder, "/", ., ".tsv", sep="")

cat(sprintf("Working on: %s \n", gse_id))
cat(sprintf("GSM folder: %s \n", gsm_folder))
cat(sprintf("Gene annotation: %s \n", args[3]))
cat(sprintf("Aggregating %i GSM \n", length(gsm_files)))

print(args)
print(gsm_files)

# getting Ensamble to Entrez mapping
# removing all genes without Entrez annotation
geneAnnot <-
  args[3] %>%
  read.csv(stringsAsFactors = F,
           sep = "\t",
           header = F) %>%
  filter(V3!="NONE")

colnames(geneAnnot) <- c("GENE", "SYMBOL", "ENTREZ")


#############################FUNCTIONS############################

max2 <- function(array) {
  n = length(array)
  sort(array, partial = n - 1)[n - 1]
}

# annotates GSE with Entrez gene annotation
# aggregates Entrez probes with same name by max value
annotate_genes <- function(gse, geneAnnot) {
  gse <-
    merge(geneAnnot, gse, by="GENE") %>%
    select(-c("GENE", "SYMBOL"))
  # select maximum expressing probes
  gse <-
    gse %>%
    aggregate(.~ENTREZ, ., sum)
  return(gse)
}

# aggregates GSE by "est_counts" column
aggregate_gse <- function(gsm_files, geneAnnot) {
  gse_tpm <-
    gsm_files[1] %>%
    read.csv(sep = "\t")

  # select(gene)
  gse_cpm <- gse_tpm

  for (gsm_file in gsm_files) {
    gsm_id <- gsm_file %>% str_extract("GSM\\d+")
      print(gsm_id)
    gsm <-
      gsm_file %>%
      read.csv(sep = "\t",
               stringsAsFactors = F)
    # from "est_counts" to cpm
    gsm$cpm <- (gsm$est_counts/sum(gsm$est_counts))*10^6
    gse_tpm <-
      gse_tpm %>%
      cbind(gsm$tpm)
    colnames(gse_tpm)[ncol(gse_tpm)] <- gsm_id
    gse_cpm <-
      gse_cpm %>%
      cbind(gsm$cpm)
    colnames(gse_cpm)[ncol(gse_cpm)] <- gsm_id
  }
  gse_tpm <- annotate_genes(gse_tpm, geneAnnot)
  rownames(gse_tpm) <- gse_tpm$ENTREZ
  gse_tpm$entrez <- NULL

  gse_cpm <- annotate_genes(gse_cpm, geneAnnot)
  rownames(gse_cpm) <- gse_cpm$entrez
  gse_cpm$entrez <- NULL

  #get genes sorted
  gse_tpm$max2 <- apply(gse_tpm, 1, FUN = max2)

  gse_tpm <-
    gse_tpm[order(gse_tpm$max2, decreasing = T),]

  gse_tpm$max2 <- NULL

  gse_cpm <- gse_cpm[rownames(gse_tpm),]

  gse_cpm$entrez <- rownames(gse_cpm)

  return(gse_cpm)
}

##################EXECUTION##########################

# getting filtered list of gsm for GSE
# filtering list by QC

cat("Aggregating... \n")

gse <- aggregate_gse(gsm_files, geneAnnot)

# rounding for memory efficiency
gse[, unlist(lapply(gse, is.numeric))] <-
  round(gse[, unlist(lapply(gse, is.numeric))], 3)

# write GSE
write.table(x=gse,
            outFile,
            sep = "\t",
            row.names = F,
            quote = F
            )

cat("Done.\n")

