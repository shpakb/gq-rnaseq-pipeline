suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)


gene_anotation_file <- args[1]
gse_file <- args[2]
gse_id <- gse_file %>% str_extract("GSE\\d+")
gsm_files <- args[3:length(args)]

cat(sprintf("Working on: %s \n", gse_id))
cat(sprintf("Output file: %s \n", gse_file))
cat(sprintf("Gene annotation: %s \n", gene_anotation_file))
cat(sprintf("Aggregating %i GSM \n", length(gsm_files)))

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
  count <- 0
  for (gsm_file in gsm_files) {
    gsm_id <- gsm_file %>% str_extract("GSM\\d+")
    gsm <-
      gsm_file %>%
      read.csv(sep = "\t",
               stringsAsFactors = F)
    # from "est_counts" to cpm
    gsm$CPM <- (gsm$EST_COUNTS/sum(gsm$EST_COUNTS))*10^6

    gse_tpm <-
      gse_tpm %>%
      cbind(gsm$TPM)
    colnames(gse_tpm)[ncol(gse_tpm)] <- gsm_id

    gse_cpm <-
      gse_cpm %>%
      cbind(gsm$CPM)
    colnames(gse_cpm)[ncol(gse_cpm)] <- gsm_id
  }

  print("Annotating tpm table...")
  print(head(gsm_tpm))
  gse_tpm <- annotate_genes(gse_tpm, geneAnnot)
  rownames(gse_tpm) <- gse_tpm$ENTREZ
  gse_tpm$ENTREZ <- NULL
  print("Anotation cpm table...")
  gse_cpm <- annotate_genes(gse_cpm, geneAnnot)
  rownames(gse_cpm) <- gse_cpm$ENTREZ
  gse_cpm$ENTREZ <- NULL

  print("Sorting genes according to tpm...")
  gse_tpm$MAX2 <- apply(gse_tpm, 1, FUN = max2)

  gse_tpm <-
    gse_tpm[order(gse_tpm$MAX2, decreasing = T),]

  gse_tpm$MAX2 <- NULL

  gse_cpm <- gse_cpm[rownames(gse_tpm),]

  gse_cpm$ENTREZ <- rownames(gse_cpm)

  return(gse_cpm)
}

##################EXECUTION##########################

# getting filtered list of gsm for GSE
# filtering list by QC
print("Aggregating exp table...")
gse <- aggregate_gse(gsm_files, geneAnnot)

print("Rounding expression table...")
gse[, unlist(lapply(gse, is.numeric))] <-
  round(gse[, unlist(lapply(gse, is.numeric))], 3)

print("Writing expression table...")
write.table(x=gse,
            gse_file,
            sep = "\t",
            row.names = F,
            quote = F
            )

cat("Done.\n")

