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
  gene_anotation_file %>%
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
  gse <-
    gsm_files[1] %>%
    read.csv(sep = "\t")

  for (gsm_file in gsm_files) {
    gsm_id <- gsm_file %>% str_extract("GSM\\d+")
    gsm <-
      gsm_file %>%
      read.csv(sep = "\t",
               stringsAsFactors = F)

    gse <-
      gse %>%
      cbind(gsm$EST_COUNTS)
    colnames(gse)[ncol(gse)] <- gsm_id
  }

  print("Anotating genes...")
  gse <- annotate_genes(gse, geneAnnot)

  gse <- gse %>% select(-c(EST_COUNTS, TPM))

  return(gse)
}

##################EXECUTION##########################

# getting filtered list of gsm for GSE
# filtering list by QC
print("Aggregating exp table...")
gse <- aggregate_gse(gsm_files, geneAnnot)

print("Writing expression table...")


write.table(x=gse,
            gse_file,
            sep = "\t",
            row.names = F,
            quote = F
            )

cat("Done.\n")

