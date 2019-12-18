suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

cat(sprintf("Output GSM file: %s \n", args[1]))
cat(sprintf("Gene mapping: %s \n", args[2]))
cat(sprintf("Number of SRR files to aggregate: %i \n", args[3:length(args)]))

gsmFile <- args[1]

gsm_id <-
  gsmFile %>%
  str_extract("GSM\\d+")

# getting gene_mapping
gene_mapping <-
  as.character(args[2]) %>%
  read.csv(stringsAsFactors = F,
           sep = "\t",
           header=F)

colnames(gene_mapping) <- c("GENE", "TRANSCRIPT")

srr_list <- args[3:length(args)]

##############FUNCS#################

aggregate_gsm <- function(srr_list) {
  gsm <-
    srr_list[1] %>%
    read.table(sep = "\t",
               header = T,
               stringsAsFactors = F)
  
  # converts nan tpm to 0
  gsm$tpm <- as.double(gsm$tpm)
  gsm[is.na(gsm)] <- 0
  
  gsm <-
    gsm %>%
    select("target_id",	"est_counts", "tpm")
  
  if (length(srr_list) > 1) {
    for (srr_file in srr_list[2:length(srr_list)]) {
      srr <-
        srr_file
        read.table(sep = "\t",
                   header = T,
                   stringsAsFactors = F)
      
      # converts nan tpm to 0
      srr$tpm <- as.double(srr$tpm)
      srr[is.na(srr)] <- 0
      
      gsm$tpm <- gsm$tpm + srr$tpm
      
      gsm$est_counts <- srr$est_counts + gsm$est_counts
    }
  }
  
  return(gsm)
}

collapse_transcripts <- function(gsm, gene_mapping) {
    gsm <-
        merge(gsm,
              gene_mapping,
              by.x = "target_id",
              by.y = "TRANSCRIPT",
              all = T) %>%
        select(-"target_id") %>%
        aggregate(. ~ GENE, ., FUN = sum)

    return(gsm)
}

##################EXECUTION##################


gsm <- aggregate_gsm(srr_list)

gsm <- collapse_transcripts(gsm, gene_mapping)

colnames(gsm) <- c("GENE", "EST_COUNTS", "TPM")

write.table(
  gsm,
  gsmFile,
  row.names = F,
  quote = F,
  sep = "\t"
)
