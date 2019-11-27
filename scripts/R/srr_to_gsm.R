library(tidyverse)

args <- commandArgs(TRUE)

cat(sprintf("Working on: %s \n", as.character(args[1])))
cat(sprintf("Gene mapping: %s \n", as.character(args[2])))
cat(sprintf("SRR to GSE mapping: %s \n", as.character(args[3])))
cat(sprintf("Kallisto out dir: %s \n", as.character(args[4])))
cat(sprintf("Out dir: %s \n", as.character(args[5])))


gsm_id <-
  as.character(args[1])

# getting gene_mapping
gene_mapping <-
  as.character(args[2]) %>%
  read.csv(stringsAsFactors = F,
           sep = "\t",
           header=F)
colnames(gene_mapping) <- c("gene", "transcript")

srr_to_gsm <-
  as.character(args[3]) %>%
  read.csv(stringsAsFactors = F, sep="\t")

kallistoDir <- as.character(args[4])
cat(kallistoDir)

outDir <- as.character(args[5])

##############FUNCS#################

aggregate_gsm <- function(srr_list, inDir) {
  gsm <-
    paste(inDir, "/", srr_list[1], "/abundance.tsv", sep = "") %>%
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
    for (i in 2:length(srr_list)) {
      srr <-
        paste(inDir, "/", srr_list[1], "/abundance.tsv", sep = "") %>%
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
          by.y = "transcript",
          all = T) %>%
    select(-"target_id") %>%
    aggregate(. ~ gene, ., FUN = sum)
  
  return(gsm)
}

##################EXECUTION##################

srr_list <-
  srr_to_gsm %>%
  filter(gsm == gsm_id) %>%
  select("srr") %>%
  unlist %>%
  unique

gsm <- aggregate_gsm(srr_list, kallistoDir)

gsm <- collapse_transcripts(gsm, gene_mapping)

write.table(
  gsm,
  paste(outDir, "/", gsm_id, ".tsv", sep = ""),
  row.names = F,
  quote = F,
  sep = "\t"
)
