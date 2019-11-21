#########################################################################
# Takes list of GSM
# Returns GSE with gene expression if GQ passed. Empty file otherwise
#########################################################################
library(tidyverse)

args <- commandArgs(TRUE)

gse_file <-
  as.character(args[1])
gse_id <-
  gse_file %>%
  str_extract("GSE\\d+")

cat(sprintf("Working on: %s \n", gse_id))
cat(sprintf("Gene annotation: %s \n", as.character(args[2])))
cat(sprintf("GSM to GSE mapping: %s \n", as.character(args[3])))
cat(sprintf("GSM input dir: %s \n", as.character(args[4])))
cat(sprintf("GC \n Minimum expressd genes per GSM: %s \n", as.character(args[5])))
cat(sprintf("Mimimum GSM per GSE: %s \n", as.character(args[6])))

# gets Ensamble to Entrez mapping
# removes all genes without Entrez annotation
geneAnnot <-
  as.character(args[2]) %>%
  read.csv(stringsAsFactors = F,
           sep = "\t",
           header = F) %>%
  filter(V3!="NONE")

colnames(geneAnnot) <- c("gene", "symbol", "entrez")

gsm_to_gse <-
  as.character(args[3]) %>%
  read.csv(stringsAsFactors=F,
           sep="\t")

inDir <- as.character(args[4])

# QC:
min_exp_genes <- as.integer(args[4])
min_gsm <- as.integer(args[5])

#############################FUNCTIONS############################

# takes list of gsm for gse and gc params
# returns filtered list of GSM
# TRUE if QC passed by GSM, FALSE otherwise
gsm_list_qc <- function(gsm_id_list, inDir, min_exp_genes) {
  filter_vec <- c()
  for (gsm_id in gsm_id_list){
    gsm <-
      gsm_id %>%
      paste(inDir, "/", ., ".tsv", sep="\t") %>%
      read.csv(sep = "\t")
    filter_vec <- c(filter_vec, sum(gsm_id$est_counts == 0) >= min_exp_genes)
  }
  return(gsm_id_list[filter_vec])
}

max2 <- function(array) {
  n = length(array)
  sort(array, partial = n - 1)[n - 1]
}

# annotates GSE with Entrez gene annotation
# aggregates Entrez probes with same name by max value
annotate_genes <- function(gse, geneAnnot) {
  gse <- 
    merge(geneAnnot, gse, by="gene") %>%
    select(-c("gene", "symbol"))
  # select maximum expressing probes
  gse <- 
    gse %>% 
    aggregate(.~entrez, ., sum) 
  return(gse)
}

# takes gsm_list filtered by qc
#
# aggregates GSE by "est_counts" column 
aggregate_gse <- function(gsm_list, inDir, geneAnnot) {
  gse_tpm <-
    gsm_list[1] %>%
    paste(inDir, "/", ., ".tsv") %>%
    read.csv(sep = "\t") %>%
    select(gene)
  
  gse_cpm <- gse_tpm
  
  for (gsm_id in gsm_list) {

    gsm <-
      gsm_id %>%
      paste(inDir, "/", ., ".tsv") %>%
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
  rownames(gse_tpm) <- gse_tpm$entrez
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
gsm_id_list <-
  gsm_to_gse %>%
  read.table(sep="\t") %>%
  filter(gse==gse_id) %>%
  select(gsm) %>%
  unlist %>%
  as.character() %>%
  unique %>%
  gsm_list_qc(inDir, min_exp_genes)

# creating either empty GSE or aggregated GSE
if (length(gsm_list) < min_gsm){
  file.create(gse_file)
} else {
  gse <- aggregate_gse(gsm_id_list, geneAnnot)

  gse[, unlist(lapply(gse, is.numeric))] <-
    round(gse[, unlist(lapply(gse, is.numeric))], 3)

  # write GSE
  write.table(x=gse,
              gse_file,
              sep = "\t",
              row.names = F,
              quote = F
              )
}
