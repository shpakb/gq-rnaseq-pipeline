suppressMessages(library(tidyverse))
suppressMessages(library(matrixStats))

logDataset <- function(ge) {
  if (is_logscale(ge))
    return(ge)
  return(log2(ge + 1))
}

is_logscale <- function(x) {
  qx <- quantile(as.numeric(x), na.rm = T)
  if (qx[5] - qx[1] > 100 || qx[5] > 100) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

args <- commandArgs(TRUE)

exp_file <- args[1]
pca_out_file <- args[2]
n_genes <- as.integer(args[3])
logscale_f <- (args[4] == "log")

cat(sprintf("Expression matrix file: %s \n", exp_file))
cat(sprintf("PCA out file: %s \n", pca_out_file))
cat(sprintf("Number of genes considered: %s \n", n_genes))
cat(sprintf("Logscale: %s \n", as.character(logscale_f)))

exp <- read.csv(exp_file, sep="\t")

rownames(exp) <- exp$entrez
exp$entrez <- NULL

exp <- as.matrix(exp)

if(logscale_f) {
  exp <- logDataset(exp)
}

vars <- rowVars(exp)
exp <- exp[which(vars!=0),]

print("Cutting genes with largest VMR...")
vmr <- rowVars(exp)/rowMeans(exp)
exp <- exp[order(vmr, decreasing = T),]
n_genes <- min(nrow(exp), n_genes)
exp <- exp[1:n_genes, ]

print("Doing PCA...")
pca <- prcomp(t(exp), scale=T, center=T)

saveRDS(pca, file = pca_out_file)
print("Done.")