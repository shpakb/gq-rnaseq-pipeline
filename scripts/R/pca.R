suppressMessages(library(tidyverse))
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

rownames(exp) <- exp$ENTREZ
exp$ENTREZ <- NULL

exp <- as.matrix(exp)

if(logscale_f) {
  exp <- logDataset(exp)
}

# remove genes with 0 variance
filtered_genes_indexes <- which(apply(t(exp), 2, var)!=0)
exp <- exp[filtered_genes_indexes,]

# take to n_genes by max2 expression.
n_genes <- min(nrow(exp), n_genes)
exp <- exp[1:n_genes, ]

pca <- prcomp(t(exp), scale=T, center=T)

saveRDS(pca, file = pca_out_file)