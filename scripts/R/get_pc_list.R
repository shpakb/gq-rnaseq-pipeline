suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

pc_list_out_file <- args[1]
max_comp <- args[2]
explained_var_threshold <- args[3]
pca_files <- args[4:length(args)]

cat(sprintf("Output file: %s \n", pc_list_out_file))
cat(sprintf("Max number of componets: %s \n", max_comp))
cat(sprintf("Minimum explained veriance threshold: %s persent \n", explained_var_threshold))
cat(sprintf("Number of PCA to concatinate: %s \n", length(pca_files)))


result <- list()

for (pca_file in pca_files) {
  print(pca_file)
  pca <- readRDS(pca_file)
  tag <- pca_file %>% str_extract("GSE\\d+-GPL\\d+|GSE\\d+")

  # % of explained variance vector for all PC
  pev <- (pca$sdev)^2 / sum(pca$sdev^2) * 100

  # chose how many first PC take from pca.
  n_comp <- min(max_comp, ncol(pca$x), sum(pev > explained_var_threshold))

  componets <-
    as.data.frame(pca$rotation)[,1:n_comp] %>%
    as.list()

  names(componets) <- paste0(tag, "_", names(componets))

  for (n in names(componets)){
    names(componets[[n]]) <- rownames(pca$rotation)
  }

  result <- c(result, componets)
}

saveRDS(result, pc_list_out_file)