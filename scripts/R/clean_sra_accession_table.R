suppressMessages(library(tidyverse))

args <- commandArgs(TRUE)

print("Some extra cleaning...")

srr_df <-
    args[1] %>%
    read.csv(sep=" ")

srr_df$GSM <-
    srr_df$GSM %>%
    str_extract("GSM\\d+")

write.table(srr_df, args[1], col.names = T, row.names = F, sep = "\t", quote=F)

print("Done.")