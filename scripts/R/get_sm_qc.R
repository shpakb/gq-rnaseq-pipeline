suppressMessages(library(tidyverse))

gse_df_file <- args[1]
sm_qc_df_out_file <- args[2]


cat(sprintf("GSE df input file: %s \n", gse_df_file))


supported_gpl_list <-
  list.files(gpl_input_folder) %>%
  str_extract("GPL\\d+")

cat(qc_values)

sm_qc_df <-
  read.csv(gse_df_file, sep = "\t") %>%
  select(GSE, NUMBER_GSM, GPL, SUPER_SERIES_OF)

sm_qc_df$TAG <- "NA"
sm_qc_df$HAS_EXP_MAT <- "NA"
sm_qc_df$LOGAV <- "NA"
sm_qc_df$LINMAX <- "NA"
sm_qc_df$LOGMAX <- "NA"
sm_qc_df$N_GENES <- "NA"
sm_qc_df$PASSED <- FALSE

write(qc_values, qc_values_out_file)