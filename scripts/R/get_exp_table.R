####################################ARGUMENTS#########################
suppressMessages(library(tidyverse))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))

args <- commandArgs(TRUE)

sm_input_file <- args[1]
gpl_input_folder <- args[2]
exp_mat_out_file <- args[3] # put file in folder with name of GSE. Give file proper name.????
qc_out_file <- args[4]

cat(sprintf("SM input file: %s \n", sm_input_file))
cat(sprintf("GPL annotations input folder: %s \n", gpl_input_folder))
cat(sprintf("Expression matrices output file: %s \n", exp_mat_out_file))
cat(sprintf("QC output file: %s \n", qc_out_file))

####################################UTILS###############################################################################
linearizeDataset <- function(ge) {
  if (is_logscale(ge))
    return(2 ^ ge - 1)
  return(ge)
}

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

max2 <- function(array) {
  n = length(array)
  sort(array, partial = n - 1)[n - 1]
}

eigenorm <- function (val, min, max) {
  if (min == 0 | max == 0)
    stop ("Line minimum or a maximum equal to 0")
  if (val <= 0) {
    return(round(val / abs(min), digits = 3))
  } else {
    return(round(val / max, digits = 3))
  }
}

exvector <- function (v) {
  min = min(v)
  max = max(v)
  for (i in 1:length(v)) {
    if (v[i] <= 0) {
      v[i] = round(v[i] / abs(min), digits = 3)
    } else {
      v[i] = round(v[i] / max, digits = 3)
    }
  }
  return(v)
}

replaceWithNa <- function(x, ...) {
  if (length(x) == 0)
    return(NA)
  return(x)
}

readGPLTable <- function(gpl) {
  df <-
    read.table(
      gzfile(gpl),
      sep = "\t",
      comment.char = "",
      stringsAsFactors = F
    )

  column <- df$V3 %>% str_extract_all("(\\d+)")
  indx <-
    column %>%
    lapply(as.integer) %>%
    lapply(which.min) %>%
    lapply(replaceWithNa) %>%
    unlist

  df$V3 <-
    lapply(1:length(column), function(x)
      column[[x]][indx[x]]) %>% unlist

  column <- strsplit(df$V2, "///")
  df$V2 <-
    lapply(1:length(column), function(x)
      column[[x]][indx[x]]) %>% unlist

  df[is.na(df)] <- "NONE"

  rownames(df) <- df$V1
  colnames(df) <- c("Probe_id", "Symbol", "Entrez_ID")

  df
}

############################################EXECUTION###################################################################
#TAG N_GSM IS_SUPER_SERIES HAS_EXP_MAT GPL LOGAV LOGMAX LINMAX N_GENES
qc_df <-
  data.frame(
  TAG=sm_input_file %>% str_extract("GSE\\d+_GPL\\d+|GSE\\d+"),
  N_GSM="NA",
  GPL="NA",
  HAS_EXP_MAT="NA",
  LOGAV="NA",
  LINMAX="NA",
  LOGMAX="NA",
  N_GENES="NA",
  PROCESSED=FALSE
)

tryCatch(
{
  in.con <- gzfile(sm_input_file)
  wholeFile <- readLines(in.con)
  gplId <- which(grepl("!Sample_platform_id", wholeFile))
  gplId <- gsub(".*(GPL\\d+).*", "\\1", wholeFile[gplId])
  gpl <-  paste0(gpl_input_folder, "/", gplId, ".3col.gz")
  qc_df$GPL <- gplId

  tableStart <-
    which(grepl("!series_matrix_table_begin", wholeFile))
  tableEnd <-
    which(grepl("!series_matrix_table_end", wholeFile))

  if ((tableEnd - tableStart) < 10 ){
    print("No expression table found.")
    qc_df$HAS_EXP_MAT <- FALSE
  } else{
    qc_df$HAS_EXP_MAT <- TRUE
    qc_df$N_GSM <- ncol(exp)
  }

  expressionTable <-
    read.table(
      in.con,
      sep = "\t",
      header = 1,
      row.names = 1,
      comment.char = "",
      skip = tableStart,
      nrows = tableEnd - tableStart - 2,
      stringsAsFactors = F
    )

  expressionTable <- as.matrix(expressionTable)

  exp <- linearizeDataset(expressionTable)
  explog <- logDataset(expressionTable)

  qc_df$LOGAV <- mean(explog) %>% round(1)
  qc_df$LOGMAX <- max(explog) %>% round(1)
  qc_df$LINMAX <- max(exp) %>% round()

  gplTable <- readGPLTable(gpl)

  exp <- cbind(gplTable[rownames(exp),], exp)

  exp <-
    exp[!(row.names(exp) %in% c("NONE", "NA", "none", "NULL")),]

  exp <-
    collapseRows(exp[4:ncol(exp)], exp$Entrez_ID, rownames(exp))

  exp <- exp$datETcollapsed

  exp <-
    exp[!(row.names(exp) %in% c("NONE", "NA", "none", "NULL")), ]

  exp <- as.data.frame(exp)

  qc_df$N_GENES <- nrow(exp)

  rownames(exp) <- exp$Entrez_ID
  exp$Entrez_ID <- NULL

  exp <- exp %>% round(3)

  exp$max2 <- apply(exp, 1, FUN = max2)
  exp <- exp[order(exp$max2, decreasing = T), ]
  exp$max2 <- NULL

  exp$ENTREZ <- rownames(exp)

  write.table(
    exp,
    file = exp_mat_out_file,
    quote = F,
    col.names = T,
    row.names = T,
    sep = "\t"
  )

  qc_df$PROCESSED <- TRUE

}, error = function(e)
    {
      print(e$message)
      write("", exp_mat_out_file)
    }
)

write.table(
  qc_df,
  file = qc_out_file,
  quote = F,
  col.names = T,
  row.names = F,
  sep = "\t"
    )
