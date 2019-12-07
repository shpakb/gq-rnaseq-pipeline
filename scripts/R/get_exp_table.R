####################################ARGUMENTS#########################
suppressMessages(library(tidyverse))
suppressMessages(library(WGCNA))
suppressMessages(library(stringr))

args <- commandArgs(TRUE)

gse_df_file <- args[1]
sm_input_folder <- args[2]
gpl_input_folder <- args[3]
exp_mat_out_dir <- args[4]

qc_df_out_file <- args[5]
qc_values_out_file <- args[6]

min_n_gsm <- as.integer(args[7])
max_n_gsm <- as.integer(args[8])

logav_min <- as.double(args[9])
logav_max <- as.double(args[10])
linmax_max <- as.double(args[11])
logmax_max <- as.double(args[12])
min_genes <- as.integer(args[13])

cat(sprintf("GSE df input file: %s \n", gse_df_file))
cat(sprintf("SM input folder: %s \n", sm_input_folder))
cat(sprintf("GPL annotations input folder: %s \n", gpl_input_folder))
cat(sprintf("Expression matrices output dir: %s \n", exp_mat_out_dir))
cat(sprintf("QC values output file: %s \n", qc_df_out_file))
cat(sprintf("QC df out file: %s \n", qc_values_out_file))

dir.create(exp_mat_out_dir)
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
supported_gpl_list <-
  list.files(gpl_input_folder) %>%
  str_extract("GPL\\d+")

qc_values <- sprintf(
"Min GSM: %g
Max GSM: %g
Log average min: %g
Log average max: %g
Linear max max: %g
Log max max: %g
Min number of genes: %i
Supported GPL: %s \n",
  min_n_gsm,
  max_n_gsm,
  logav_min,
  logav_max,
  linmax_max,
  logmax_max,
  min_genes,
  paste(supported_gpl_list, collapse = " ")
)

cat(qc_values)

sm_qc_df <-
  read.csv(gse_df_file, sep = "\t") %>%
  select(GSE, NUMBER_GSM, GPL, SUPER_SERIES_OF)

sm_qc_df$IS_SUPER_SERIES <-
  !is.na(sm_qc_df$SUPER_SERIES_OF)
sm_qc_df$SUPER_SERIES_OF <- NULL

filtered_sm <-
  sm_qc_df %>%
  filter(NUMBER_GSM > min_n_gsm & NUMBER_GSM < max_n_gsm) %>%
  filter(GPL %in% supported_gpl_list) %>%
  filter(!IS_SUPER_SERIES) %>%
  select(GSE) %>%
  unlist %>%
  unique

sm_qc_df$TAG <- "NA"
sm_qc_df$HAS_EXP_MAT <- "NA"
sm_qc_df$LOGAV <- "NA"
sm_qc_df$LINMAX <- "NA"
sm_qc_df$LOGMAX <- "NA"
sm_qc_df$N_GENES <- "NA"
sm_qc_df$PASSED <- FALSE

for(sm in filtered_sm){
  tryCatch(
  {
    print(sm)
    sm_file <- paste0(sm_input_folder, "/", sm, "_series_matrix.txt.gz")
    in.con <- gzfile(sm_file)
    wholeFile <- readLines(in.con)
    gplId <- which(grepl("!Sample_platform_id", wholeFile))
    gplId <- gsub(".*(GPL\\d+).*", "\\1", wholeFile[gplId])
    gpl <-  paste0(gpl_input_folder, "/", gplId, ".3col.gz")
    # it is necessary that tag has gpl
    tag <- sm %>% str_extract("GSE\\d+")
    tag <- paste0(tag, "_", gplId)

    sm_qc_df[sm_qc_df$GSE == sm,]$tag <- tag

    if (!file.exists(gpl))
    {
      cat(sprintf("GPL %s not supported", gplId))
      sm_qc_df[sm_qc_df$GSE == sm,]$GPL <- gplId
      next
    }

    tableStart <-
      which(grepl("!series_matrix_table_begin", wholeFile))
    tableEnd <-
      which(grepl("!series_matrix_table_end", wholeFile))

    if ((tableEnd - tableStart) < 10 ){
      print("No expression table found.")
      next
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

    logav <- mean(explog) %>% round(1)
    logmax <- max(explog) %>% round(1)
    linmax <- max(exp) %>% round()

    sm_qc_df[sm_qc_df$GSE == sm,]$LOGAV <- logav
    sm_qc_df[sm_qc_df$GSE == sm,]$LOGMAX <- logmax
    sm_qc_df[sm_qc_df$GSE == sm,]$LINMAX <- linmax

    if (!(!any(is.na(c(
          logav, logmax, linmax
        ))) &&
        logav >= 4 &&
        logav <= 10 &&
        linmax >= 2000 &&
        logmax < 20))
    {
      cat("Failed sanity check.")
      next
    }

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

    sm_qc_df[sm_qc_df$GSE == sm,]$N_GENES <- nrow(exp)

    if (nrow(exp) < min_genes){
      cat("Failed. Number of genes less than threshold.")
      next
    }

    rownames(exp) <- exp$Entrez_ID
    exp$Entrez_ID <- NULL

    exp$max2 <- apply(exp, 1, FUN = max2)
    exp <- exp[order(exp$max2, decreasing = T), ]
    exp$max2 <- NULL

    exp$ENTREZ <- rownames(EXP)

    write.table(
      exp,
      file = paste0(exp_mat_out_dir, "/", tag, ".tsv"),
      quote = F,
      col.names = T,
      row.names = T,
      sep = "\t"
    )

    sm_qc_df[sm_qc_df$GSE == sm,]$PASSED <- TRUE

  }, error = function(e)
      {
        print(e$message)
      }
  )

}

write(qc_values, qc_values_out_file)


write.table(
      sm_qc_df,
      file = qc_df_out_file,
      quote = F,
      col.names = T,
      row.names = F,
      sep = "\t"
    )
