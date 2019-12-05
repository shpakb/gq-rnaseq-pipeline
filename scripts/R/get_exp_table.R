####################################ARGUMENTS#########################
library(WGCNA)
library(stringr)
library(tidyverse)
library(reshape)
library(svglite)
library(limma)

source("/home/shpakb/gq/scripts/processing/utils.R")

args <- commandArgs(TRUE)

inFile <- as.character(args[1])
outDir <- as.character(args[2])
gplDir <- as.character(args[3])
nGenes <- as.integer(args[4])
CT <- as.character(args[5])
NT <- as.character(args[6]) %>% gsub("_", " ", .)
logScaleF <- as.logical(args[7])
quantileF <- as.logical(args[8])

# Outs needed
expMatF <- as.logical(args[9])
modulesF <- as.logical(args[10])
eigengenesF <- as.logical(args[11])
svgF <- as.logical(args[12])

cat(sprintf("Processing: %s\n", inFile))
cat(sprintf(
  "Top %d most expressed genes will be considered\n corType = %s \n network type: %s\n",
  nGenes,
  CT,
  NT
))
cat(
  sprintf(
    "Expression Matrix: %s; Modules: %s; Eigengenes: %s; SVG-plots: %s\n",
    expMatF,
    modulesF,
    eigengenesF,
    svgF
  )
)

############################################EXECUTION###########################
tryCatch({
  in.con <- gzfile(inFile)
  wholeFile <- readLines(in.con)
  processingStatus <- "proceed"
  
  gplId <- which(grepl("!Sample_platform_id", wholeFile))
  gplId <- gsub(".*(GPL\\d+).*", "\\1", wholeFile[gplId])
  gpl <-  paste(c(gplDir, "/", gplId, ".3col.gz"), collapse = "")
  
  # it is necessary that tag has gpl
  tag <- inFile %>% str_extract("GSE\\d+")
  tag <- paste(c(tag, "_", gplId), collapse = "")
  
  if (!file.exists(gpl))
  {
    processingStatus <-
      gsub("gplId", gplId, "Failed. gplId is not supported.")
    writeProcessingStatus(outDir, inFile, tag, processingStatus, NA, NA, NA)
  } else {
    isSuperSeries <- grepl("SuperSeries of", wholeFile) %>% any()
    
    if (isSuperSeries)
    {
      processingStatus <-
        gsub("tag", tag, "Failed. tag is a superseries.")
      writeProcessingStatus(outDir, inFile, tag, processingStatus, NA, NA, NA)
    } else {

      
      tableStart <-
        which(grepl("!series_matrix_table_begin", wholeFile))
      tableEnd <-
        which(grepl("!series_matrix_table_end", wholeFile))
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
      colnames(expressionTable) <-
        mapply(
          colnames(expressionTable),
          sampleTitles,
          FUN = function(x, y)
            sprintf("%s_%s", x, y)
        )
      
  if (ncol(expressionTable) < 8 |
      ncol(expressionTable) > 200)
  {
    print("Failed due to number of samples been out of range")
  } else {
    exp <- linearizeDataset(expressionTable)
    explog <- logDataset(expressionTable)

    logav <- mean(explog)
    logmax <- max(explog)
    linmax <- max(exp)
        
  if (!(!any(is.na(c(
    logav, logmax, linmax
  ))) &&
  logav >= 4 &&
  logav <= 10 &&
  linmax >= 2000 &&
  logmax < 20))
  {
    processingStatus <- "Failed sanity check"
    writeProcessingStatus(outDir, inFile, tag, processingStatus, NA, NA, NA)
  } else {
          
  gplTable <- readGPLTable(gpl)

  exp <- cbind(gplTable[rownames(exp), ], exp)

  exp <-
    exp[!(row.names(exp) %in% c("NONE", "NA", "none", "NULL")), ]

  exp <-
    collapseRows(exp[4:ncol(exp)], exp$Entrez_ID, rownames(exp))

  exp <- exp$datETcollapsed

  exp <-
    exp[!(row.names(exp) %in% c("NONE", "NA", "none", "NULL")),]

  exp <- apply(exp, 2, function(x) (x/sum(x))*10^6)
          

  if (expMatF)
  {
    write.table(
      exp,
      file = paste(outDir, "/matrices/", tag, ".tsv", sep = ""),
      quote = F,
      col.names = T,
      row.names = T,
      sep = "\t"
    )
  }