## Libraries ---------------------------

library(minfi)
library(wateRmelon)
library(readr)
library(magrittr)
library(optparse)
library(dplyr)

## Arguments ---------------------------

option_list <- list(
  make_option(c("-r", "--mset"), action="store", type = "character", default=NULL,
              help="Path to Mset (rds file)"),
  make_option(c("-s", "--removesamples"), action="store", type = "character", default=NULL,
              help="List of samples to remove"),
  make_option(c("-k", "--keepsamples"), action="store", type = "character", default=NULL,
              help="List of samples to keep"),
  make_option(c("-p", "--removeprobes"), action="store", type = "character", default=NULL,
              help="List of probes to remove"),
  make_option(c("-t", "--type"), action="store", type = "character", default="dasen",
              help="Which normalization method to use, dasen/nanes"),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Path to write the normalized (dasen) intensities toe")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Print options
cat("Following options are specified:\n\n")
for(i in names(opt)) {
  if(i != "help") {
    cat(sprintf("--%s %s\n", i, opt[[i]]))
  }
}
cat("\n")

if(!is.null(opt$removesamples)) {
  remove_samples <- read_lines(opt$removesamples)
} else {
  remove_samples <- c()
}

if(!is.null(opt$removeprobes)) {
  remove_probes <- read_lines(opt$removeprobes)
} else {
  remove_probes <- c()
}

## Load data ---------------------------
mset <- readRDS(opt$mset)

if(!is.null(opt$keepsamples)) {
  keep_samples <- read_lines(opt$keepsamples)
} else {
  keep_samples <- colnames(mset)
}

## Remove failed samples ---------------------------
keep <- setdiff(keep_samples, remove_samples)
ncol1 <- ncol(mset)
mset <- mset[,(colnames(mset) %in% keep)]
ncol2 <- ncol(mset)
cat(sprintf("\nRemoved %s samples from the mset!", ncol1-ncol2))
cat(sprintf("\nKept %s out of %s samples!", ncol2, ncol1))

## Get normalized intensities  ---------------------------
cat(sprintf("\nRunning %s function..", opt$type))
if(opt$type == "dasen") {
  mset <- dasen(mset)
} else if(opt$type == "nanes") {
  mset <- nanes(mset)
} else {
  stop("type: %s not supported", opt$type)
}

M <- mset@assayData$methylated
U <- mset@assayData$unmethylated
rm(mset)
Total <- M + U

cat("\nDone!")

## Save intensities  ---------------------------
saveRDS(M, file = sprintf("%s_methylated.rds", opt$out))
saveRDS(U, file = sprintf("%s_unmethylated.rds", opt$out))
saveRDS(Total, file = sprintf("%s_total.rds", opt$out))

## SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
