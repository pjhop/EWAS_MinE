## Create ORM and perform PCA using OSCA

## Libraries ---------------------------
library(tidyverse)
library(optparse)

## Data ---------------------------
option_list <- list(
  make_option(c("-b", "--befile"), action="store", type = "character", default=NULL,
              help="Path to binary files"),
  make_option(c("-o", "--orm"), action="store", type = "character", default=NULL,
              help="Path to where orm should be written"),
  make_option(c("--log"), action="store", type = "character", default="logfile.txt",
              help="Logfile"),
  make_option(c("-r", "--remove"), action="store", type = "character", default=NULL,
              help="txt-file containing samples to remove"),
  make_option(c("-k", "--keep"), action="store", type = "character", default=NULL,
              help="txt-file containing samples to keep"),
  make_option(c("-a", "--removeprobes"), action="store", type = "character", default=NULL,
              help="txt-file containing probes to remove"),
  make_option(c("--keepprobes"), action="store", type = "character", default=NULL,
              help="txt-file containing probes to keep"),
  make_option(c("-n", "--npcs"), action="store", type = "numeric", default=20,
              help="Nr of PCs to calculate"),
  make_option(c("-p", "--pca"), action="store", type = "character", default=NULL,
              help="Path to were PCs should be written)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat(sprintf("Analysis started at: %s\n\n", Sys.time()))

# Print options
cat("Following options are specified:\n\n")
for(i in names(opt)) {
    if(i != "help") {
        cat(sprintf("--%s %s\n", i, opt[[i]]))
     }
}
cat("\n")

# Remove/keep samples
if(!is.null(opt[["remove"]])) {
	remove <- readr::read_lines(opt$remove)
  cat(sprintf("%s indivuals will be kept from %s\n", dplyr::n_distinct(remove), opt$remove))
  readr::write_tsv(tibble(IID = remove, FID = remove), path = sprintf("%s_remove", opt$orm),
                   col_names = FALSE)
      }

if(!is.null(opt[["keep"]])) {
	keep <- readr::read_lines(opt[["keep"]])
  cat(sprintf("%s individuals will be kept from %s\n", dplyr::n_distinct(keep), opt[["keep"]]))
  readr::write_tsv(tibble(IID = keep, FID = keep), path = sprintf("%s_keep", opt$orm),
                   col_names = FALSE)
}

## Load opi-file
opi <- readr::read_tsv(paste0(opt$befile, ".opi"), col_names = FALSE)

if(!is.null(opt[["removeprobes"]])) {
  remove_probes <- readr::read_lines(opt[["removeprobes"]])
  cat(sprintf("%s probes will be removed from %s\n", 
              sum(unique(remove_probes) %in% opi$X2), opt[["removeprobes"]]))
} else {
  remove_probes <- c()
}

if(!is.null(opt[["keepprobes"]])) {
  keep_probes <- readr::read_lines(opt[["keepprobes"]])
  cat(sprintf("%s probes will be kept from %s\n", 
              sum(unique(keep_probes) %in% opi$X2), opt[["keepprobes"]]))
} else {
  keep_probes <- opi$X2
}

## Total keep list 
keep_probes <- keep_probes[!keep_probes %in% remove_probes]
cat(sprintf("In total, %s out of %s probes will be kept\n", 
            sum(unique(keep_probes) %in% opi$X2), length(opi$X2)))
readr::write_lines(keep_probes, path = sprintf("%s_keep_probes", opt$orm))

## Calculate ORM
cat("Calculating ORM..\n")
if(!is.null(opt[["remove"]]) && !is.null(opt[["keep"]])) {
  system(sprintf("~/osca_Linux --befile %s --make-orm --keep %s --remove %s --extract-probe %s --out %s > %s",
                 opt$befile, sprintf("%s_keep", opt$orm), sprintf("%s_remove", opt$orm), sprintf("%s_keep_probes", opt$orm),  opt$orm, opt$log))
} else if(!is.null(opt[["remove"]]) && is.null(opt[["keep"]])) {
  system(sprintf("~/osca_Linux --befile %s --make-orm --remove %s --extract-probe %s --out %s > %s",
                 opt$befile, sprintf("%s_remove", opt$orm), sprintf("%s_keep_probes", opt$orm), opt$orm, opt$log))
} else if(is.null(opt[["remove"]]) && !is.null(opt[["keep"]])) {
  system(sprintf("~/osca_Linux --befile %s --make-orm --keep %s --extract-probe %s --out %s > %s",
                 opt$befile, sprintf("%s_keep", opt$orm), sprintf("%s_keep_probes", opt$orm), opt$orm, opt$log))
} else {
  system(sprintf("~/osca_Linux --befile %s --make-orm --extract-probe %s --out %s > %s",
                 opt$befile, sprintf("%s_keep_probes", opt$orm), opt$orm, opt$log))
}

## Perform PCA
cat("Performing PCA..\n")
system(sprintf("~/osca_Linux --pca %s --orm %s --out %s >> %s",
               opt$npcs, opt$orm, opt$pca, opt$log))

## SessionInfo
cat(sprintf("Analysis finished at: %s\n\n", Sys.time()))
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
