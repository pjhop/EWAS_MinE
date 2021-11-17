## Perform functional normalization using the meffil package, optionally including random effects

## Libraries ---------------------------

library(meffil)
library(tidyverse)
library(magrittr)
library(optparse)
library(stringr)

## Data ---------------------------
option_list <- list(
  make_option(c("-q", "--qcobjects"), action="store", type = "character", default=NULL,
              help="Path to qc objects (rds file)"),
  make_option(c("-r", "--removesamples"), action="store", type = "character", default=NULL,
              help="Txt file containing samples that should be removed"),
  make_option(c("-p", "--removeprobes"), action="store", type = "character", default=NULL,
              help="Txt file containing probes that should be removed"),
  make_option(c("-v", "--variables"), action="store", type = "character", default=NULL,
              help="Additional variables that should be included as random effects. If more than one, use commas to separate the variable names"),
  make_option(c("-n", "--npcs"), action="store", type = "numeric", default=NULL,
              help="Number of PCs to use for functional normalization"),
  make_option(c("-c", "--cores"), action="store", type = "numeric", default=NULL,
              help="Number of cores that should be used"),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Path to directory where to store the plot and data")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
options(mc.cores = opt$cores)

## Print options
cat("Following options are specified:\n\n")
for(i in names(opt)) {
  if(i != "help") {
    cat(sprintf("--%s %s\n", i, opt[[i]]))
  }
}
cat("\n")

qc.objects <- readRDS(opt$qcobjects)

## Normalization ---------------------------

# Remove failed samples/probes
if(!is.null(opt$removesamples)) {
  remove_samples <- read_lines(opt$removesamples)
  cat(sprintf("Removing %s samples.\n", length(unique(remove_samples))))
  qc.objects.clean <- meffil.remove.samples(qc.objects, remove_samples)
} else {
  remove_samples <- c()
  qc.objects.clean <- qc.objects
}

if(!is.null(opt$removeprobes)) {
  remove_probes <- read_lines(opt$removeprobes)
} else {
  remove_probes <- c()
}

rm(qc.objects); gc()

## Normalization  (no fixed effects option at the moment)
if(!is.null(opt$variables)) {
  variables <- str_trim(unlist(str_split(opt$variables, pattern = ",")))
  cat(sprintf("Normalizing using %s as random effect..\n", paste(variables, collapse = " and ")))
  norm.objects <- meffil.normalize.quantiles(qc.objects.clean, random.effects = variables, fixed.effects = NULL, number.pcs = opt$npcs)
  cat("Done!")
} else {
  cat("Normalizing without including fixed or random effects..\n")
  norm.objects <- meffil.normalize.quantiles(qc.objects.clean, random.effects = NULL, fixed.effects = NULL, number.pcs = opt$npcs)
  cat("Done!")
}

## Save
saveRDS(norm.objects , file = sprintf("%s.objects.rds", opt$out))

## Normalized beta matrix
cat(sprintf("\nCreating a normalized beta matrix, removing %s probes..\n\n", length(remove_probes)))
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = remove_probes)
cat("Done!")

## Save
cat(sprintf("\n\nSaving the normalized beta-matrix to: %s\n", opt$out))
saveRDS(norm.beta, file = sprintf("%s.beta.rds", opt$out))

## SessionInfo
cat(sprintf("Analysis finished at: %s\n\n", Sys.time()))
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()

