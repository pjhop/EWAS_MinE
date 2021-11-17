options(stringsAsFactors = FALSE)
## Get qc.objects from the raw iDat files (meffil R package)

## Libraries ---------------------------

library(meffil)
library(minfi)
library(optparse)

## Arguments ---------------------------

option_list <- list(
  make_option(c("-c", "--cores"), action="store", type = "numeric", default=1,
              help="Number of cores to use"),
  make_option(c("-s", "--samplesheet"), action="store", type = "character", default=NULL,
              help="Path to samplesheet (stored as rds-file)."),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Path to write the qc objects to (rds extension)")
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

## Data ---------------------------
samplesheet <- readRDS(opt$samplesheet)

## Get qc.objects ---------------------------
cat("\nReading the iDat-files...\n")
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069 complete", verbose = TRUE)
cat("\n\n Done!\n")
cat(sprintf("Created a qc.object file for %s samples!", length(qc.objects)))

## Save ---------------------------
cat(sprintf("\n\nSaving the qc-objects to: %s\n", opt$out))
saveRDS(qc.objects, file = opt$out)

## SessionInfo
cat(sprintf("Analysis finished at: %s\n\n", Sys.time()))
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
