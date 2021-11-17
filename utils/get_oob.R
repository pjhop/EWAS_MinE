## Get normalized OOB betas using the DNAmCrosshyb R package

## Libraries ---------------------------
library(optparse)
library(minfi)
library(DNAmCrosshyb)

## Arguments ---------------------------
option_list <- list(
  make_option(c("-s", "--samplesheet"), action="store", type = "character", default=NULL,
              help="Path to samplesheet (stored as rds-file). The samplesheet should containg a 'Sample_Name' column"),
  make_option(c("-r", "--rgset"), action="store", type = "character", default=NULL,
              help="Path to RGset"),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Path to write OOB betas to")
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

# Load rgset
rgset <- readRDS(opt$rgset)
samplesheet <- readRDS(opt$samplesheet)

# Keep
rgset <- rgset[,samplesheet$Sample_Name]

# Get OOB values
beta_oob <- DNAmCrosshyb::get_OOB(rgset, keep = "beta", normalized = TRUE)
beta_oob <- beta_oob$beta

# Save
saveRDS(beta_oob, opt$out)