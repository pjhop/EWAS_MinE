### Convert beta-values to m-values

## Libraries
library(optparse)

## Data ---------------------------
option_list <- list(
  make_option(c("-b", "--beta"), action="store", type = "character", default=NULL,
              help="Path to beta matrix (RDS file)"),
  make_option(c("-m", "--mvalues"), action="store", type = "character", default=NULL,
              help="Path where to store mvalue matrix (RDS file)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cat(sprintf("Analysis started at: %s\n\n", Sys.time()))

# Load betas
beta <- readRDS(opt$beta)

## Transform ---------------------------

# Check if there are any 0 or 1 values
onevalues <- beta == 1
zerovalues <- beta == 0

# Set these values to the maximum (or minimum) value non-one/non-zero value
if(sum(onevalues, na.rm = TRUE) > 0 || sum(zerovalues, na.rm = TRUE) > 0 ) {

  if(sum(onevalues, na.rm = TRUE) > 0) {
    beta[onevalues] <- NA
    beta[onevalues] <- max(beta, na.rm = TRUE)
  }
  if(sum(zerovalues, na.rm = TRUE) > 0)  {
    beta[zerovalues] <- NA
    beta[zerovalues] <- min(beta, na.rm = TRUE)
  }
}

# Transform
beta <- log2(beta/(1-beta))

## Save ---------------------------

saveRDS(beta, file = opt$mvalues)

cat(sprintf("Analysis finished at: %s\n\n", Sys.time()))
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
