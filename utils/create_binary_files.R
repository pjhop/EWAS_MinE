## Create binary files (OSCA) from a beta (or mvalues) matrix

## Libraries ---------------------------
library(tidyverse)
library(magrittr)
library(optparse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## Data ---------------------------
option_list <- list(
  make_option(c("-b", "--beta"), action="store", type = "character", default=NULL,
              help="Path to beta matrix (RDS file)"),
  make_option(c("-m", "--merge"), action="store", type = "character", default=NULL,
              help="(Optional) Txt-file containing paths to multiple beta/mvalue files that should be merged"),
  make_option(c("-t", "--type"), action="store", type = "character", default="beta",
              help="Type of methylation values (either beta or mvalue)"),
  make_option(c("--log"), action="store", type = "character", default="logfile.txt",
              help="Logfile"),
  make_option(c("-a", "--array"), action="store", type = "character", default="450k",
              help="Beadchip (450k or EPIC)"),
  make_option(c("-n", "--opi"), action="store", type = "character", default=NULL,
               help="Path to opi-file"),
  make_option(c("-r", "--removesamples"), action="store", type = "character", default=NULL,
              help="Filepath to list of samples that should be removed (optional)"),
  make_option(c("-p", "--removeprobes"), action="store", type = "character", default=NULL,
              help="Filepath to list of probes that should be removed (optional)"),
  make_option(c("-s", "--sex"), action="store", type = "character", default="keep",
              help="What to do with sex-chromosomal probes? Three options: keep, remove or sep (sep will store these probes in a separate file"),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Path to write the binary files to")
             )

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

## Load betas/mvalues

# If txt-file with multiple beta-matrices is specified, these are loaded and merged
if(!is.null(opt$merge)) {
	filepaths <- readr::read_lines(opt$merge)
	cat(sprintf("\nLoading %s files..", length(filepaths)))
	fun <- function(path) {
		cat(sprintf("\nLoading %s", path))
		file <- readRDS(path)
	}
	files <- purrr::map(filepaths, fun)

	# Find intersecting probes
	get_probes <- function(file) {
		probes <- rownames(file)
		probes
	}

	all_probes <- purrr::map(files, get_probes)
	probes_intersect <- Reduce(intersect, all_probes)
	cat(sprintf("\n%s probes overlap between the methylation files.", length(probes_intersect)))

	cat("\nMerging data..")
	files <- purrr::map(files, .f = function(x) x[probes_intersect, ])
	beta <- do.call(cbind,files)
	cat("\nDone!")
	rm(files); gc(verbose = FALSE)


} else if(!is.null(opt$beta)) {
	cat("\nLoading beta/mvalue matrix..")
	beta <- readRDS(opt$beta)
} else{
	stop("Path to a beta/mvalue matrix or path to list of files should be specified!")
}


## Remove samples and/or probes
if(!is.null(opt$removesamples)) {
	removesamples <- readr::read_lines(opt$removesamples)
	cat(sprintf("\nRemoving %s samples from the data", n_distinct(removesamples)))
	beta <- beta[,!(colnames(beta) %in% removesamples)]
} else {
	cat("\n--removesamples argument not specified, not removing any samples.")
}

if(!is.null(opt$removeprobes)) {
	removeprobes <- readr::read_lines(opt$removeprobes)
	cat(sprintf("\nRemoving %s probes from the data", n_distinct(removeprobes)))
	beta <- beta[!(rownames(beta) %in% removeprobes),]
} else {
	cat("\n--removeprobes argument not specified, not removing any probes")
}


# XY probes
if(opt$array == "450k") {
  anno <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  xy <- anno %>% dplyr::filter(chr %in% c("chrX", "chrY")) %$% Name
} else {
  anno <- data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  xy <- anno %>% dplyr::filter(chr %in% c("chrX", "chrY")) %$% Name
}


## Write betas to txt-file  ---------------------------
samples <- colnames(beta)

if(opt$sex == "remove") {
	cat(sprintf("\nRemoving %s sex-chromosomal probes.\n\n", sum((rownames(beta) %in% xy))))
	beta <- beta[!(rownames(beta) %in% xy),]
} else if(opt$sex == "sep") {
	cat(sprintf("\nCreating a separate file for %s sex-chromosomal probes\n\n", length(xy)))
	beta_xy <- beta[rownames(beta) %in% xy,]
	beta <- beta[!(rownames(beta) %in% xy),]
	beta_xy <- t(beta_xy)
	beta_xy <- dplyr::as_tibble(beta_xy)
	beta_xy <- dplyr::bind_cols(data.frame(samples), data.frame(samples), beta_xy)
	colnames(beta_xy)[1] <- 'IID'
	colnames(beta_xy)[2] <- 'FID'
  data.table::fwrite(beta_xy, file = sprintf("%s_xy.txt", opt$out), sep = "\t", na = "NA", quote = FALSE)
}

# Transpose
beta <- t(beta)
beta <- dplyr::as_tibble(beta)

beta <- dplyr::bind_cols(data.frame(samples), data.frame(samples), beta)
colnames(beta)[1] <- 'IID'
colnames(beta)[2] <- 'FID'

# Write betas to a txt file (this txt file is removed after creating the binary files)
cat("\nWriting betas to txt-file..")
#readr::write_tsv(beta, path = sprintf("%s.txt", opt$out), col_names = TRUE)
data.table::fwrite(beta, file = sprintf("%s.txt", opt$out), sep = "\t", na = "NA", quote = FALSE)
cat("\nDone!")

name <- tail(stringr::str_split(opt$out, pattern = "/")[[1]], n = 1)

# Make binary files
if(opt$type == "beta") {
  call <- sprintf("~/osca_Linux\
                  --efile %1$s\
                  --methylation-beta\
                  --make-bod\
                  --out %2$s > %3$s", sprintf("%s.txt", opt$out), opt$out, opt$log)
  call <- gsub("\n", " ", call)
} else if(opt$type == "mvalues") {
  call <- sprintf("~/osca_Linux\
                  --efile %1$s\
                  --methylation-m\
                  --make-bod\
                  --out %2$s > %3$s", sprintf("%s.txt", opt$out), opt$out, opt$log)
  call <- gsub("\n", " ", call)
} else {
  call <- sprintf("~/osca_Linux\
                  --efile %1$s\
                  --make-bod\
                  --out %2$s > %3$s", sprintf("%s.txt", opt$out), opt$out, opt$log)
  call <- gsub("\n", " ", call)
}

# Run
cat("\nCreating binary files (OSCA)..")
system(call)
cat("\nDone!")

# Update .opi file
cat("\nUpdating OPI file..")
if(opt$array == "450k") {
	system(sprintf("~/osca_Linux --befile %s --update-opi %s",opt$out, opt$opi), ignore.stdout = TRUE, ignore.stderr = TRUE)
} else if(opt$array == "EPIC") {
	system(sprintf("~/osca_Linux --befile %s --update-opi %s",opt$out, opt$opi ), ignore.stdout = TRUE, ignore.stderr = TRUE)
}
cat("\nDone!..")

## Save sex-chromomal probes separately (if opt$sex = 'sep')
if(opt$sex == "sep") {
	if(opt$type == "beta") {
	call <- sprintf("~/osca_Linux\
--efile %1$s \
--methylation-beta \
--make-bod \
--out %2$s >> %3$s", sprintf("%s_xy.txt", opt$out), sprintf("%s_xy",opt$out), opt$log)
	call <- gsub("\n", " ", call)
	cat("\nCreating binary files for xy probes..")
	system(call)
	if(opt$array == "450k") {
		system(sprintf("~/osca_Linux --befile %s_xy --update-opi %s",opt$out, opt$opi ), ignore.stdout = TRUE, ignore.stderr = TRUE)
	} else if(opt$array == "EPIC") {
		system(sprintf("~/osca_Linux --befile %s_xy --update-opi %s",opt$out, opt$opi ), ignore.stdout = TRUE, ignore.stderr = TRUE)
	}
	cat("\nDone!")

} else if (opt$type == "mvalues") {
	call <- sprintf("~/osca_Linux\
--efile %1$s\
--methylation-m\
--make-bod\
--out %2$s >> %3$s", sprintf("%s_xy.txt", opt$out), sprintf("%s_xy",opt$out), opt$log)
	call <- gsub("\n", " ", call)
	cat("\nCreating binary files for xy probes..")
	system(call)
	if(opt$array == "450k") {
		system(sprintf("~/osca_Linux --befile %s_xy --update-opi %s",opt$out, opt$opi), ignore.stdout = TRUE, ignore.stderr = TRUE)
	} else if(opt$array == "EPIC") {
		system(sprintf("~/osca_Linux --befile %s_xy --update-opi %s",opt$out, opt$opi ), ignore.stdout = TRUE, ignore.stderr = TRUE)
	}
	cat("\nDone!")
} else {
  call <- sprintf("~/osca_Linux \
--efile %1$s \
--make-bod \
--out %2$s >> %3$s", sprintf("%s_xy.txt", opt$out), sprintf("%s_xy",opt$out), opt$log)
  call <- gsub("\n", " ", call)
  cat("\nCreating binary files for xy probes..")
  system(call)
  if(opt$array == "450k") {
    system(sprintf("~/osca_Linux --befile %s_xy --update-opi %s",opt$out, opt$opi), ignore.stdout = TRUE, ignore.stderr = TRUE)
  } else if(opt$array == "EPIC") {
    system(sprintf("~/osca_Linux --befile %s_xy --update-opi %s",opt$out, opt$opi ), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
  cat("\nDone!")
}
}

cat(sprintf("\n\nAnalysis finished at: %s\n\n", Sys.time()))
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
