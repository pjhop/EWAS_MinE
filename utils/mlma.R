## R-wrapper around OSCA (linux executable)

## Libraries ---------------------------

library(tidyverse)
library(magrittr)
library(optparse)
library(Matrix)

## Parse parameters ---------------------------

option_list <- list(
  make_option(c("-b", "--befile"), action="store", type = "character", default=NULL,
              help="Path to binary files"),
  make_option(c("-s", "--samplesheet"), action="store", type = "character", default=NULL,
              help="Path to samplesheet (RDS file)"),
  make_option(c("--log"), action="store", type = "character", default="logfile.txt",
              help="Logfile."),
  make_option(c("-r", "--removesamples"), action="store", type = "character", default=NULL,
              help="Txt file containing samples that should be removed"),
  make_option(c("--keepsamples"), action="store", type = "character", default=NULL,
              help="Txt file containing samples that should be mept"),
  make_option(c("-a", "--removeprobes"), action="store", type = "character", default=NULL,
              help="Txt file containing probes that should be removed"),
  make_option(c("-x", "--removexy"),
  			  action="store", type = "character", default="n",
              help="Remove xy probes? (y/n)"),
  make_option(c("--filtervariable"), action="store", type = "character", default=NULL,
              help="Variable to filter on"),
  make_option(c("--filtervariablekeep"), action="store", type = "character", default=NULL,
              help="Value of filter variable to keep"),
  make_option(c("--filtervariabledrop"), action="store", type = "character", default=NULL,
              help="Value of filter variable to drop"),
  make_option(c("-p", "--pheno"), action="store", type = "character", default=NULL,
              help="Name of column that contains phenotype, should be present in samplesheet"),
  make_option(c("-l", "--levels"), action="store", type = "character", default=NULL,
              help="If phenotype is binary, specifiy the order of levels. Seperate levels by a comma"),
  make_option(c("-c", "--covariates"), action="store", type = "character", default=NULL,
              help="Categorical covariates to include, should be present in samplesheet"),
  make_option(c("--missingness"), action="store", type = "character", default=NULL,
              help="Sparse matrix containing missingness"),
  make_option(c("--missthreshold"), action="store", type = "numeric", default=1,
              help="Probes with missingness > threshold will be removed"),
  make_option(c("-q", "--qcovariates"), action="store", type = "character", default=NULL,
              help="Quantative c ovariates to include, should be present in samplesheet"),
  make_option(c("-n", "--nmatrix"), action="store", type = "numeric", default=NULL,
              help="Number of covariates to include from matrix specified in --matrix"),
  make_option(c("-d", "--method"), action="store", type = "character", default="moa",
              help="Which method? Should be moa/moment/loco/linear/logistic"),
  make_option(c("-v", "--cutofform"), action="store", type = "numeric", default=1,
              help="Cutoff for orm"),
  make_option(c("--orm"), action="store", type = "character", default=NULL,
              help="Precomputed orm. If not specified, will compute the orm in this script.
              Note that using a variance cutoff to filter probe before calculating the ORM, seems to be only possible when using
              a precomputed ORM"),
  make_option(c("--momentpercent"), action="store", type = "numeric", default=NULL,
              help="Fraction of probes to fit into one ORM"),
  make_option(c("--momentwind"), action="store", type = "numeric", default=NULL,
              help="MOMENT flanking region"),
  make_option(c("--osca"), action="store", type = "character", default="~/osca_Linux",
              help="Path to osca binary"),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Filename of output")
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

## Load data -------------------------------------------------------------------
samplesheet <- as.data.frame(readRDS(opt$samplesheet))
covariates <- str_trim(unlist(str_split(opt$covariates, pattern = ",")))
qcovariates <- str_trim(unlist(str_split(opt$qcovariates, pattern = ",")))
pheno <- opt$pheno
levels <- str_trim(unlist(str_split(opt$levels, pattern = ",")))

if(!is.null(opt$filtervariablekeep)) {
  filtervariablekeep <- str_trim(unlist(str_split(opt$filtervariablekeep, pattern = ",")))
}

if(!is.null(opt$filtervariabledrop)) {
  filtervariabledrop <- str_trim(unlist(str_split(opt$filtervariabledrop, pattern = ",")))
}

if(!is.null(opt$removesamples)) {
	remove_samples <- read_lines(opt$removesamples)
} else {
	remove_samples <- c()
}

if(!is.null(opt$keepsamples)) {
  keep_samples <- readr::read_lines(opt$keepsamples)
}

if(!is.null(opt$removeprobes)) {
	remove_probes <- read_lines(opt$removeprobes)
} else {
	remove_probes <- c()
}

if(!is.null(opt$missingness)) {
  missingness <- readRDS(opt$missingness)
}

if(opt$removexy == "y") {
	library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
	anno <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
	xy <- anno %>% filter(chr %in% c("chrX", "chrY")) %$% Name
	anno <- data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
	xy_EPIC <- anno %>% filter(chr %in% c("chrX", "chrY")) %$% Name
	xy <- unique(c(xy, xy_EPIC))
	rm(anno)
} else {
	xy <- c()
}

# Check if phenotype and covariates are present in samplesheet
if(!all(covariates %in% colnames(samplesheet))) {
	stop("Not all categorical covariates are present in samplesheet")
}

if(!all(qcovariates %in% colnames(samplesheet))) {
	stop("Not all quantitative covariates are present in samplesheet")
}
if(!(pheno %in% colnames(samplesheet))) {
	stop("Phenotype is not present in samplesheet")
}

## Filter data --------------------------------------------------------

## Samples
cat("Samples----------\n\n")
cat(sprintf("Removing %s samples specified by user.\n", sum(remove_samples %in% samplesheet$Sample_Name)))
samplesheet <- samplesheet[!(samplesheet$Sample_Name %in% remove_samples),]

if(!is.null(opt$keepsamples)) {
  cat(sprintf("Keeping %s samples as specified by user.\n", sum(keep_samples %in% samplesheet$Sample_Name)))
  samplesheet <- samplesheet[(samplesheet$Sample_Name %in% keep_samples),]
}

## Filter on a specific variable
if(!is.null(opt$filtervariable)) {
  if(!is.null(opt$filtervariablekeep)) {
    cat(sprintf("Keeping samples where %s == %s\n", opt$filtervariable, filtervariablekeep))
    samplesheet <- samplesheet %>% dplyr::filter(!!as.name(opt$filtervariable) %in% filtervariablekeep)
  }
  if(!is.null(opt$filtervariabledrop)) {
    cat(sprintf("Dropping samples where %s != %s\n", opt$filtervariable, filtervariabledrop))
    samplesheet <- samplesheet %>% dplyr::filter(!(!!as.name(opt$filtervariable) %in% filtervariabledrop))
  }
}

## Load oii file
oii <- readr::read_tsv(sprintf("%s.oii", opt$befile), col_names = FALSE)
cat(sprintf("%s out of %s samples in the samplesheet are present in the binary files.\n",
            sum(samplesheet$Sample_Name %in% oii$X1), nrow(samplesheet)))
samplesheet <- samplesheet[samplesheet$Sample_Name %in% oii$X1,, drop = FALSE]

## Load opi file
opi <- readr::read_tsv(sprintf("%s.opi", opt$befile), col_names = FALSE)

## Phenotype levels
if(!is.null(opt$levels)) {
	cat(sprintf("%s will be coded as 1, and %s will be coded as 0.\n", levels[2], levels[1]))
	cat(sprintf("%s column contains %s sample(s) that have missing data for this variable, these samples will be removed.\n",
		pheno, sum(is.na(samplesheet[,pheno]))))
	cat(sprintf("%s column contains %s samples(s) with a level that is different than the specified levels, these samples will be removed.\n",
		pheno, sum(!(samplesheet[,pheno] %in% levels))))
	samplesheet <- samplesheet[!is.na(samplesheet[,pheno]) & (samplesheet[,pheno] %in% levels),]
	samplesheet[,pheno] <- factor(samplesheet[,pheno], levels = levels)
} else {
	cat(sprintf("%s column contains %s sample(s) that have missing data for this variable, these samples will be removed.\n",
		pheno, sum(is.na(samplesheet[,pheno]))))
	samplesheet <- samplesheet[!is.na(samplesheet[,pheno]),]

}

# Covariates - remove missings
cat(sprintf("%s samples have missing data for one of the categorical covariates, these samples will be removed.\n",
	sum(!complete.cases(samplesheet[,covariates]))))
samplesheet <- samplesheet[complete.cases(samplesheet[,covariates]),]

cat(sprintf("%s samples have missing data for one of the quantative covariates, these samples will be removed.\n\n",
	sum(!complete.cases(samplesheet[,qcovariates]))))
samplesheet <- samplesheet[complete.cases(samplesheet[,qcovariates]),]

## Remove probes based on missingness / user-specified
cat("Probes----------\n\n")
cat(sprintf("%s probes specified by user will be removed.\n", sum(opi$X2 %in% remove_probes)))
cat(sprintf("%s probes on the sex chromosomes will be removed.\n", sum(opi$X2 %in% xy)))

if(!is.null(opt$missingness)) {
  missingness <- missingness[ ,samplesheet$Sample_Name]
  prop_missingness <- Matrix::rowMeans(missingness)
  prop_missingness <- names(prop_missingness[prop_missingness > opt$missthreshold])
  cat(sprintf("%s probes with missingness > %s will be removed\n", sum(opi$X2 %in% prop_missingness), opt$missthreshold))
} else {
  prop_missingness <- c()
}

## Write temporary files ---------------------------

# Excluded probes
remove_probes_total <- unique(c(remove_probes, xy, prop_missingness))

# Quantative covariates
if(!is.null(opt$qcovariates)) {
	qcovar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, samplesheet[,qcovariates])
}

# Categorical covariates
if(!is.null(opt$covariates)) {
	covar <- apply(samplesheet[,covariates,drop=FALSE], 2, FUN = function(x) as.numeric(as.factor(x)))
	nr.levels <- apply(covar, 2, FUN =  dplyr::n_distinct)
	nr.levels <- nr.levels[nr.levels == 1]
	if(length(nr.levels) > 0) {
	  cat("The following covariates have only one level and will be dropped:\n")
	  for(i in length(nr.levels)) {
	    cat(names(nr.levels)[i], "\n")
	  }
	  covar <- covar[,!colnames(covar) %in% names(nr.levels), drop = FALSE]
	  if(ncol(covar) == 0) {
	    cat("\nThere are no qualitative covariates left after removing covariates with only one level.\n")
	    opt$covariates <- NULL
	  } else {
	    covar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, covar)
	  }
	} else {
	  covar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, covar)
	}
}

# Phenotype
phenotype <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, ifelse(samplesheet[,pheno] == levels[2], 1, 0))

## Create a directory where to put these files
dir.create(path = sprintf("%s_input", opt$out))
write_lines(remove_probes_total, path = sprintf("%s_input/remove_probes.txt", opt$out))
if(!is.null(opt$qcovariates)) {
	write_tsv(qcovar, path = sprintf("%s_input/qcovar.txt", opt$out), col_names=FALSE)
}
if(!is.null(opt$covariates)) {
	write_tsv(covar, path = sprintf("%s_input/covar.txt", opt$out), col_names=FALSE)
}

write_tsv(phenotype, path = sprintf("%s_input/phenotype.txt", opt$out), col_names=FALSE)


## System call ---------------------------

cat("\nOSCA----------\n")
cat("Fitting..\n")

type <- case_when(!is.null(opt$qcovariates) & !is.null(opt$covariates) ~ "covar_qcovar",
				  !is.null(opt$qcovariates) ~ "qcovar",
				  !is.null(opt$covariates) ~ "covar",
				  is.null(opt$covariates) & is.null(opt$qovariates) ~ "None"
				)

name <- tail(str_split(opt$out, pattern = "/")[[1]], n = 1)

if(opt$method == "moa") {
	if(!is.null(opt$orm)) {
		call <- switch(type,
	qcovar = sprintf("%6$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm %3$s \
--orm-cutoff %4$s \
--out %2$s > %5$s", opt$befile, opt$out,opt$orm, opt$cutofform, opt$log, opt$osca),
covar = sprintf("%6$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm %3$s \
--orm-cutoff %4$s \
--out %2$s > %5$s", opt$befile, opt$out,opt$orm, opt$cutofform, opt$log, opt$osca),
	covar_qcovar = sprintf("%6$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm %3$s \
--orm-cutoff %4$s \
--out %2$s > %5$s", opt$befile, opt$out,opt$orm, opt$cutofform, opt$log, opt$osca),
	None = sprintf("%6$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm %3$s \
--orm-cutoff %4$s \
--out %2$s > %5$s", opt$befile, opt$out,opt$orm, opt$cutofform, opt$log, opt$osca)
	)
	call <- gsub("\n", " ", call)
	system(call)
	} else {
		call <- switch(type,
	qcovar =  sprintf("%5$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform,opt$log, opt$osca),
	covar = sprintf("%5$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform,opt$log, opt$osca),
	covar_qcovar = sprintf("%5$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--qcovar %2$s_input/qcovar.txt \
--covar %2$s_input/covar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform,opt$log, opt$osca),
	None = sprintf("%5$s --moa \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform,opt$log, opt$osca))
		call <- gsub("\n", " ", call)
		system(call) }
	} else if(opt$method %in%  c("loco", "moment", "moment2")){
	  argname <- dplyr::case_when(
	    opt$method == "loco" ~ "--mlma-loco",
	    opt$method == "moment" ~ "--moment",
	    opt$method == "moment2" ~ "--moment2-beta"
	  )
		call <- switch(type,
	qcovar = sprintf("%6$s %5$s \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
%7$s %8$s \
%9$s %10$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform, opt$log, argname, opt$osca, 
	                 if(!is.null(opt$momentpercent)) "--moment-percent" else "",
	                 if(!is.null(opt$momentpercent)) opt$momentpercent else "",
	                 if(!is.null(opt$momentwind)) "--moment-wind" else "",
	                 if(!is.null(opt$momentwind)) opt$momentwind else ""
	                 ),
          	
	covar = sprintf("%6$s %5$s \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
%7$s %8$s \
%9$s %10$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform, opt$log, argname, opt$osca,
	                if(!is.null(opt$momentpercent)) "--moment-percent" else "",
	                if(!is.null(opt$momentpercent)) opt$momentpercent else "",
	                if(!is.null(opt$momentwind)) "--moment-wind" else "",
	                if(!is.null(opt$momentwind)) opt$momentwind else ""
	                ),
	covar_qcovar = sprintf("%6$s %5$s \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
%7$s %8$s \
%9$s %10$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform, opt$log, argname, opt$osca,
	                       if(!is.null(opt$momentpercent)) "--moment-percent" else "",
	                       if(!is.null(opt$momentpercent)) opt$momentpercent else "",
	                       if(!is.null(opt$momentwind)) "--moment-wind" else "",
	                       if(!is.null(opt$momentwind)) opt$momentwind else ""
	                       ),
	None = sprintf("%6$s %5$s \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--orm-cutoff %3$s \
%7$s %8$s \
%9$s %10$s \
--out %2$s > %4$s", opt$befile, opt$out, opt$cutofform, opt$log, argname, opt$osca,
	               if(!is.null(opt$momentpercent)) "--moment-percent" else "",
	               if(!is.null(opt$momentpercent)) opt$momentpercent else "",
	               if(!is.null(opt$momentwind)) "--moment-wind" else "",
	               if(!is.null(opt$momentwind)) opt$momentwind else ""
	               ),
	)
	call <- gsub("\n", " ", call)
	# Temporary fix -> for loco use v0.41
	# if(opt$method == "loco") {
	#   call <- gsub("~/osca_Linux", "~/osca_Linux0.41", call)
	# }
	system(call)
} else if(opt$method == "linear") {
  call <- switch(type,
                 qcovar = sprintf("%4$s --linear \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
 covar = sprintf("%4$s --linear \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
 covar_qcovar = sprintf("%4$s --linear \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
 None = sprintf("%4$s --linear \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
                 )
  call <- gsub("\n", " ", call)
  system(call)
} else if(opt$method == "logistic") {
  call <- switch(type,
                 qcovar = sprintf("%4$s --logistic \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
                 covar = sprintf("%4$s --logistic \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
                 covar_qcovar = sprintf("%4$s --logistic \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--covar %2$s_input/covar.txt \
--qcovar %2$s_input/qcovar.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
                 None = sprintf("%4$s --logistic \
--befile %1$s \
--pheno %2$s_input/phenotype.txt \
--exclude-probe %2$s_input/remove_probes.txt \
--out %2$s > %3$s", opt$befile, opt$out, opt$log, opt$osca),
  )
  call <- gsub("\n", " ", call)
  system(call)
}

cat("Done!\n")

# If method != Linear and !is.null(missingnessmatrix), add a 'NMISS column'
if(!opt$method %in% c("linear", "logistic") && !is.null(opt$missingness)) {
  filename <- dplyr::case_when(
    opt$method == "loco" ~ paste0(opt$out, ".loco.mlma"),
    opt$method %in% c("moment", "moment2") & opt$osca != "~/osca_Linux0.41" ~ paste0(opt$out, ".moment"),
    opt$method == "moment" & opt$osca == "/home/hers_en/phop/osca_Linux0.41" ~ paste0(opt$out, ".mlma"),
    opt$method == "moa" ~ paste0(opt$out, ".moa")
  )
  if(file.exists(filename)) {
    cat("Adding NMISS column")
    ewas <- readr::read_tsv(filename)
    missingness <- missingness[ewas$Probe,samplesheet$Sample_Name]
    ewas$NMISS <- nrow(samplesheet) - rowSums(missingness)
    readr::write_tsv(ewas, path = filename)
  }
}

cat(sprintf("\n\nAnalysis finished at: %s\n\n", Sys.time()))
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()

