## Perform linear regression using limma

## Libraries ---------------------------

library(limma)
library(tidyverse)
library(magrittr)
library(optparse)

## Parse parameters ---------------------------

option_list <- list(
  make_option(c("-b", "--beta"), action="store", type = "character", default=NULL,
              help="Path to beta matrix (RDS file)"),
  make_option(c("-s", "--samplesheet"), action="store", type = "character", default=NULL,
              help="Path to samplesheet (RDS file)"),
  make_option(c("-r", "--removesamples"), action="store", type = "character", default=NULL,
              help="Txt file containing samples that should be removed"),
  make_option(c("--keepsamples"), action="store", type = "character", default=NULL,
              help="Txt file containing samples that should be kept"),
  make_option(c("-a", "--removeprobes"), action="store", type = "character", default=NULL,
              help="Txt file containing probes that should be removed"),
  make_option(c("-x", "--removexy"), action="store", type = "character", default="y",
              help="Remove xy probes? (y/n)"),
  make_option(c("-p", "--pheno"), action="store", type = "character", default=NULL,
              help="Name of column that contains phenotype, should be present in samplesheet"),
  make_option(c("-l", "--levels"), action="store", type = "character", default=NULL,
              help="If phenotype is binary, specifiy the order of levels. Seperate levels by a comma"),
  make_option(c("-c", "--covariates"), action="store", type = "character", default=NULL,
              help="Qualitative ovariates to include, should be present in samplesheet"),
  make_option(c("--qcovariates"), action="store", type = "character", default=NULL,
              help="Quantiative ovariates to include, should be present in samplesheet"),
  make_option(c("--array"), action="store", type = "character", default=NULL,
              help="Which array is used? (450k/EPIC)"),
  make_option(c("--missingness"), action="store", type = "character", default=NULL,
              help="Sparse matrix containing missingness"),
  make_option(c("--missthreshold"), action="store", type = "numeric", default=1,
              help="Probes with missingness > threshold will be removed"),
  make_option(c("-q", "--qqplot"), action="store", type = "character", default="n",
              help="Whether to make a qqplot (y/n)"),
  make_option(c("-o", "--out"), action="store", type = "character", default=NULL,
              help="Filename of output")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Print options
cat("Following options are specified:\n\n")
for(i in names(opt)) {
  if(i != "help") {
    cat(sprintf("--%s %s\n", i, opt[[i]]))
  }
}
cat("\n")

## Load data -------------------------------------------------------------------
beta <- readRDS(opt$beta)
samplesheet <- as.data.frame(readRDS(opt$samplesheet))
if(!is.null(opt$covariates)) {
  covariates <- str_trim(unlist(str_split(opt$covariates, pattern = ",")))
} else {
  covariates <- NULL
}

if(!is.null(opt$qcovariates)) {
  qcovariates <- str_trim(unlist(str_split(opt$qcovariates, pattern = ",")))
} else {
  qcovariates <- NULL
}

pheno <- opt$pheno
levels <- str_trim(unlist(str_split(opt$levels, pattern = ",")))

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
  if(opt$array == "450k") {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    anno <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
    xy <- anno %>% filter(chr %in% c("chrX", "chrY")) %$% Name
  } else {
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    anno <- data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
    xy <- anno %>% filter(chr %in% c("chrX", "chrY")) %$% Name
  }
}
  
# Check if phenotype and covariates are present in samplesheet
if(!is.null(covariates) && !all(covariates %in% colnames(samplesheet))) {
  stop("Not all covariates are present in samplesheet")
}

if(!is.null(qcovariates) && !all(qcovariates %in% colnames(samplesheet))) {
  stop("Not all qcovariates are present in samplesheet")
}

if(!(pheno %in% colnames(samplesheet))) {
  stop("Phenotype is not present in samplesheet")
}

if(!all(colnames(beta) %in% samplesheet$Sample_Name)) {
  cat("Not all samples in DNAm data are present in samplesheet, subsetting the beta-matrix")
  ncol1 <- ncol(beta)
  beta <- beta[,colnames(beta) %in% as.character(samplesheet$Sample_Name)]
  ncol2 <- ncol(beta)
  cat(sprintf("Removed %s out of %s samples from the beta-matrix", ncol1 - ncol2, ncol1))
}

if(!all(samplesheet$Sample_Name %in% colnames(beta))) {
  cat("Not all samples in DNAm data are present in samplesheet, subsetting the samplesheet")
  ncol1 <- nrow(samplesheet)
  samplesheet <- samplesheet %>% dplyr::filter(Sample_Name %in% colnames(beta))
  ncol2 <- nrow(samplesheet)
  cat(sprintf("Removed %s out of %s samples from the samplesheet", ncol1 - ncol2, ncol1))
}

## Filter data ---------------------------

## Samples
cat("Samples----------\n\n")
cat(sprintf("Removing %s samples specified by user.\n", length(remove_samples)))
samplesheet <- samplesheet[!(samplesheet$Sample_Name %in% remove_samples),]

if(!is.null(opt$keepsamples)) {
  cat(sprintf("Keeping %s samples as specified by user.\n", sum(keep_samples %in% samplesheet$Sample_Name)))
  samplesheet <- samplesheet[(samplesheet$Sample_Name %in% keep_samples),]
}

if(!is.null(opt$levels)) {
  cat(sprintf("%s will be coded as 1, and %s will be coded as 0.\n", levels[2], levels[1]))
  cat(sprintf("%s column contains %s sample(s) that have missing data for this variable, these samples will be removed.\n",
              pheno, sum(is.na(samplesheet[,pheno]))))
  cat(sprintf("%s column contains %s samples with level(s) that are different than the specified levels, these samples will be removed.\n",
              pheno, sum(!(samplesheet[,pheno] %in% levels))))
  samplesheet <- samplesheet[!is.na(samplesheet[,pheno]) & (samplesheet[,pheno] %in% levels),]
  samplesheet[,pheno] <- factor(samplesheet[,pheno], levels = levels)
} else {
  cat(sprintf("%s column contains %s sample(s) that have missing data for this variable, these samples will be removed.\n",
              pheno, sum(is.na(samplesheet[,pheno]))))
  samplesheet <- samplesheet[!is.na(samplesheet[,pheno]),]
  
}

# Covariates - remove missings
if(!is.null(covariates)) {
  cat(sprintf("%s samples have missing data for one of the covariates, these samples will be removed.\n\n",
              sum(!complete.cases(samplesheet[,covariates]))))
  samplesheet <- samplesheet[complete.cases(samplesheet[,covariates]),]  
}

if(!is.null(qcovariates)) {
  cat(sprintf("%s samples have missing data for one of the qcovariates, these samples will be removed.\n\n",
              sum(!complete.cases(samplesheet[,qcovariates]))))
  samplesheet <- samplesheet[complete.cases(samplesheet[,qcovariates]),]
}


## Remove probes based on missingness / user-specified
if(!is.null(opt$missingness)) {
  missingness <- missingness[ ,samplesheet$Sample_Name]
  prop_missingness <- Matrix::rowMeans(missingness)
  prop_missingness <- names(prop_missingness[prop_missingness > opt$missthreshold])
  cat(sprintf("%s probes with missingness > %s will be removed\n", sum(rownames(beta) %in% prop_missingness), opt$missthreshold))
} else {
  prop_missingness <- c()
}

cat("Probes----------\n\n")
cat(sprintf("Removing %s probes (missingness+user-specified probes)\n", sum(c(remove_probes, prop_missingness) %in% rownames(beta))))
beta  <- beta[!(rownames(beta) %in% c(remove_probes, prop_missingness)),, drop = FALSE]

if(opt$removexy == "y") {
  cat(sprintf("Removing %s probes located on the sex chromosomes.\n", sum(rownames(beta) %in% xy)))
  beta  <- beta[!(rownames(beta) %in% xy),, drop = FALSE]
}

# Match beta file
beta <- beta[,as.character(samplesheet$Sample_Name)]

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
    covariates <- covariates[!covariates %in% names(nr.levels)]
    if(length(covariates) == 0) {
      cat("\nThere are no qualitative covariates left after removing covariates with only one level.\n")
      covariates <- NULL
    } 
  } 
}

## Model matrix ---------------------------

if (!is.null(covariates) && !is.null(qcovariates)) {
  formula <- sprintf("~%s + %s", pheno, paste(c(covariates, qcovariates), collapse = " + "))
  cat(sprintf("\n The following regression is used:\n%s\n\n", formula))
  design <- model.matrix(as.formula(formula), data = samplesheet)
} else if(!is.null(covariates)) {
  formula <- sprintf("~%s + %s", pheno, paste(c(covariates), collapse = " + "))
  cat(sprintf("\n The following regression is used:\n%s\n\n", formula))
  design <- model.matrix(as.formula(formula), data = samplesheet)
} else if(!is.null(qcovariates)) {
  formula <- sprintf("~%s + %s", pheno, paste(c(qcovariates), collapse = " + "))
  cat(sprintf("\n The following regression is used:\n%s\n\n", formula))
  design <- model.matrix(as.formula(formula), data = samplesheet)
} else {
  formula <- sprintf("~%s", pheno)
  cat(sprintf("\n The following regression is used:\n%s\n\n", formula))
  design <- model.matrix(as.formula(formula), data = samplesheet)
}

## lmFit ---------------------------

cat("Fitting..\n\n")
fit <- lmFit(beta, design)
cat("Done!\n\n")

## Test statistics ---------------------------

ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
pval <- 2*pt(-abs(ordinary.t), fit$df.residual)
pval <- pval[,2]

stats <- tibble(Probe = names(pval),
                b = fit$coefficients[,2][as.character(names(pval))],
                t = ordinary.t[,2],
                se = b / t,
                p = pval, 
                logp = -log10(pval))

# Calculate lambda:
chisq <- qchisq(1 - stats$p,1)
inflation <- median(chisq)/qchisq(0.5,1)
cat(sprintf("Lambda: %s\n\n", inflation))

if(opt$qqplot == "y") {
  # QQ-plot
  png(sprintf("%s_qqplot.png", opt$out), width = 7, height = 7, units = "in", res = 900)
  qqman::qq(stats$p)
  maxlogp <- -log10(0.05/nrow(stats))
  adjust <- maxlogp * 0.08
  inflation2 <- signif(inflation,3)
  text(adjust, maxlogp - adjust, bquote(lambda == .(inflation2)))
  dev.off()
  
  # Make plot with truncated logp-values
  png(sprintf("%s_qqplot_truncated.png", opt$out) , width = 7, height = 7, unit = "in", res = 900)
  stats_truncated <- stats %>% mutate(p = ifelse(p < 0.05/nrow(.), 0.05/nrow(.), p))
  qqman::qq(stats_truncated$p)
  maxlogp <- -log10(0.05/nrow(stats))
  adjust <- maxlogp * 0.08
  inflation2 <- signif(inflation,3)
  text(adjust, maxlogp - adjust, bquote(lambda == .(inflation2)))
  dev.off()
}

## Add NMISS and annotation
stats$NMISS <- rowSums(!is.na(beta))

if(opt$array == "450k" && opt$removexy != "y") {
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  anno <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
} else if(opt$array == "EPIC" && opt$removexy != "y") {
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  anno <- data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
}

stats <- stats %>%
  dplyr::left_join(anno[,c("Name", "chr", "pos", "strand")], by = c("Probe" = "Name")) %>%
  dplyr::rename(
    Chr = chr,
    bp = pos,
    Orientation = strand
  ) %>%
  dplyr::mutate(Gene = NA) %>%
  dplyr::select(Chr, Probe, bp, Gene, Orientation, b, se, t, p, logp, NMISS)

## Save data -------------------------------------------------------------------

saveRDS(fit, file = sprintf("%s_fit.rds", opt$out))
write_tsv(stats, path = sprintf("%s_stats.tsv", opt$out))