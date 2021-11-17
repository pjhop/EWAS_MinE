## Libraries
library(dplyr)
library(magrittr)

## beta filepaths
betas <- c(
  "",
  "",
  "",
  ""
)
names(betas) <- c("mine_450k", "mine_epic", "AUS_batch1", "AUS_batch2")

## McCartney et al. predictors 
McCartney_bmi <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=1)
McCartney_smoking <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=2)
McCartney_alc <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=3)
McCartney_education <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=4)
McCartney_totalchol <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=5)
McCartney_HDL <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=6)
McCartney_LDL <- readxl::read_excel("../data/13059_2018_1514_MOESM1_ESM.xlsx", sheet=7)

## BMI (Hamilton et al.)
bmi <- readxl::read_excel("../data/41366_2018_262_MOESM1_ESM.xlsx",
                          sheet = 8, skip = 1)
## Gadd et al.
coefs_Gadd <- readxl::read_excel("../data/media-5.xlsx", sheet = 3, skip=2)

## Liu alcohol scores
library(dnamalci)

## CRP
crp <- tibble(
  Probe = c("cg10636246", "cg17501210", "cg18608055",
            "cg03957124", "cg04987734", "cg04523589", 
            "cg17980786", "cg02341197"),
  coef = c(-0.0069, -0.0065, -0.0043, -0.0030,
           0.0041, 0.0022, 0.0026, 0.0030
  )
)

## Function to get the scores
get_scores <- function(study) {
  print(study)
  beta <- readRDS(betas[study])
  beta[is.na(beta)] <- 0
  
  ## McCartney et al. 
  McCartney_bmi_ <- McCartney_bmi %>% dplyr::filter(CpG %in% rownames(beta))
  McCartney_smoking_ <- McCartney_smoking %>% dplyr::filter(CpG %in% rownames(beta))
  McCartney_alc_ <- McCartney_alc %>% dplyr::filter(CpG %in% rownames(beta))
  McCartney_education_ <- McCartney_education %>% dplyr::filter(CpG %in% rownames(beta))
  McCartney_totalchol_ <- McCartney_totalchol %>% dplyr::filter(CpG %in% rownames(beta))
  McCartney_HDL_ <- McCartney_HDL %>% dplyr::filter(CpG %in% rownames(beta))
  McCartney_LDL_ <- McCartney_LDL %>% dplyr::filter(CpG %in% rownames(beta))
  
  scores_McCartney <- tibble(
    Sample_Name = colnames(beta),
    BMI_McCartney = (t(beta[McCartney_bmi_$CpG,]) %*% McCartney_bmi_$Beta)[,1],
    smoking_McCartney = (t(beta[McCartney_smoking_$CpG,]) %*% McCartney_smoking_$Beta)[,1],
    alcohol_McCartney = (t(beta[McCartney_alc_$CpG,]) %*% McCartney_alc_$Beta)[,1],
    education_McCartney = (t(beta[McCartney_education_$CpG,]) %*% McCartney_education_$Beta)[,1],
    totalchol_McCartney = (t(beta[McCartney_totalchol_$CpG,]) %*% McCartney_totalchol_$Beta)[,1],
    HDLchol_McCartney = (t(beta[McCartney_HDL_$CpG,]) %*% McCartney_HDL_$Beta)[,1],
    LDLchol_McCartney = (t(beta[McCartney_LDL_$CpG,]) %*% McCartney_LDL_$Beta)[,1]
  )
  
  ## BMI 
  intercept_bmi <- bmi %>% filter(CpG == "(Intercept)") %$% Coefficient 
  bmi_ <- bmi %>%
    dplyr::filter(CpG %in% rownames(beta))
  BMI_Hamilton <- tibble(
    Sample_Name = colnames(beta),
    BMI_Hamilton = intercept_bmi+((t(beta[bmi_$CpG,]) %*% bmi_$Coefficient)[,1])
  )
    
  ## Alcohol 
  models <- dnamalci.models()
  lst <- list()
  for(model in models) {
    lst[[model]] <- dnamalci(beta, model = model)
    lst[[model]]$Sample_Name <- colnames(beta)
  }
  alcohol_Liu <- tibble(
    Sample_Name = lst$dnamalc.5cpg$Sample_Name,
    dnamalc.5cpg_Liu = lst$dnamalc.5cpg$score,
    dnamalc.23cpg_Liu = lst$dnamalc.23cpg$score,
    dnamalc.78cpg_Liu =   lst$dnamalc.78cpg$score,
    dnamalc.144cpg_Liu  = lst$dnamalc.144cpg$score
  )
  
  ## CRP
  crp_ <- crp %>% dplyr::filter(Probe %in% rownames(beta))
  CRP <- t(beta[crp_$Probe,]) %*% crp_$coef
  CRP <- tibble(Sample_Name = rownames(CRP), CRP_Ligthart = CRP[,1])
  
  # Gadd et al. predictors
  get_score <- function(protein, beta, coefs) {
    coef <- coefs %>% dplyr::filter(`Proxy protein` == protein)
    beta_ <- beta[rownames(beta) %in% coef$`CpG Site`,,drop=FALSE]
    coef <- coef %>% dplyr::filter(`CpG Site` %in% rownames(beta_))
    score <- tibble(Sample_Name = colnames(beta_))
    score[[protein]] <- (t(beta_[coef$`CpG Site`,]) %*% coef$Coefficient)[,1]
    score
  }
  scores_Gadd <- purrr::map(unique(coefs_Gadd$`Proxy protein`),
                       .f = get_score,
                       beta = beta,
                       coefs = coefs_Gadd)
  scores_Gadd_joined <- scores_Gadd[[1]] 
  for(i in 2:length(scores_Gadd)) {
    scores_Gadd_joined <- scores_Gadd_joined %>%
      dplyr::left_join(scores_Gadd[[i]],by="Sample_Name")
  }
  
  
  # Combine all
  scores_all <- scores_McCartney %>%
    dplyr::left_join(BMI_Hamilton, by = "Sample_Name") %>%
    dplyr::left_join(alcohol_Liu, by = "Sample_Name") %>% 
    dplyr::left_join(CRP, by = "Sample_Name") %>%
    dplyr::left_join(scores_Gadd_joined, by = "Sample_Name") 
  scores_all
}

predictions <- purrr::map_df(names(betas), .f = get_scores)

saveRDS(predictions, file = "")