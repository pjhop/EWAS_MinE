Scripts accompanying the manuscript 'Genome-wide study of DNA methylation in ALS implicates metabolic, inflammatory and cholesterol pathways'.  

**EWAS_main.Rmd**: All EWAS analyses (including sensitivity analyses), invokes several scripts in the **utils/** folder. 
Analyses are submitted on a slurm hpc.  
**EWAS_finetuningPCs.Rmd**: Fine-tuning the number of PCs for the LB algorithm.  
**EWAS_sensitivity.Rmd**: Plots sensitivity analyses + run leave-one-out analyses.  
**EWAS_plots.Rmd**: Main figure 1 of the manuscript.  
**EWASdb_enrichments.Rmd**: EWASdb enrichments + plots (including main figure 2).  
**getPMS.R**: Generate PMSs using published coefficients.  
**incremental_R2.Rmd**: PMS validation.  
**PMS_main.Rmd**: PMS case/con and survival analyses + plots (including main figure 3).  
**Survival_analyses.Rmd**: Survival analyses of DMPs.  

**utils/**: several wrapper scripts invoked in EWAS_main.Rmd.
