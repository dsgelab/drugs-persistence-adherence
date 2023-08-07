# Merge persistence files in a unique phenotype file for GWAS analyses

rm(list=ls())
gc()

library(data.table)
library(dplyr)

covs <- fread('/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz') %>% 
  select(FINNGENID = IID, everything())

drugs <- c("statins", "blood_pressure", "breast_cancer", "clopi_dipy", "doac")

for (d in drugs){
  
  print(d)
  
  # Calculate discontinuation
  # Read full trajectories
  p <- fread(paste0('/home/ivm/drugs/data/R10_',d,'_persistence.txt')) %>% 
    select(FINNGENID, age_first_purch, persistent)
  p[[d]] <- p$persistent
  p[[paste0("age_first_",d)]] <- p$age_first_purch
  p <- p %>% select(-age_first_purch, -persistent)
    
  covs <- covs %>%
    left_join(p, by = "FINNGENID")

}  

covs <- covs %>% select(FID, IID = FINNGENID, everything())
fwrite(covs, '/home/ivm/drugs/data/R10_cov_pheno_persistence.txt', sep = '\t', quote = F, na = "NA")
