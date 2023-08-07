# Merge adherence files in a unique phenotype file for GWAS analyses

rm(list=ls())
gc()

library(data.table)
library(dplyr)

covs <- fread('/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz') %>% 
  select(FINNGENID = IID, everything())

# # # Statins
stt <- fread('/home/ivm/drugs/data/R10_statins_summarized.txt')

# Merge phenotypes with covariates for GWAS
ad <- stt %>%
  mutate(statins = adherence_std,
         age_first_statins = age_first_purch) %>%
  select(FINNGENID,statins,age_first_statins)

covs <- covs %>%
  left_join(ad)


# # # BP meds
bpp <- fread('/home/ivm/drugs/data/R10_blood_pressure_summarized.txt')

# Merge phenotypes with covariates for GWAS
ad <- bpp %>%
  mutate(blood_pressure = adherence_std,
         age_first_blood_pressure = age_first_purch) %>%
  select(FINNGENID,blood_pressure,age_first_blood_pressure)

covs <- covs %>%
  left_join(ad)


# # # Clopi&asa+dip
cldi <- fread('/home/ivm/drugs/data/R10_clopi_dipy_summarized.txt')

# Merge phenotypes with covariates for GWAS
ad <- cldi %>%
  mutate(clopi_dipy = adherence_std,
         age_first_clopi_dipy = age_first_purch) %>%
  select(FINNGENID,clopi_dipy,age_first_clopi_dipy)

covs <- covs %>%
  left_join(ad)


# # # Breast cancer
bcc <- fread('/home/ivm/drugs/data/R10_breast_cancer_summarized.txt')

# Merge phenotypes with covariates for GWAS
ad <- bcc %>%
  mutate(breast_canc = adherence_std,
         age_first_breast_canc = age_first_purch) %>%
  select(FINNGENID,breast_canc,age_first_breast_canc)

covs <- covs %>%
  left_join(ad)

# # # DOAC
do <- fread('/home/ivm/drugs/data/R10_doac_summarized.txt')

# # # 
# Merge phenotypes with covariates for GWAS
ad <- do %>%
  mutate(doac = adherence_std,
         age_first_doac = age_first_purch) %>%
  select(FINNGENID,doac,age_first_doac)

covs <- covs %>%
  left_join(ad)

# Rename columns FID IID for REGENIE
covs <- covs %>% select(FID, IID = FINNGENID, everything())
fwrite(covs, '/home/ivm/drugs/data/R10_cov_pheno_adherence.txt', sep = '\t', quote = F, na = "NA")

#Covs in REGENIE json
c <- "SEX_IMPUTED,AGE_AT_DEATH_OR_END_OF_FOLLOWUP,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,IS_AFFY_V2,IS_AFFY_V2P2,IS_AFFY_V3,BATCH_DS1_BOTNIA_Dgi_norm,BATCH_DS10_FINRISK_Palotie_norm,BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm,BATCH_DS12_FINRISK_Summit_norm,BATCH_DS13_FINRISK_Bf_norm,BATCH_DS14_GENERISK_norm,BATCH_DS15_H2000_Broad_norm,BATCH_DS16_H2000_Fimm_norm,BATCH_DS17_H2000_Genmets_norm_relift,BATCH_DS18_MIGRAINE_1_norm_relift,BATCH_DS19_MIGRAINE_2_norm,BATCH_DS2_BOTNIA_T2dgo_norm,BATCH_DS20_SUPER_1_norm_relift,BATCH_DS21_SUPER_2_norm_relift,BATCH_DS22_TWINS_1_norm,BATCH_DS23_TWINS_2_norm_nosymmetric,BATCH_DS24_SUPER_3_norm,BATCH_DS25_BOTNIA_Regeneron_norm,BATCH_DS3_COROGENE_Sanger_norm,BATCH_DS4_FINRISK_Corogene_norm,BATCH_DS5_FINRISK_Engage_norm,BATCH_DS6_FINRISK_FR02_Broad_norm_relift,BATCH_DS7_FINRISK_FR12_norm,BATCH_DS8_FINRISK_Finpcga_norm,BATCH_DS9_FINRISK_Mrpred_norm,age_first_statins"
cc <- strsplit(c,",")[[1]]

setdiff(cc,colnames(covs))
colnames(covs)[startsWith(colnames(covs),"BATCH_DS6")]


phenolist <- c('statins','blood_pressure','clopi_dipy','breast_canc', 'doac')
fwrite(list(phenolist), '/home/ivm/drugs/data/R10_adherence_phenolist.txt', col.names = F)