rm(list=ls())
gc()

library(data.table)
library(dplyr)

covs <- fread('/home/ivm/drugs/data/finngen_R8_cov.txt')

# # # Statins
stt <- fread('/home/ivm/drugs/data/R8_statins_summarized.txt')

# Merge phenotypes with covariates for GWAS
ad <- stt %>%
  mutate(statins = adherence_std,
         age_first_statins = age_first_purch) %>%
  select(FINNGENID,statins,age_first_statins)

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence.txt', sep = '\t', quote = F)

# # # BP meds
bpp <- fread('/home/ivm/drugs/data/R8_blood_pressure_summarized.txt')

# # # 
# Merge phenotypes with covariates for GWAS
ad <- bpp %>%
  mutate(blood_pressure = adherence_std,
         age_first_blood_pressure = age_first_purch) %>%
  select(FINNGENID,blood_pressure,age_first_blood_pressure)

covs <- fread('/home/ivm/drugs/data/R8_cov_pheno_adherence.txt')

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence.txt', sep = '\t', quote = F)


# # # Clopi&asa+dip
cldi <- fread('/home/ivm/drugs/data/R8_clopi_dipy_summarized.txt')

# # # 
# Merge phenotypes with covariates for GWAS
ad <- cldi %>%
  mutate(clopi_dipy = adherence_std,
         age_first_clopi_dipy = age_first_purch) %>%
  select(FINNGENID,clopi_dipy,age_first_clopi_dipy)

covs <- fread('/home/ivm/drugs/data/R8_cov_pheno_adherence.txt')

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence.txt', sep = '\t', quote = F)


# # # Breast cancer
bcc <- fread('/home/ivm/drugs/data/R8_breast_cancer_summarized.txt')

# Merge phenotypes with covariates for GWAS
covs <- fread('/home/ivm/drugs/data/R8_cov_pheno_adherence.txt')

ad <- bcc %>%
  mutate(breast_canc = adherence_std,
         age_first_breast_canc = age_first_purch) %>%
  select(FINNGENID,breast_canc,age_first_breast_canc)

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence.txt', sep = '\t', quote = F)


# # # DOAC
do <- fread('/home/ivm/drugs/data/R8_doac_summarized.txt')

# # # 
# Merge phenotypes with covariates for GWAS
ad <- do %>%
  mutate(doac = adherence_std,
         age_first_doac = age_first_purch) %>%
  select(FINNGENID,doac,age_first_doac)

covs <- fread('/home/ivm/drugs/data/R8_cov_pheno_adherence.txt')

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence.txt', sep = '\t', quote = F)


# # # Glauc
gl <- fread('/home/ivm/drugs/data/R8_glauc_summarized.txt')

# # # 
# Merge phenotypes with covariates for GWAS
ad <- gl %>%
  mutate(glauc = adherence_std,
         age_first_glauc = age_first_purch) %>%
  select(FINNGENID,glauc,age_first_glauc)

covs <- fread('/home/ivm/drugs/data/R8_cov_pheno_adherence.txt')

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence.txt', sep = '\t', quote = F)


phenolist <- c('statins','blood_pressure','clopi_dipy','breast_canc', 'doac', 'glauc')
fwrite(list(phenolist), '/home/ivm/drugs/data/adherence_phenolist.txt', col.names = F)

phenolist <- c('glauc')
fwrite(list(phenolist), '/home/ivm/drugs/data/adherence_glauc_phenolist.txt', col.names = F)


# Create combined adherence phenotype
statins <- stt %>%
  select(FINNGENID, tot_purch, age_first_purch, adherence)

blood_pressure <- bpp %>%
  select(FINNGENID, tot_purch, age_first_purch, adherence)

clopi_dipy <- cldi %>%
  select(FINNGENID, tot_purch, age_first_purch, adherence)

breast <- bcc %>%
  select(FINNGENID, tot_purch, age_first_purch, adherence)

doac <- do %>%
  select(FINNGENID, tot_purch, age_first_purch, adherence)

glauc <- gl %>%
  select(FINNGENID, tot_purch, age_first_purch, adherence)

combo <- statins %>%
  full_join(blood_pressure, by = c("FINNGENID"), suffix = c('.statins', '.blood_pressure')) %>%
  full_join(clopi_dipy, by = c("FINNGENID"), suffix = c('.blood_pressure', '.clopi_dipy')) %>%
  full_join(breast, by = c("FINNGENID"), suffix = c('.clopi_dipy', '.breast')) %>%
  full_join(doac, by = c("FINNGENID"), suffix = c('.breast', '.doac')) %>%
  full_join(glauc, by = c("FINNGENID"), suffix = c('.doac', '.glauc'))

combo <- combo %>%
  mutate_at(
    vars(starts_with("adherence")),
    ~ case_when(is.na(.) ~ 0, TRUE ~ .))

combo <- combo %>%
  mutate_at(
    vars(starts_with("tot_purch")),
    ~ case_when(is.na(as.numeric(.)) ~ 0, TRUE ~ as.numeric(.)))

combo <- combo %>%
  mutate(combo_adherence = (tot_purch.statins*adherence.statins + tot_purch.blood_pressure*adherence.blood_pressure + 
                             tot_purch.clopi_dipy*adherence.clopi_dipy + tot_purch.breast*adherence.breast + 
                             tot_purch.doac*adherence.doac + tot_purch.glauc*adherence.glauc)/
                            (tot_purch.statins + tot_purch.blood_pressure + 
                            tot_purch.clopi_dipy + tot_purch.breast + 
                            tot_purch.doac + tot_purch.glauc))
combo <- combo %>%
  mutate(combo_adherence_std = as.numeric(scale(combo_adherence)))

fwrite(combo, '/home/ivm/drugs/data/R8_combo_adherence.txt', sep = '\t', quote = F)

ad <- combo %>%
  select(FINNGENID,combo_adherence_std)
  
covs <- fread('/home/ivm/drugs/data/finngen_R8_cov.txt')

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_combo_adherence.txt', sep = '\t', quote = F)

phenolist <- c('combo_adherence_std')
fwrite(list(phenolist), '/home/ivm/drugs/data/adherence_combo_phenolist.txt', col.names = F)