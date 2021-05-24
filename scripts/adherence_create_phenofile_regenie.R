rm(list=ls())

library(data.table)
library(dplyr)

pheno_adherence_saige <- fread('/home/cordioli/drugs/data/R7_cov_pheno_adherence_150.txt')
pheno_adherence_comb_saige <- fread('/home/cordioli/drugs/data/R7_cov_pheno_combo_adherence.txt')

pheno_finngen_regenie <- fread('/home/cordioli/drugs/data/R7_COV_PHENO_V2.FID.txt.gz')

pheno_finngen_regenie <- pheno_finngen_regenie %>% 
  select(FID:DEATH)

pheno_adherence_saige <- pheno_adherence_saige %>% 
  select(IID = FINNGENID, statins:age_first_glauc)

pheno_adherence_regenie <- pheno_finngen_regenie %>% 
  left_join(pheno_adherence_saige, by = "IID")


pheno_adherence_comb_saige <- pheno_adherence_comb_saige %>% 
  select(IID = FINNGENID, combo_adherence_std)

pheno_adherence_comb_regenie <- pheno_finngen_regenie %>% 
  left_join(pheno_adherence_comb_saige, by = "IID")

fwrite(pheno_adherence_regenie, '/home/cordioli/drugs/data/R7_cov_pheno_adherence_regenie.txt.gz', sep = '\t', quote = F, na = "NA", compress = "gzip")
fwrite(pheno_adherence_comb_regenie, '/home/cordioli/drugs/data/R7_cov_pheno_combo_adherence_regenie.txt.gz', sep = '\t', quote = F, na = "NA", compress = "gzip")

system('gsutil -m cp /home/cordioli/drugs/data/R7_cov_pheno_adherence_regenie.txt.gz /home/cordioli/drugs/data/R7_cov_pheno_combo_adherence_regenie.txt.gz gs://mattia/drugs/')

colnames(pheno_adherence_regenie)

# count shared and non-shared missingness
m <- pheno_adherence_regenie %>% 
  filter( !is.na(statins) & !is.na(blood_pressure) & !is.na(clopi_dipy) & !is.na(breast_canc) & !is.na(glauc) & !is.na(doac))

sum(is.na(pheno_adherence_regenie$statins))
