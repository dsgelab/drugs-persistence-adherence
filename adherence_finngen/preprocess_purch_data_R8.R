rm(list=ls())

library(R.utils)
library(data.table)
library(dplyr)

dat <- fread('/finngen/library-red/finngen_R8/phenotype_3.0/data/finngen_R8_detailed_longitudinal.txt.gz')

cov_pheno <- fread('/finngen/library-red/finngen_R8/analysis_covariates/finngen_R8_cov_1.0.txt.gz')

vnr <- fread('/home/ivm/drugs/data/vnr_mapping_6_5.tsv')

# Filter medication purchases
dat <- dat %>%
  filter(SOURCE == "PURCH") %>%
  # cast event_day to Date
  # cast CODE3 to integer for merging with VNRs
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY),
         vnr = as.integer(CODE3)) %>%
  # filter purchases from 1998 on
  filter(format(APPROX_EVENT_DAY,'%Y') >= 1998) %>%
  # filter individuals in the cov_pheno file, for which the GWAS will be run
  filter(FINNGENID %in% cov_pheno$FINNGENID) %>%
  # merge with VNRs
  left_join(vnr) %>%
  # discard not needed columns (ICD version/category)
  select(-ICDVER, -CATEGORY)

# Write to file and compress
fwrite(dat, '/home/ivm/drugs/data/finngen_R8_purch_vnr_98.gz', sep = '\t', quote = F, compress = "gzip", na = "NA")

print("* * * * * WRITTEN finngen_R8_purch_vnr_98.gz * * * * * ")

# Preprocess cov_pheno file - save only covariates
cov_pheno <- cov_pheno[,1:grep("PC20", colnames(cov_pheno))]
fwrite(cov_pheno, '/home/ivm/drugs/data/finngen_R8_cov.txt', sep = '\t', quote = F, na = "NA")

print("* * * * * WRITTEN finngen_R8_cov.txt * * * * * ")

dat <- fread('/home/ivm/drugs/data/finngen_R8_purch_vnr_98.gz')

# # # SOME BASIC STATS / QC

# individuals
length(unique(dat$FINNGENID))

# purchases
nrow(dat)

length(which(is.na(dat$CODE1)))
missing_atc <- dat %>% 
  filter(is.na(CODE1))

summary(missing_atc)

# n purchases with kela remb code 
length(which(!is.na(dat$CODE2)))


# Number of NA numb. packages
length(which(is.na(dat$CODE4)))

# Number packages = 0
length(which(dat$CODE4==0))