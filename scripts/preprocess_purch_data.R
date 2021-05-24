rm(list=ls())

library(R.utils)
library(data.table)
library(dplyr)

# Mount finngen data bucket, 
# Read in finngen files and VNRs mapping file

# cov_pheno and endpoint file:
# gs://finngen-production-library-red/finngen_R7/phenotype_4.0/data/finngen_R7_cov_1.0.txt.gz
# gs://finngen-production-library-red/finngen_R7/phenotype_4.0/data/finngen_R7_endpoint.txt.gz

# detailed longitudinal file:
# gs://finngen-production-library-red/finngen_R7/phenotype_2.0/data/finngen_R7_detailed_longitudinal.txt.gz

system("gcsfuse --only-dir finngen_R7/phenotype_2.0/data/ --file-mode 444 finngen-production-library-red /home/cordioli/mount_fg/")
dat <- fread('/home/cordioli/mount_fg/finngen_R7_detailed_longitudinal.txt.gz')
system("")

print("* * * * * LOADED finngen_R7_detailed_longitudinal.txt.gz * * * * * ")

system("gcsfuse --only-dir finngen_R7/phenotype_4.0/data/ --file-mode 444 finngen-production-library-red /home/cordioli/mount_fg/")
cov_pheno <- fread('/home/cordioli/mount_fg/finngen_R7_cov_1.0.txt.gz')
system("fusermount -u /home/cordioli/mount_fg")

print("* * * * * LOADED finngen_R7_cov_1.0.txt.gz * * * * * ")

vnr <- fread('/home/cordioli/drugs/data/vnr_mapping_6_5.tsv')

print("* * * * * LOADED vnr_mapping_6_5.tsv * * * * * ")

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
fwrite(dat, '/home/cordioli/drugs/data/finngen_R7_purch_vnr_98.gz', sep = '\t', quote = F, compress = "gzip")

print("* * * * * WRITTEN finngen_R7_purch_vnr_98.gz * * * * * ")

# Preprocess cov_pheno file - save only covariates
cov_pheno <- cov_pheno[,1:grep("PC20", colnames(cov_pheno))]
fwrite(cov_pheno, '/home/cordioli/drugs/data/finngen_R7_cov.txt', sep = '\t', quote = F)

print("* * * * * WRITTEN finngen_R7_cov.txt * * * * * ")

# dat <- fread('/home/cordioli/drugs/data/finngen_R7_purch_vnr_98.gz')

# # # DATA DESCRIPTION
#
#   Field |Register        |Register long name         |
#   ------|----------------|---------------------------|
#   SOURCE|PURCH           |Kela drug purchase register|
#   
#   Field |Code in register|Code description.         |
#   ------|---------------|--------------------------|
#   CODE1 | ATC_CODE      | ATC code                  |
#   CODE2 | SAIR 	        | Kela reimbursement code   |
#   CODE3 | VNRO	        | Product number		        |
#   CODE4 | PLKM          | Number of packages        |

# VNR mapping
# vahvuus: strength, mg/pill | mg-microg/ml | international units (IU)
# pkoko: size of the package, e.g. “100 FOL” 100 pills in foil package

# # # SOME BASIC STATS / QC

# purchases
nrow(dat)
# [1] 54140783
# R7: 63084035

length(which(is.na(dat$CODE1)))
# [1] 0

# n purchases with kela remb code 
length(which(!is.na(dat$CODE2)))
# [1] 16513633
# R7: 19053983

# Number of NA numb. packages
length(which(is.na(dat$CODE4)))
# [1] 0

# Number packages = 0
length(which(dat$CODE4==0))
# 1126
# R7: 1325
# comes from the registry like this