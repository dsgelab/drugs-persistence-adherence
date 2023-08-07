# # # Extract and pre-process drug purchase data

rm(list=ls())
gc()

library(R.utils)
library(data.table)
library(dplyr)
library(lubridate)

dat <- fread('/finngen/library-red/finngen_R10/phenotype_1.0/data/finngen_R10_detailed_longitudinal_1.0.txt.gz')

cov <- fread('/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz')

vnr <- fread('/home/ivm/drugs/data/vnr_mapping_6_5.tsv')

# Define index population: alive in 2010

# First calculate birth date: take first event in longitudinal file and
# calculate as aprox_event_day - event_age
bday <- dat %>%
  arrange(APPROX_EVENT_DAY) %>%
  distinct(FINNGENID, .keep_all = TRUE) %>% 
  mutate(bday = lubridate::as_date(APPROX_EVENT_DAY) - lubridate::dyears(EVENT_AGE) ) %>%
  select(FINNGENID, bday)

# Calculate death date
cov <- left_join(cov,bday, by = c("IID" = "FINNGENID")) %>%
  mutate(dday = bday + lubridate::dyears(DEATH_AGE))

alive.2010 <- cov %>% 
  filter(dday > "2010-01-01") %>% 
  select(FINNGENID = IID) %>% 
  pull(FINNGENID)

# Write list of index people
fwrite(list(alive.2010), '/home/ivm/drugs/data/R10_index_people.txt')

# Keep only medication purchases and index person
purch <- dat %>%
  filter(SOURCE == "PURCH",
         FINNGENID %in% alive.2010) %>%
  # cast event_day to Date
  # cast CODE3 to integer for merging with VNRs
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY),
         vnr = as.integer(CODE3)) %>%
  # filter purchases from 1998 on
  filter(format(APPROX_EVENT_DAY,'%Y') >= 1998) %>%
  # filter individuals in the cov_pheno file, for which the GWAS will be run
  filter(FINNGENID %in% cov$FID) %>%
  # merge with VNRs
  left_join(vnr) %>%
  # discard not needed columns (ICD version/category)
  select(-ICDVER, -CATEGORY)

# Write to file and compress
fwrite(purch, '/home/ivm/drugs/data/finngen_R10_purch_vnr_98_alive_2010.gz', sep = '\t', quote = F, compress = "gzip", na = "NA")
print("* * * * * WRITTEN finngen_R10_purch_vnr_98.gz * * * * * ")



dat <- fread('/home/ivm/drugs/data/finngen_R10_purch_vnr_98_alive_2010.gz')

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

# individuals
length(unique(dat$FINNGENID))
# 338183
# R10, index person only: 392428

# purchases
nrow(dat)
# [1] 54140783
# R7: 63084035
# R8: 68826654
# R10: 92806680

sum(is.na(dat$CODE1))
# [1] 0
# R8: 72012
# R10: 93954

missing_atc <- dat %>% 
  filter(is.na(CODE1))

summary(missing_atc)

# n purchases with kela remb code 
length(which(!is.na(dat$CODE2)))
# [1] 16513633
# R7: 19053983
# R8: 20654288


# Number of NA numb. packages
length(which(is.na(dat$CODE4)))
# [1] 0


# Number packages = 0
length(which(dat$CODE4==0))
# 1126
# R7: 1325
# R8: 1454
# comes from the registry like this