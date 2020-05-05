rm(list=ls())

library(R.utils)
library(data.table)
library(dplyr)

# Read finngen longitudinal file, merge with VNR mapping, save only purchases
dat <- fread('/home/cordioli/R5_pheno/finngen_R5_v3_detailed_longitudinal.gz')
vnr <- fread('/home/cordioli/drugs/vnr_mapping.tsv')
vnr$CODE3 <- as.character(vnr$vnr)

dat <- dat %>%
  filter(SOURCE == "PURCH")

# Select only drug purchases data from 1998 and merge with vnr, write to disk
dat <- dat %>%
  filter(SOURCE == "PURCH") %>%
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
  filter(format(APPROX_EVENT_DAY,'%Y') >= 1998) %>%
  left_join(vnr) %>%
  select(-ICDVER, -CATEGORY, -ATC, -vnr)

# # Write to file and compress
# fwrite(dat, '/home/cordioli/drugs/finngen_R5_v3_purch_vnr', sep = '\t', quote = F)
# gzip('/home/cordioli/drugs/finngen_R5_v3_purch_vnr', destname='/home/cordioli/drugs/R5_v3_purch_vnr_98.gz')


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
# nrow(dat)
# [1] 42603655

# VNRs in the data
vnr_dat <- dat %>%
  select(CODE3, valmiste, vahvuus, vahvuus_num, pkoko, pkoko_num) %>%
  distinct()

# VNRs with missing dosage and size
# length(which(is.na(vnr_dat$vahvuus_num) & is.na(vnr_dat$pkoko_num)))
# [1] 6870
# VNRs with missing size
# length(which(is.na(vnr_dat$pkoko_num)))
# [1] 6871
# VNRs with missing dosage
# length(which(is.na(vnr_dat$vahvuus_num)))
# [1] 7778

# Are we missing old VNRs ?
vnr_missing_years <- dat %>%
  select(CODE3, APPROX_EVENT_DAY, pkoko_num) %>%
  mutate(year = format(APPROX_EVENT_DAY, "%Y")) %>%
  select(-APPROX_EVENT_DAY) %>%
  distinct() %>%
  group_by(year) %>%
  summarise(TOT = n(),
            nNA = length(which(is.na(pkoko_num))),
            perc_NA = NNA/TOT)
# year    TOT   nNA perc_NA
# <chr> <int> <int>   <dbl>
# 1 1998   3844  1046   0.272
# 2 1999   3906  1059   0.271
# 3 2000   3811  1042   0.273
# 4 2001   3784  1144   0.302
# 5 2002   3873  1377   0.356
# 6 2003   3989  1634   0.410
# 7 2004   4011  1898   0.473
# 8 2005   4156  2239   0.539
# 9 2006   4179  2474   0.592
# 0 2007   4211  2593   0.616
# 1 2008   4441  2696   0.607
# 2 2009   4609  2635   0.572
# 3 2010   4753  2796   0.588
# 4 2011   4768  2714   0.569
# 5 2012   4784  2558   0.535
# 6 2013   4712  2382   0.506
# 7 2014   4716  2274   0.482
# 8 2015   4813  2209   0.459
# 9 2016   4918  2190   0.445
# 0 2017   5147  2237   0.435
# 1 2018   5345  2308   0.432
# 2 2019   2940  1239   0.421

# n purchases with ATC code NA
# length(which(is.na(dat$CODE1)))
# [1] 43970
# dat <- dat %>%
#   filter(SOURCE == "PURCH")
# length(which(is.na(dat$CODE1)))
# [1] 48761 IN TOTAL

# n purchases with kela remb code 
# length(which(!is.na(dat$CODE2)))
# [1] 13059311

# Number of NA numb. packages
# length(which(is.na(dat$CODE4)))
# [1] 0

# Number packages = 0
# length(which(dat$CODE4==0))
# 958
# comes from the registry like this