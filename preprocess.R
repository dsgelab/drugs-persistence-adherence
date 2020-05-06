rm(list=ls())

library(R.utils)
library(data.table)
library(dplyr)

# Read finngen longitudinal file, merge with VNR mapping, save only purchases
dat <- fread('/home/cordioli/R5_pheno/finngen_R5_v3_detailed_longitudinal.gz')
vnr <- fread('/home/cordioli/drugs/data/vnr_mapping_6_5.tsv')

# Select medication purchases data from 1998 and merge with vnr, write to disk
dat <- dat %>%
  filter(SOURCE == "PURCH") %>%
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY),
         vnr = as.integer(CODE3)) %>%
  filter(format(APPROX_EVENT_DAY,'%Y') >= 1998) %>%
  left_join(vnr) %>%
  select(-ICDVER, -CATEGORY)

# Write to file and compress
fwrite(dat, '/home/cordioli/drugs/data/finngen_R5_v3_purch_vnr', sep = '\t', quote = F)
gzip('/home/cordioli/drugs/data/finngen_R5_v3_purch_vnr', destname='/home/cordioli/drugs/data/R5_v3_purch_vnr_98.gz')

# Veeeeery weird behaviour in merging VNRs wtf
# v_map <- unique(vnr$vnr) #integer
# v_dat <- unique(dat$CODE3) #character
# v_map_chr <- as.character(v_map)
# length(intersect(v_map_chr, v_dat))
# # 7202
# v_dat_int <- as.integer(v_dat)
# length(intersect(v_map, v_dat_int))
# # 12564

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