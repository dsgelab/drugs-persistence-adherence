rm(list=ls())

library(dplyr)
library(data.table)

# Read finngen files / atc mapping

dat <- fread('/home/cordioli/R5_pheno/finngen_R5_v2_detailed_longitudinal.gz')
atc <- fread('/home/cordioli/drugs/ATC_translate_984.tsv')

atc <- atc %>%
  select('Id', 'Parent code', 'Long name', 'Longname')

# Select ongly drug purchases data and merge with vnr
dat <- subset(dat, SOURCE == "PURCH")

c <- dat %>%
  group_by(CODE1) %>%
  summarize(n())

cc <- merge(c,atc, by.x = "CODE1", by.y = "Id", all.x = T)

write.table(cc, '/home/cordioli/drugs/meds_count.tsv', sep = '\t', quote = F, row.names = F)
