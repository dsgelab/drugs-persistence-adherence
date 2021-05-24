rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

source('/home/cordioli/drugs/adherence_funs.R')

write.flag <- T

purch <- fread('/home/cordioli/drugs/data/finngen_R7_purch_vnr_98.gz')
covs <- fread('/home/cordioli/drugs/data/finngen_R7_cov.txt')


aspirin <- purch %>% 
  filter(grepl("^B01AC06", CODE1))

length(unique(aspirin$FINNGENID))

N <- aspirin %>% 
  group_by(FINNGENID) %>% 
  summarise(n = n())

summary(N$n)
