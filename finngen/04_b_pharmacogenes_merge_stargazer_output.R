# # # Merge output files (chunks of 15k samples) from stargazer

rm(list=ls())

library(data.table)
library(dplyr)

# Merge chunk results

genes <- c("abcg2", "cyp2c9", "cyp2c19", "cyp2d6", "slco1b1")

for (g in genes){
  all <- NULL
  for (i in 1:29){
    d <- fread(paste0('/home/ivm/drugs/data/pharmacogenes/stargazer_call/',g,'/chunk_',i,'/report.tsv'))
    all <- bind_rows(all,d)
  }
  print(length(unique(all$Sample)))
  fwrite(all, paste0('/home/ivm/drugs/data/pharmacogenes/stargazer_call/',g,'_report.tsv'), sep = "\t")
}

# # # calculate frequency for each diplotype in FinnGen
tot <- NULL
for (g in genes){
  d <- fread(paste0('/home/ivm/drugs/data/pharmacogenes/stargazer_call/',g,'_report.tsv'))
  diplotypes <- sort(unique(d$Diplotype))
  
  for (di in diplotypes){
    d <- d %>% mutate(has_diplo = ifelse(Diplotype == di, 1, 0))
    
    su <- data.frame(gene = g,
                     diplotype = di,
                     count = sum(d$has_diplo),
                     freq = sum(d$has_diplo)/nrow(d))
    
    tot <- bind_rows(tot, su)
  }
}

fwrite(tot, '/home/ivm/drugs/data/pharmacogenes/stargazer_call/diplotype_freq.tsv', sep = "\t")
fwrite(tot, '/home/ivm/drugs/results/manuscript/pharmacogenes/diplotype_freq.tsv', sep = "\t")
