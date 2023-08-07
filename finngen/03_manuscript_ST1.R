rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyr)
library(broom)

# Get N of index population
N <- nrow(fread('/home/ivm/drugs/data/R10_index_people.txt', header = F))

drugs <- c('statins', 'blood_pressure', 'clopi_dipy', 'breast_cancer', 'doac')

summary_tot <- NULL
IDs_tot <- NULL

for (drug in drugs) {
  
  adh <- fread(paste0('/home/ivm/drugs/data/R10_',drug,'_summarized.txt'))
  per <- fread(paste0('/home/ivm/drugs/data/R10_',drug,'_persistence.txt'))
  
  summary <- adh %>% 
    summarise(drug = drug,
              n_users = length(unique(c(adh$FINNGENID, per$FINNGENID))),
              prop_population = n_users/N,
              mean_adherence = mean(adherence),
              sd_adherence = sd(adherence),
              mean_treatment_years = mean(tot_days)/365.25,
              max_treatment_years = max(tot_days)/365.25,
              n_per = sum(per$persistent),
              n_non_per = sum(per$persistent == 0),
              prop_persistent = sum(per$persistent)/nrow(per),
              tot_purch = sum(tot_purch) + n_non_per)
  
  summary_tot <- bind_rows(summary_tot, summary)
  IDs_tot <- unique(c(IDs_tot, adh$FINNGENID, per$FINNGENID))
}

summary_tot <- summary_tot %>% 
  arrange(-n_users) %>%
  mutate(across(3:8, round, 3)) %>% 
  mutate(tot_unique_users = length(unique(IDs_tot)))

# N unique individuals included
length(unique(IDs_tot))
# [1] 190076

# N tot purchases
sum(summary_tot$tot_purch)
# [1] 10708878

fwrite(summary_tot, '/home/ivm/drugs/results/manuscript/table1_FG.tsv', sep = "\t")
fwrite(list(IDs_tot), '/home/ivm/drugs/data/R10_all_drugs_users.txt')

