rm(list=ls())
gc()

library(data.table)
library(dplyr)

# Create phenofile to run analysis: early stopping (1-2-3 purchases) VS good adherers (at least 4 purch, at least 80% adherence)

covs <- fread('/home/ivm/drugs/data/finngen_R8_cov.txt')

drugs <- c("statins", "blood_pressure", "breast_cancer", "clopi_dipy", "doac", "glauc")

for (d in drugs){
  
  print(d)
  
  # Calculate discontinuation
  # Read full trajectories
  traj <- fread(paste0('/home/ivm/drugs/data/R8_',d,'_stopping.txt')) %>%
    left_join(covs, by = "FINNGENID")
  
  # Get IDs last purchase was at least 2 years before end of follow up
  to.keep <- traj %>%
    group_by(FINNGENID) %>%
    # keep last row per each group
    filter(row_number()==n()) %>% 
    filter(EVENT_AGE < (AGE_AT_DEATH_OR_END_OF_FOLLOWUP-2)) %>% 
    pull(FINNGENID)
  
  traj <- traj %>% 
    filter(FINNGENID %in% to.keep)

  min_age <- ifelse(d %in% c("breast_cancer", "glauc"), 0, 18)
  
  traj_sum <- traj %>% 
    group_by(FINNGENID) %>% 
    summarise(tot_purch = n(),
              tot_days = sum(days_next_purch, na.rm = T),
              age_first := first(EVENT_AGE)) %>% 
    filter(tot_purch < 4,
           age_first >= min_age)

  # "Good" adherers 
  adh <- fread(paste0('/home/ivm/drugs/data/R8_',d,'_summarized.txt')) %>% 
    select(FINNGENID, adherence, age_first = age_first_purch)
  
  good <- adh %>%
    filter(adherence >= 0.8) %>% 
    pull(FINNGENID)
  
  covs <- covs %>%
    left_join(traj_sum, by = "FINNGENID") %>%
    left_join(adh, by = "FINNGENID") %>% 
    mutate("{d}_stop_1" := case_when(tot_purch == 1 ~ 1,
                                     FINNGENID %in% good ~ 0),
           "{d}_stop_2" := case_when(tot_purch == 2 ~ 1,
                                     FINNGENID %in% good ~ 0),
           "{d}_stop_3" := case_when(tot_purch == 3 ~ 1,
                                     FINNGENID %in% good ~ 0),
           good = ifelse(FINNGENID %in% good, 1, 0)) %>% 
    mutate("age_first_{d}" := ifelse(good == 1, age_first.y, age_first.x),
           "age_first_{d}_stop_1" := ifelse(good == 1, age_first.y, age_first.x),
           "age_first_{d}_stop_2" := ifelse(good == 1, age_first.y, age_first.x),
           "age_first_{d}_stop_3" := ifelse(good == 1, age_first.y, age_first.x)) %>% 
    select(-age_first.x, -age_first.y, -tot_purch, -tot_days, -adherence, -good)
    
}

# add age_first squared
squared <- function(x) {
  return(x**2)
}

covs <- covs %>% 
  mutate(across(starts_with("age_first"), squared, .names="{.col}2"))

fwrite(covs, '/home/ivm/drugs/data/R8_cov_pheno_adherence_stop.txt', sep = '\t', quote = F)
fwrite(list(c(paste0(drugs,"_stop_1"),paste0(drugs,"_stop_2"),paste0(drugs,"_stop_3"))),
       '/home/ivm/drugs/data/R8_stop_phenolist.txt', col.names = F)

fwrite(covs, '/finngen/red/mcordioli/drugs-gwas/R8/R8_cov_pheno_adherence_stop.txt', sep = '\t', quote = F)

fwrite(list(c(paste0(drugs,"_stop_1"),paste0(drugs,"_stop_2"),paste0(drugs,"_stop_3"))),
       '/finngen/red/mcordioli/drugs-gwas/R8/R8_stop_phenolist.txt', col.names = F)

colnames(covs)

# count cases and controls
res <- NULL
for (d in drugs){
  for (i in c(1,2,3)){
    col <- paste0(d, "_stop_", i)
    
    res <- rbind(res, c(col,table(covs[[col]])["0"],table(covs[[col]])["1"]))
  }

  age.col <- paste0("age_first_",d)
  print( max(covs[[age.col]], na.rm = T) )
}
