rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyr)

# Load files with the predictors we want to use for the model
# Minimal phenotype file:
# - sex
# - year of birth
# - mother tongue
demo <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv', select = c("FINREGISTRYID", "date_of_birth", "sex", "mother_tongue")) %>% 
    mutate(date_of_birth = as.Date(date_of_birth, format = "%d-%m-%Y"),
           yob = as.integer(format(date_of_birth, "%Y")),
           mother_tongue_fin_sv = ifelse(mother_tongue %in% c("fi","sv"), 1, 0)) %>% 
    select(HETU = FINREGISTRYID, date_of_birth, yob, sex, mother_tongue_fin_sv)

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
drug <- 'statins'

for (drug in drugs) {
  # Load all predictors for adherence
  cci <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.CCI.',drug,'.rds'))
  prev <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.prevention.',drug,'.rds'))
  edu <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.edu.',drug,'.rds'))
  liv <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.living.',drug,'.rds'))
  soc <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.social.',drug,'.rds'))
  
  d <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
    left_join(demo, by = "HETU") %>% 
    left_join(cci, by = "HETU") %>% 
    left_join(prev, by = "HETU") %>%
    left_join(edu, by = "HETU") %>% 
    left_join(liv, by = "HETU") %>%
    left_join(soc, by = "HETU") %>% 
    mutate(age_first_purch_scaled = as.numeric(scale(age_first_purch)),
           CCI_score_scaled = as.numeric(scale(CCI_score)))
  
  saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors_combined.adherence.',drug,'.rds'))
  
  if (drug == "blood_pressure"){
    d <- d %>% select(-secondary)
  }
  
  d_complete <- d %>% 
    filter(complete.cases(.))
    
  saveRDS(d_complete, paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.adherence.',drug,'.rds'))

  # Load all predictors for persistence
  cci <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.CCI.',drug,'.rds'))
  prev <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.prevention.',drug,'.rds'))
  edu <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.edu.',drug,'.rds'))
  liv <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.living.',drug,'.rds'))
  soc <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.social.',drug,'.rds'))
  
  d <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
    left_join(demo, by = "HETU") %>% 
    left_join(cci, by = "HETU") %>% 
    left_join(prev, by = "HETU") %>%
    left_join(edu, by = "HETU") %>% 
    left_join(liv, by = "HETU") %>%
    left_join(soc, by = "HETU") %>% 
    mutate(age_first_purch_scaled = as.numeric(scale(age_first_purch)),
           CCI_score_scaled = as.numeric(scale(CCI_score)))
  
  saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors_combined.persistence.',drug,'.rds'))
  
  if (drug == "blood_pressure"){
    d <- d %>% select(-secondary)
  }

  d_complete <- d %>% 
    filter(complete.cases(.))
  
  saveRDS(d_complete, paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.persistence.',drug,'.rds'))
  
}