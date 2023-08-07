rm(list=ls())
gc()

library(ICCI)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

source('/data/projects/project_adherence/scripts/00_adherence_funs_finregistry.R')

files <- list.files('/data/projects/project_adherence/data/longitudinal/', pattern = 'all.drugs.*')
long <- rbindlist(lapply(files, function(name) readRDS(paste0('/data/projects/project_adherence/data/longitudinal/', name))))

print('Loaded long file')

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac', 'glaucoma')

define_CCI <- function(long, drug) {
  print(drug)
  
  # Adherence
  adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
    select(HETU, age_first_purch)
  
  long2 <- inner_join(long, adh, by = c("FINREGISTRYID" = "HETU")) %>% 
    mutate(ICDVER = ifelse(ICDVER == "O3", "10", ICDVER)) %>% 
    select(ID = FINREGISTRYID, primary_ICD = CODE1, ICD_version = ICDVER, Event_age = EVENT_AGE, age_first_purch) %>% 
    filter(ICD_version != "8")
  
  score_data <- ICCI::calc_cci(long2, exp_end = long2$age_first_purch)
  
  # Some indivuduals won't have records in the longitudinal files, hence no comorbidities and CCI == 0
  d <- left_join(adh, score_data, by = c("HETU" = "ID")) %>% 
    mutate(CCI_score = ifelse(is.na(CCI_score), 0, CCI_score)) %>% 
    select(HETU, CCI_score)
  
  saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.CCI.', drug,'.rds'))
  #rm(long2, score_data)
  
  # Persistence
  per <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
    select(HETU, age_first_purch)
  
  long2 <- inner_join(long, per, by = c("FINREGISTRYID" = "HETU")) %>% 
    mutate(ICDVER = ifelse(ICDVER == "O3", "10", ICDVER)) %>% 
    select(ID = FINREGISTRYID, primary_ICD = CODE1, ICD_version = ICDVER, Event_age = EVENT_AGE, age_first_purch) %>% 
    filter(ICD_version != "8")
  
  score_data <- ICCI::calc_cci(long2, exp_end = long2$age_first_purch)
  
  d <- left_join(per, score_data, by = c("HETU" = "ID")) %>% 
    mutate(CCI_score = ifelse(is.na(CCI_score), 0, CCI_score)) %>% 
    select(HETU, CCI_score)
  
  saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.CCI.', drug,'.rds'))
  rm(long2, score_data)
}

for (drug in drugs){
  define_CCI(long, drug)
}