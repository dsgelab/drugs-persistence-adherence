rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# 176min to run in parallel - 303min user

source('/data/projects/project_adherence/scripts/00_adherence_funs_finregistry.R')

cov <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv') %>% 
  select(FINREGISTRYID, date_of_birth)

ATCs <- c('^C10AA', '^C0(2|3|8|9)', '^B01AC(04|30)', '^L02B(A01|G04|G06|G03)', '^B01AF(01|02|03|07)', '^S01E(E|D)')
drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

path <- '/data/projects/project_adherence/data/purch_per_year_filt/'

compute_trajectories_adherence <- function(path, drug) {
  files <- list.files(path, pattern = paste0('*.',drug,'.index.rds'))
  print(files)
  
  purch <- rbindlist(lapply(files, function(name) readRDS(paste0(path, name))))
  
  # Extract complete trajectories
  if (drug %in% c('statins', 'blood_pressure', 'breast_cancer')) {
    t <- getTrajectoriesDates(purch)
  } else if (drug %in% c('antiplatelet', 'doac')) {
    ATCs <- data.frame(atc = c("B01AC04","B01AC30","B01AF01","B01AF02","B01AF03","B01AE07"),
                       dose = c(1,2,1,2,1,2),
                       stringsAsFactors = F)
    t <- getTrajectoriesDates2(purch,ATCs)
  } else if (drug == "glaucoma") {
    t <- getTrajectoriesDates(purch, 0.1)
  }
  saveRDS(t, file = paste0('/data/projects/project_adherence/data/trajectories/trajectories.full.',drug,'.rds'))
  
  # Summarise trajectories
  # Add event age first
  t <- left_join(t, cov, by = c("HETU" = "FINREGISTRYID")) %>% 
    mutate(EVENT_AGE = as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25)
  
  tt <- summarizeTrajectories(t, 18)
  
  saveRDS(tt, file = paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) 
}

my.cluster <- parallel::makeCluster(
  6, 
  type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% compute_trajectories_adherence(path, drugs[i])
