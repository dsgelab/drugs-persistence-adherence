rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# real    93m36.642s
# user    135m39.318s

source('/data/projects/project_adherence/scripts/00_adherence_funs_finregistry.R')

cov <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv') %>% 
  select(FINREGISTRYID, date_of_birth, death_date, emigration_date)

ATCs <- c('^C10AA', '^C0(2|3|8|9)', '^B01AC(04|30)', '^L02B(A01|G04|G06|G03)', '^B01AF(01|02|03|07)')
drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

path <- '/data/projects/project_adherence/data/purch_per_year_filt/'

compute_trajectories_persistence <- function(path, drug) {
  files <- list.files(path, pattern = paste0('*.',drug,'.index.rds'))
  purch <- rbindlist(lapply(files, function(name) readRDS(paste0(path, name))))
  adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds'))
  
  # Extract persistent, discontinuation and combine them
  if (drug %in% c('statins', 'blood_pressure', 'breast_cancer')) {
    p <- getTrajectoriesPersistence(adh)
    d <- getTrajectoriesDiscontinuation(purch, cov, 1, 18)
    t <- bind_rows(p,d)
  } else if (drug %in% c('antiplatelet', 'doac')) {
    ATCs <- data.frame(atc = c("B01AC04","B01AC30","B01AF01","B01AF02","B01AF03","B01AE07"),
                       dose = c(1,2,1,2,1,2),
                       stringsAsFactors = F)
    p <- getTrajectoriesPersistence(adh)
    d <- getTrajectoriesDiscontinuation2(purch, cov, ATCs, 18)
    t <- bind_rows(p,d)
  }
  saveRDS(t, file = paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds'))
}

my.cluster <- parallel::makeCluster(
  length(drugs), 
  type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% compute_trajectories_persistence(path, drugs[i])