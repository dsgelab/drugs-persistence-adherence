rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

source('/data/projects/project_adherence/scripts/00_adherence_funs_finregistry.R')

cov <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv') %>% 
  select(FINREGISTRYID, date_of_birth)

ATCs <- c('^C10AA', '^C0(2|3|8|9)', '^B01AC(04|30)', '^L02B(A01|G04|G06|G03)', '^B01AF(01|02|03|07)', '^S01E(E|D)')
drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

drug <- 'breast_cancer'

tot_N <- NULL
tot_over <- NULL

compute_n_exclusion <- function(drug, cov) {
  # Read in full trajectory
  d <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.full.',drug,'.rds'))

  # Add event age first
  d <- left_join(d, cov, by = c("HETU" = "FINREGISTRYID")) %>% 
    mutate(EVENT_AGE = as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25)
  
  # Id breast cancer or DOAC, keep only people buying after 2nd prevention
  if (drug %in% c('breast_cancer', 'doac')) {
    IDs <- unique(readRDS(paste0('/data/projects/project_adherence/data/total_IDs.',drug,'.filt_age_prev.rds')) %>% pull(HETU))
    d <- d %>% filter(HETU %in% IDs)
  }
  
  # Summarise trajectories without filtering for adherence
  t <- summarizeTrajectories(d, 18, Inf)
  saveRDS(t, file = paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.w_overbuying.rds')) 
}

my.cluster <- parallel::makeCluster(
  5, 
  type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% compute_n_exclusion (drugs[i], cov)

res <- NULL
unique_eligible <- c()
unique_adh <- c()

for (drug in drugs){
  d <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.w_overbuying.rds'))

  n1year <- d %>% filter(tot_days >= 365) %>% nrow()
  n1.1 <- d %>% filter(adherence <= 1.1) %>% nrow()

  id1year <- d %>% filter(tot_days >= 365) %>% pull(HETU) %>% unique()
  id1.1 <- d %>% filter(adherence <= 1.1) %>% pull(HETU) %>% unique()

  r <- data.frame(drug = drug, N_eligible = n1year, N_adh = n1.1)
  res <- bind_rows(res,r)

  unique_eligible <- unique(c(unique_eligible, id1year))
  unique_adh <- unique(c(unique_adh, id1.1))
}

length(unique_eligible)
length(unique_adh)

unique_adh <- c()
unique_per <- c()

for (drug in drugs){
  a <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.adherence.',drug,'.rds'))
  if(drug != 'antiplatelet') {
    p <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.persistence.',drug,'.rds')) %>% 
      filter(persistent == 0)
  }
  unique_adh <- unique(c(unique_adh, a$HETU))
  unique_per <- unique(c(unique_per, p$HETU))
}

length(unique_adh)
length(unique_per)
length(unique(c(unique_adh, unique_per)))

drug <- 'antiplatelet'
a <-  readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.adherence.',drug,'.rds'))
