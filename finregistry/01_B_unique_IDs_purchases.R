rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

source('/data/projects/project_adherence/scripts/00_adherence_funs_finregistry.R')

getTrajectoriesDiscontinuation <- function(purch,cov,dose=1,min_age) {
  # Define early stopping users
  df <- purch %>%
    #left_join(cov, by = c("HETU" = "FINREGISTRYID")) %>%
    group_by(HETU) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(PKOKO))) %>%
    # for each individual, order by event_date
    arrange(OSTOPV, .by_group = T)  %>%
    filter(n() <= 3) %>% 
    # keep trajectories with last purchase at least 2 years before end of followup (EOF = death OR emigration OR 2020-01-01)
    mutate(END_OF_FOLLOWUP = min(death_date, emigration_date, as.Date('2020-01-01'), na.rm = T),
           years_to_EOF = as.numeric(difftime(as.Date(END_OF_FOLLOWUP), as.Date(OSTOPV), units = "days"))/365.25,
           pills = PKOKO*PLKM/dose) %>%
    summarise(tot_pills = sum(pills),
              tot_purch = n(),
              last_years_to_EOF = last(years_to_EOF),
              age_first_purch = first(as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25)) %>%
    mutate(age_bin = cut(age_first_purch, breaks = c(18, 48, 60, 80, Inf), labels = c("18-39", "40-59", "60-79", "80+"), right= FALSE),
           persistent = 0) %>%
    filter(age_first_purch >= min_age,
           last_years_to_EOF > 2) %>%
    select(HETU, age_first_purch, age_bin, tot_pills, tot_purch, last_years_to_EOF, persistent)
  return(df)
}

getTrajectoriesDiscontinuation2 <- function(purch,cov,ATCs,min_age) {
  # ATCs as a df with daily dose of each ATC
  df <- purch %>%
    inner_join(ATCs, by = c("ATC"="atc")) %>%
    #left_join(cov, by = c("HETU" = "FINREGISTRYID")) %>%
    group_by(HETU) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(PKOKO))) %>%
    # for each individual, order by event_date
    arrange(OSTOPV, .by_group = T)  %>%
    filter(n() <= 3) %>% 
    # keep trajectories with last purchase at least 2 years before end of followup (EOF = death OR emigration OR 2020-01-01)
    mutate(END_OF_FOLLOWUP = min(death_date, emigration_date, as.Date('2020-01-01'), na.rm = T),
           years_to_EOF = as.numeric(difftime(as.Date(END_OF_FOLLOWUP), as.Date(OSTOPV), units = "days"))/365.25,
           pills = PKOKO*PLKM/dose) %>%
    summarise(tot_pills = sum(pills),
              tot_purch = n(),
              last_years_to_EOF = last(years_to_EOF),
              age_first_purch = first(as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25)) %>%
    mutate(age_bin = cut(age_first_purch, breaks = c(18, 48, 60, 80, Inf), labels = c("18-39", "40-59", "60-79", "80+"), right= FALSE),
           persistent = 0) %>%
    filter(age_first_purch >= min_age,
           last_years_to_EOF > 2) %>%
    select(HETU, age_first_purch, age_bin, tot_pills, tot_purch, last_years_to_EOF, persistent)
  return(df)
}

cov <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv') %>%
  select(FINREGISTRYID, date_of_birth, death_date, emigration_date)
ep <- fread('/data/processed_data/endpointer/R8/densified_first_events_DF8_all_endpoints_2021-09-04.txt')

ATCs <- c('^C10AA', '^C0(2|3|8|9)', '^B01AC(04|30)', '^L02B(A01|G04|G06|G03)', '^B01AF(01|02|03|07)')
#drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
#drugs <- c('blood_pressure')
#drugs <- c('antiplatelet')
#drug <- drugs[1]

path <- '/data/projects/project_adherence/data/purch_per_year_filt/'

count_ids_purchases <- function(path, drug) {
  files <- list.files(path, pattern = paste0('*.',drug,'.index.rds'))
  purch <- rbindlist(lapply(files, function(name) readRDS(paste0(path, name))))
  
  # keep only age >= 18 and only patients with complete info for the whole trajectory
  purch <- purch %>%
    left_join(cov, by = c("HETU" = "FINREGISTRYID")) %>% 
    mutate(EVENT_AGE = as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25) %>% 
    group_by(HETU) %>% 
    filter(all(!is.na(PKOKO)))
  
  # extract age FIRST purchase
  fp18 <- purch %>% 
    group_by(HETU) %>% 
    arrange(EVENT_AGE, .by_group = T) %>% 
    summarise(age_first_purch = first(EVENT_AGE)) %>% 
    filter(age_first_purch >= 18) %>% 
    pull(HETU)
  
  # filter purchase dataset
  purch <- purch %>% filter(HETU %in% fp18)
  
  # keep only secondary prevention for breast cancer and doac
  if (drug %in% c('breast_cancer', 'doac')) {
    
    if (drug == 'breast_cancer') {
      ep_secondary <- c('C3_BREAST')
    } else {
      ep_secondary <- c('I9_DVTANDPULM','I9_VTE','I9_AF')
    }

    age_first_ep <- ep %>%
      filter(FINNGENID %in% purch$HETU,
             ENDPOINT %in% ep_secondary) %>%
      mutate(AGE_FIRST_EP = AGE) %>%
      group_by(FINNGENID) %>%
      arrange(AGE_FIRST_EP) %>%
      filter(row_number()==1) %>%
      select(FINREGISTRYID = FINNGENID, FIRST_EP = ENDPOINT, AGE_FIRST_EP)
    
    secondary <- purch %>% 
      group_by(HETU) %>% 
      arrange(EVENT_AGE, .by_group = T) %>% 
      summarise(age_first_purch = first(EVENT_AGE)) %>% 
      left_join(age_first_ep, by = c("HETU" = "FINREGISTRYID")) %>% 
      mutate(secondary = ifelse( (AGE_FIRST_EP > age_first_purch | is.na(AGE_FIRST_EP)), 0, 1)) %>% 
      filter(secondary == 1) %>% 
      pull(HETU)
    
    purch <- purch %>% filter(HETU %in% secondary)
  }
  
  saveRDS(purch %>% select(HETU), file = paste0('/data/projects/project_adherence/data/total_IDs.',drug,'.filt_age_prev.rds'))

}

my.cluster <- parallel::makeCluster(
  length(drugs),
  type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% count_ids_purchases(path, drugs[i])

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
#drugs <- c('statins', 'antiplatelet', 'breast_cancer', 'doac')
IDs <- NULL

for (drug in drugs){
  d <- readRDS(paste0('/data/projects/project_adherence/data/total_IDs.',drug,'.filt_age_prev.rds'))
  IDs <-
    bind_rows(IDs, d %>% distinct()) %>%
    distinct()
}

IDs %>% distinct() %>% nrow()
(IDs %>% distinct() %>% nrow()) / 5339804

res <- NULL
for (drug in drugs){
  d <- readRDS(paste0('/data/projects/project_adherence/data/total_counts.',drug,'.filt_age_prev.rds'))
  res <- bind_rows(res,d)
}


# drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
# #drugs <- c('statins', 'antiplatelet', 'breast_cancer', 'doac')
# IDs <- NULL
# 
# for (drug in drugs){
#   d <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.max3purch.filt_age_prev.rds'))
#   
#   print(drug)
#   print(table(d$tot_purch))
# 
# }