rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# Add years of education for each education level and calculate age
birth_date <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv', select = c("FINREGISTRYID", "date_of_birth"))
soc <- readRDS('/data/projects/project_adherence/data/soc_assistance_date.rds') %>% 
    left_join(birth_date) %>% 
    mutate(age = as.numeric(difftime(date, date_of_birth, units = "days"))/365.25)

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
drug <- 'doac'
define_social <- function(soc, drug) {
    # Adherence
    adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    soc_prev_year <- left_join(adh, soc, by = c("HETU" = "FINREGISTRYID")) %>% 
        filter(age >= age_first_purch-1,
               age <= age_first_purch) %>% 
        group_by(HETU) %>% 
        summarise(soc_prev_year = ifelse(sum(received) == 0, 0, 1))
    
    # Assume those are not in the register have never received social assistance, but keep track of who is in the register or not
    d <- left_join(adh, soc_prev_year) %>% 
        mutate(in_soc_reg = ifelse(is.na(soc_prev_year), 0, 1),
               soc_prev_year = ifelse(is.na(soc_prev_year), 0, soc_prev_year)) %>% 
        select(-age_first_purch)
    

    saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.social.', drug,'.rds'))
    rm(d)
    
    # Persistence
    per <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    soc_prev_year <- left_join(per, soc, by = c("HETU" = "FINREGISTRYID")) %>% 
        filter(age >= age_first_purch-1,
               age <= age_first_purch) %>% 
        group_by(HETU) %>% 
        summarise(soc_prev_year = ifelse(sum(received) == 0, 0, 1))
    
    # Assume those are not in the register have never received social assistance, but keep track of who is in the register or not
    d <- left_join(per, soc_prev_year) %>% 
        mutate(in_soc_reg = ifelse(is.na(soc_prev_year), 0, 1),
               soc_prev_year = ifelse(is.na(soc_prev_year), 0, soc_prev_year))%>% 
        select(-age_first_purch)
    
    saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.social.', drug,'.rds'))
    rm(d)
}

my.cluster <- parallel::makeCluster(
    6, 
    type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% define_social(soc, drugs[i])