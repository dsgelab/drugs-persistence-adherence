rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# Get living information
birth_date <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv', select = c("FINREGISTRYID", "date_of_birth"))
living <- fread('/data/processed_data/dvv/dvv_living_extended/dvv_ext_core.csv', select = c("FINREGISTRYID", "Start_of_residence", "End_of_residence", "Luokka")) %>% 
    left_join(birth_date) %>% 
    mutate(start_residence_age = as.numeric(difftime(Start_of_residence, date_of_birth, units = "days"))/365.25,
           end_residence_age = as.numeric(difftime(End_of_residence, date_of_birth, units = "days"))/365.25)

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac', 'glaucoma')
drug <- 'statins'

define_living <- function(living, drug) {
    adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)

    adh_living <- left_join(adh, living, by = c("HETU" = "FINREGISTRYID"))
    
    # Individuals with complete info at baseline
    d1 <- adh_living %>% 
        filter((age_first_purch >= start_residence_age),
               (age_first_purch <= end_residence_age | is.na(end_residence_age)),
               Luokka != "") %>% 
        mutate(urban = case_when(grepl("K", Luokka) ~ 1,
                                 grepl("M", Luokka) ~ 0,
                                 TRUE ~ NA_real_)) %>% 
        group_by(HETU) %>% 
        arrange(start_residence_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    # Indivudals with Luokka info and end of residence after baseline, but NA start of residence
    d2 <- adh_living %>% 
        filter(!(HETU %in% d1$HETU),
               is.na(start_residence_age),
               age_first_purch <= end_residence_age,
               Luokka != "") %>% 
        mutate(urban = case_when(grepl("K", Luokka) ~ 1,
                                 grepl("M", Luokka) ~ 0,
                                 TRUE ~ NA_real_)) %>% 
        group_by(HETU) %>% 
        arrange(start_residence_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    # If not info at baseline, just get the first record available
    d3 <- adh_living %>%
        filter(Luokka != "",
               !(HETU %in% c(d1$HETU, d2$HETU))) %>% 
        mutate(urban = case_when(grepl("K", Luokka) ~ 1,
                                 grepl("M", Luokka) ~ 0,
                                 TRUE ~ NA_real_)) %>% 
        group_by(HETU) %>% 
        arrange(start_residence_age, .by_group = T) %>% 
        filter(row_number()==n())
        
    saveRDS(bind_rows(d1,d2,d3) %>% select(HETU, urban), paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.living.', drug,'.rds'))
    
    per <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    per_living <- left_join(per, living, by = c("HETU" = "FINREGISTRYID"))
    
    d1 <- per_living %>% 
        filter((age_first_purch >= start_residence_age),
               (age_first_purch <= end_residence_age | is.na(end_residence_age)),
               Luokka != "") %>% 
        mutate(urban = case_when(grepl("K", Luokka) ~ 1,
                                 grepl("M", Luokka) ~ 0,
                                 TRUE ~ NA_real_)) %>% 
        group_by(HETU) %>% 
        arrange(start_residence_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    # Indivudals with Luokka info and end of residence after baseline, but NA start of residence
    d2 <- per_living %>% 
        filter(!(HETU %in% d1$HETU),
               is.na(start_residence_age),
               age_first_purch <= end_residence_age,
               Luokka != "") %>% 
        mutate(urban = case_when(grepl("K", Luokka) ~ 1,
                                 grepl("M", Luokka) ~ 0,
                                 TRUE ~ NA_real_)) %>% 
        group_by(HETU) %>% 
        arrange(start_residence_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    # If not info at baseline, just get the first record available
    d3 <- per_living %>%
        filter(Luokka != "",
               !(HETU %in% c(d1$HETU,d2$HETU))) %>% 
        mutate(urban = case_when(grepl("K", Luokka) ~ 1,
                                 grepl("M", Luokka) ~ 0,
                                 TRUE ~ NA_real_)) %>% 
        group_by(HETU) %>% 
        arrange(start_residence_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    saveRDS(bind_rows(d1,d2,d3) %>% select(HETU, urban), paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.living.', drug,'.rds'))
}

my.cluster <- parallel::makeCluster(
    5, 
    type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% define_living(living, drugs[i])