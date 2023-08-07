rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# Add years of education for each education level and calculate age
edu <- fread('/data/processed_data/sf_socioeconomic/tutkinto_u1442_a.csv.finreg_IDsp', select = c("FINREGISTRYID", "vuosi", "kaste_t2")) %>% 
    # NAs represent primary or lower secondary education, which corresponds to either 7 or 10 years of education. We'll set NAs to 8.5
    mutate(edu_level = case_when(is.na(kaste_t2) ~ 1,
                                TRUE ~ as.numeric(substr(kaste_t2, 1, 1))))

birth_year <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv', select = c("FINREGISTRYID", "date_of_birth")) %>% 
    mutate(year_of_birth = as.integer(format(date_of_birth, "%Y"))) %>% 
    select(FINREGISTRYID, year_of_birth)

edu_years <- data.frame(edu_level = c(1,3,4,5,6,7,8),
                        edu_years = c(8.5, 13, 15, 16, 16, 18, 22))

# calculate age and set age to 13 for the primary education records
edu <- left_join(edu, birth_year) %>% 
    mutate(edu_age = ifelse(is.na(vuosi), 13, vuosi - year_of_birth))

# add years of education
edu <- left_join(edu, edu_years) %>% 
    mutate(edu_years_scaled = as.numeric(scale(edu_years))) %>% 
    select(FINREGISTRYID, edu_age, edu_level, edu_years, edu_years_scaled)

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
# drug <- 'statins'

define_years_edu <- function(edu, drug) {
    adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    d1 <- left_join(adh, edu, by = c("HETU" = "FINREGISTRYID")) %>% 
        filter(edu_age <= age_first_purch) %>% 
        group_by(HETU) %>% 
        arrange(edu_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    # if not at baseline, just take the first available record
    d2 <- left_join(adh, edu, by = c("HETU" = "FINREGISTRYID")) %>% 
        filter(!(HETU %in% d1$HETU)) %>% 
        group_by(HETU) %>% 
        arrange(edu_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    d <- bind_rows(d1,d2) %>% 
        select(-edu_age, -age_first_purch)
    
    saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.edu.', drug,'.rds'))
    
    rm(d, d1, d2)
    
    per <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    d1 <- left_join(per, edu, by = c("HETU" = "FINREGISTRYID")) %>% 
        filter(edu_age <= age_first_purch) %>% 
        group_by(HETU) %>% 
        arrange(edu_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    d2 <- left_join(per, edu, by = c("HETU" = "FINREGISTRYID")) %>% 
        filter(!(HETU %in% d1$HETU)) %>% 
        group_by(HETU) %>% 
        arrange(edu_age, .by_group = T) %>% 
        filter(row_number()==n())
    
    d <- bind_rows(d1,d2)%>% 
        select(-edu_age, -age_first_purch)
    
    saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.edu.', drug,'.rds'))
    rm(d, d1, d2)
}

my.cluster <- parallel::makeCluster(
    5, 
    type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% define_years_edu(edu, drugs[i])