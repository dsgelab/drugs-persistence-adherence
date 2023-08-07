rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# real    46m13.778s
# user    57m16.193s

ep <- fread('/data/processed_data/endpointer/R8/densified_first_events_DF8_all_endpoints_2021-09-04.txt')

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

ep_secondary <- list(
    c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH', 'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS'),
    NA,
    c('I9_STR_SAH', 'I9_TIA', 'I9_MI'),
    c('C3_BREAST'),
    c('I9_DVTANDPULM','I9_VTE','I9_AF'))

define_prevention <- function(ep, drug, ep_secondary) {
    adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    if (drug != "blood_pressure") {
        # Get age of first occurence of any of the endpoints
        age_first_ep <- ep %>%
            filter(FINNGENID %in% adh$HETU,
                   ENDPOINT %in% ep_secondary) %>%
            mutate(AGE_FIRST_EP = AGE) %>%
            group_by(FINNGENID) %>%
            arrange(AGE_FIRST_EP) %>%
            filter(row_number()==1) %>%
            select(FINREGISTRYID = FINNGENID, FIRST_EP = ENDPOINT, AGE_FIRST_EP)
        
        d <- left_join(adh, age_first_ep, by = c("HETU" = "FINREGISTRYID")) %>%
            mutate(secondary = ifelse( (AGE_FIRST_EP > age_first_purch | is.na(AGE_FIRST_EP)), 0, 1)) %>% 
            select(HETU, secondary)
    
    } else {
        d <- adh %>% mutate(secondary = NA) %>% 
            select(HETU, secondary)
    }
    
    saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.prevention.', drug,'.rds'))
    rm(d)
    
    
    per <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
        select(HETU, age_first_purch)
    
    if (drug != "blood_pressure") {
        # Get age of first occurence of any of the endpoints
        age_first_ep <- ep %>%
            filter(FINNGENID %in% per$HETU,
                   ENDPOINT %in% ep_secondary) %>%
            mutate(AGE_FIRST_EP = AGE) %>%
            group_by(FINNGENID) %>%
            arrange(AGE_FIRST_EP) %>%
            filter(row_number()==1) %>%
            select(FINREGISTRYID = FINNGENID, FIRST_EP = ENDPOINT, AGE_FIRST_EP)
        
        d <- left_join(per, age_first_ep, by = c("HETU" = "FINREGISTRYID")) %>%
            mutate(secondary = ifelse( (AGE_FIRST_EP > age_first_purch | is.na(AGE_FIRST_EP)), 0, 1)) %>% 
            select(HETU, secondary)
    } else {
        d <- per %>% mutate(secondary = NA) %>% 
            select(HETU, secondary)
    }
    
    saveRDS(d, paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.prevention.', drug,'.rds'))
    rm(d)

}

my.cluster <- parallel::makeCluster(
    6, 
    type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(drugs), .verbose = TRUE) %dopar% define_prevention(ep, drugs[i], ep_secondary[[i]])