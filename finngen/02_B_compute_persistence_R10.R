rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

source('/home/ivm/drugs/scripts/FG_00_adherence_funs_finngen.R')

cov <- fread('/finngen/library-red/finngen_R10/analysis_covariates/R10_COV_V1.FID.txt.gz')
purch <- fread('/home/ivm/drugs/data/finngen_R10_purch_vnr_98_alive_2010.gz')
ep <- fread('/finngen/library-red/finngen_R10/phenotype_1.0/data/finngen_R10_endpoint_1.0.txt.gz')

ATCs <- c('^C10AA', '^C0(2|3|8|9)', '^B01AC(04|30)', '^L02B(A01|G04|G06|G03)', '^B01AF(01|02|03|07)')
drugs <- c('statins', 'blood_pressure', 'clopi_dipy', 'breast_cancer', 'doac')
ep_secondary <- list(
  c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH', 'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS'),
  NA,
  c('I9_STR_SAH', 'I9_TIA', 'I9_MI'),
  c('C3_BREAST'),
  c('I9_DVTANDPULM','I9_VTE','I9_AF'))

for (i in 1:length(drugs)) {
  drug <- drugs[i]
  print(drug)
  adh <- fread(paste0('/home/ivm/drugs/data/R10_',drug,'_summarized.txt'))
  
  # Extract persistent, discontinuation and combine them
  if (drug %in% c('statins', 'breast_cancer')) {
    p <- getTrajectoriesPersistence(adh)
    d <- getTrajectoriesDiscontinuation(purch, ATCs[i], cov, 1, 18)
    t <- bind_rows(p,d)
  
    # Age first related event
    age_first <- getAgeFirstEndpoint(ep, t$FINNGENID, ep_secondary[[i]])
    t <- t %>%
      left_join(age_first) %>%
      mutate(secondary = case_when(age_first_ev <= age_first_purch ~ 1,
                                   TRUE ~ 0))

    # Filter out non-secondary users for breast cancer
    if (drug == 'breast_cancer'){
      t <- t %>% filter(secondary == 1)
    } 
    
  } else if (drug == 'blood_pressure') {
    p <- getTrajectoriesPersistence(adh)
    d <- getTrajectoriesDiscontinuation(purch, ATCs[i], cov, 1, 18)
    t <- bind_rows(p,d)
    t$secondary <- NA
  
  } else if (drug == 'doac') {
    ATCs_df <- data.frame(atc = c("B01AF01","B01AF02","B01AF03","B01AE07"),
                       dose = c(1,2,1,2),
                       stringsAsFactors = F)
    p <- getTrajectoriesPersistence(adh)
    d <- getTrajectoriesDiscontinuation2(purch, ATCs_df, cov, 18)
    t <- bind_rows(p,d)
    
    # Age first related event
    age_first <- getAgeFirstEndpoint(ep, t$FINNGENID, ep_secondary[[i]])
    t <- t %>%
      left_join(age_first) %>%
      mutate(secondary = case_when(age_first_ev <= age_first_purch ~ 1,
                                   TRUE ~ 0)) %>%
      # Filter out non-secondary users
      filter(secondary == 1)
  
  } else if (drug == 'clopi_dipy') {
    ATCs_df <- data.frame(atc = c("B01AC04","B01AC30"),
                       dose = c(1,2),
                       stringsAsFactors = F)
    p <- getTrajectoriesPersistence(adh)
    d <- getTrajectoriesDiscontinuation2(purch, ATCs_df, cov, 18)
    t <- bind_rows(p,d)
    
    # Age first related event
    age_first <- getAgeFirstEndpoint(ep, t$FINNGENID, ep_secondary[[i]])
    t <- t %>%
      left_join(age_first) %>%
      mutate(secondary = case_when(age_first_ev <= age_first_purch ~ 1,
                                   TRUE ~ 0))
  }

  fwrite(t, paste0('/home/ivm/drugs/data/R10_',drug,'_persistence.txt'), sep = '\t', quote = F)
  
}