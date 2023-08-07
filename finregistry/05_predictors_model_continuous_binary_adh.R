rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(foreach)
library(doParallel)


# Load files with the predictors we want to use for the model

# Minimal phenotype file:
# - sex
# - year of birth
# - mother tongue
drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

res <- NULL

for (drug in drugs) {
  # Load all predictors for adherence
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.adherence.',drug,'.rds')) %>% 
    mutate(sex = factor(sex),
           mother_tongue_fin_sv = factor(mother_tongue_fin_sv),
           urban = factor(urban),
           soc_prev_year = factor(soc_prev_year),
           good_adherence = factor(ifelse(adherence > 0.8, 1, 0)))
  
  if (drug != 'blood_pressure'){
    d <- d %>% mutate(secondary = factor(secondary))
  }
  
  if (drug %in% c('blood_pressure', 'doac')){
    f <- 'good_adherence ~ age_first_purch + yob + sex + CCI_score + mother_tongue_fin_sv + edu_years + urban + soc_prev_year'
    k <- 8
  } else if (drug == 'breast_cancer') {
    f <- 'good_adherence ~ age_first_purch + yob +  CCI_score + mother_tongue_fin_sv + edu_years + urban + soc_prev_year'
    k <- 7
  } else {
    f <- 'good_adherence ~ age_first_purch + yob + sex + secondary + CCI_score + mother_tongue_fin_sv + edu_years + urban + soc_prev_year'
    k <- 9
  }
 
  # mcfadden r2 for the model
  model_intercept <- glm(good_adherence ~ 1, d, family="binomial")
  model <- glm(formula = f, data = d, family="binomial")
  adj_mcfadden <- as.numeric(1 - (logLik(model) - k)/logLik(model_intercept))
  
  # mcfadden r2 for the baseline model
  if (drug %in% c('blood_pressure', 'doac')){
    f <- 'good_adherence ~ age_first_purch + yob + sex'
    k <- 3
  } else if (drug == 'breast_cancer') {
    f <- 'good_adherence ~ age_first_purch + yob'
    k <- 2
  } else {
    f <- 'good_adherence ~ age_first_purch + yob + sex'
    k <- 3
  }
  
  model_base <- glm(formula = f, data = d, family="binomial")
  adj_mcfadden_base <- as.numeric(1 - (logLik(model_base) - k)/logLik(model_intercept))
  
  r <- data.frame(Drug = drug,
                  Pheno = "Good adherence",
                  Base_mcfadden_r2 = adj_mcfadden_base,
                  Model_mcfadden_r2 = adj_mcfadden)
  
  res <- bind_rows(res, r)

  # Load all predictors for persistence
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.persistence.',drug,'.rds')) %>% 
    mutate(sex = factor(sex),
           mother_tongue_fin_sv = factor(mother_tongue_fin_sv),
           urban = factor(urban),
           soc_prev_year = factor(soc_prev_year),
           persistent = factor(persistent))
  
  if (drug != 'blood_pressure'){
    d <- d %>% mutate(secondary = factor(secondary))
  }

  if (drug %in% c('blood_pressure', 'doac')){
    f <- 'persistent ~ age_first_purch + yob + sex + CCI_score + mother_tongue_fin_sv + edu_years + urban + soc_prev_year'
    k <- 8
  } else if (drug == 'breast_cancer') {
    f <- 'persistent ~ age_first_purch + yob + CCI_score + mother_tongue_fin_sv + edu_years + urban + soc_prev_year'
    k <- 7
  } else {
    f <- 'persistent ~ age_first_purch + yob + sex + secondary + CCI_score + mother_tongue_fin_sv + edu_years + urban + soc_prev_year'
    k <- 9
  }
  
  # mcfadden r2 for the model
  model_intercept <- glm(persistent ~ 1, d, family="binomial")
  model <- glm(formula = f, data = d, family="binomial")
  adj_mcfadden <- as.numeric(1 - (logLik(model) - k)/logLik(model_intercept))
  
  # mcfadden r2 for the baseline model
  if (drug %in% c('blood_pressure', 'doac')){
    f <- 'persistent ~ age_first_purch + yob + sex'
    k <- 3
  } else if (drug == 'breast_cancer') {
    f <- 'persistent ~ age_first_purch + yob'
    k <- 2
  } else {
    f <- 'persistent ~ age_first_purch + yob + sex'
    k <- 3
  }
  
  model_base <- glm(formula = f, data = d, family="binomial")
  adj_mcfadden_base <- as.numeric(1 - (logLik(model_base) - k)/logLik(model_intercept))
  
  r <- data.frame(Drug = drug,
                  Pheno = "Persistence",
                  Base_mcfadden_r2 = adj_mcfadden_base,
                  Model_mcfadden_r2 = adj_mcfadden)
  
  res <- bind_rows(res, r)
}

fwrite(res, paste0('/data/projects/project_adherence/results/manuscript/predictors.model.adherence.mcfadden.all.tsv'), sep = '\t', quote = F, na = NA)