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
drug <- 'statins'

run_models <- function(drug) {
  print(drug)
  # Load all predictors for adherence and binarize continuous ones
  print("ADHERENCE")
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.adherence.',drug,'.rds')) %>% 
    mutate(age_first_purch_60 = factor(ifelse(age_first_purch >= 60, 1, 0)),
           sex = factor(sex),
           CCI_score_severe_5 = factor(ifelse(CCI_score >= 5, 1, 0)),
           mother_tongue_fin_sv = factor(mother_tongue_fin_sv),
           urban = factor(urban),
           soc_prev_year = factor(soc_prev_year),
           edu_years_bachelor_16 = factor(ifelse(edu_years >= 16, 1, 0))
           )
  
  if (drug != 'blood_pressure'){
    d <- d %>% mutate(secondary = factor(secondary))
  }
  
  if (drug %in% c('blood_pressure', 'doac')){
    f <- 'adherence ~ age_first_purch_60 + yob + sex + CCI_score_severe_5 + mother_tongue_fin_sv + edu_years_bachelor_16 + urban + soc_prev_year'
  } else if (drug == 'breast_cancer') {
    f <- 'adherence ~ age_first_purch_60 + yob +  CCI_score_severe_5 + mother_tongue_fin_sv + edu_years_bachelor_16 + urban + soc_prev_year'
  } else {
    f <- 'adherence ~ age_first_purch_60 + yob + sex + secondary + CCI_score_severe_5 + mother_tongue_fin_sv + edu_years_bachelor_16 + urban + soc_prev_year'
  }
 
  model <- lm(formula = f, data = d)
  s <- tidy(model) %>% 
    mutate(drug = drug,
           adj.r.squared = rep(summary(model)$adj.r.squared, nrow(tidy(model))))
  fwrite(s, paste0('/data/projects/project_adherence/results/manuscript/predictors.binary.model.adherence.summary.',drug,'.tsv'), sep = '\t', quote = F, na = NA)
  
  # # # PERSISTENCE
  # Load all predictors for persistence and binarize continuous ones
  print("PERSISTENCE")
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.persistence.',drug,'.rds')) %>% 
    mutate(age_first_purch_60 = factor(ifelse(age_first_purch >= 60, 1, 0)),
           sex = factor(sex),
           CCI_score_severe_5 = factor(ifelse(CCI_score >= 5, 1, 0)),
           mother_tongue_fin_sv = factor(mother_tongue_fin_sv),
           urban = factor(urban),
           soc_prev_year = factor(soc_prev_year),
           edu_years_bachelor_16 = factor(ifelse(edu_years >= 16, 1, 0)),
           persistent = factor(persistent)
    )
  
  
  if (drug != 'blood_pressure'){
    d <- d %>% mutate(secondary = factor(secondary))
  }

  if (drug %in% c('blood_pressure', 'doac')){
    f <- 'persistent ~ age_first_purch_60 + yob + sex + CCI_score_severe_5 + mother_tongue_fin_sv + edu_years_bachelor_16 + urban + soc_prev_year'
  } else if (drug == 'breast_cancer') {
    f <- 'persistent ~ age_first_purch_60 + yob + CCI_score_severe_5 + mother_tongue_fin_sv + edu_years_bachelor_16 + urban + soc_prev_year'
  } else {
    f <- 'persistent ~ age_first_purch_60 + yob + sex + secondary + CCI_score_severe_5 + mother_tongue_fin_sv + edu_years_bachelor_16 + urban + soc_prev_year'
  }
  
  nullmod <- glm(persistent ~ 1, d, family="binomial")
  model <- glm(formula = f, data = d, family = "binomial")
  s <- tidy(model) %>% 
    mutate(drug = drug,
           mcfadden.r.squared = 1-logLik(model)/logLik(nullmod))
  fwrite(s, paste0('/data/projects/project_adherence/results/manuscript/predictors.binary.model.persistence.summary.',drug,'.tsv'), sep = '\t', quote = F, na = NA)
}

my.cluster <- parallel::makeCluster(
  5, 
  type = "FORK")

registerDoParallel(cl = my.cluster)
foreach(i=1:length(drugs), .verbose = TRUE) %dopar% run_models(drugs[i])


tot <- NULL
for (drug in drugs){
  d <- fread(paste0('/data/projects/project_adherence/results/manuscript/predictors.binary.model.adherence.summary.',drug,'.tsv'))
  tot <- bind_rows(tot, d)
}
fwrite(tot, '/data/projects/project_adherence/results/manuscript/predictors.binary.model.adherence.summary.tsv', sep = '\t', quote = F, na = NA)


tot <- NULL
for (drug in drugs){
  d <- fread(paste0('/data/projects/project_adherence/results/manuscript/predictors.binary.model.persistence.summary.',drug,'.tsv'))
  tot <- bind_rows(tot, d)
}
fwrite(tot, '/data/projects/project_adherence/results/manuscript/predictors.binary.model.persistence.summary.tsv', sep = '\t', quote = F, na = NA)

df <- fread('/data/projects/project_adherence/results/manuscript/predictors.binary.model.adherence.summary.tsv')
df2 <- fread('/data/projects/project_adherence/results/manuscript/predictors.binary.model.persistence.summary.tsv')
