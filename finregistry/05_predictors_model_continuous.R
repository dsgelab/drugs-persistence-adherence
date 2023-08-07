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
drug <- 'breast_cancer'

run_models <- function(drug) {
  print(drug)
  # Load all predictors for adherence and binarize continuous ones
  print("ADHERENCE")
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.adherence.',drug,'.rds')) %>% 
    mutate(sex = factor(sex),
           mother_tongue_fin_sv = factor(mother_tongue_fin_sv),
           urban = factor(urban),
           soc_prev_year = factor(soc_prev_year)
           )
  
  if (drug != 'blood_pressure'){
    d <- d %>% mutate(secondary = factor(secondary))
  }
  
  if (drug == 'breast_cancer') {
    base <- 'adherence ~ age_first_purch + yob'
  } else {
    base <- 'adherence ~ age_first_purch + yob + sex'
  } 
  
  model <- lm(formula = base, data = d)
  base_adj_r2 <- summary(model)$adj.r.squared
  res <- NULL

  for (v in c("secondary", "CCI_score", "mother_tongue_fin_sv", "edu_years", "urban", "soc_prev_year")){
    
    if ( (drug %in% c('blood_pressure', 'doac', 'breast_cancer')) & (v == "secondary")){
      f <- base
    } else {
      f <- paste0(base, " + ", v)
    }
    
    model <- lm(formula = f, data = d)
    res_v <- data.frame(drug = drug, predictor = v, base_adj_r2 = base_adj_r2, adj_r2 = summary(model)$adj.r.squared)
    res <- bind_rows(res, res_v)
    
  }

  fwrite(res, paste0('/data/projects/project_adherence/results/manuscript/predictors.continuous.model.adherence.r2.',drug,'.tsv'), sep = '\t', quote = F, na = NA)
  
  # # # PERSISTENCE
  # Load all predictors for persistence and binarize continuous ones
  print("PERSISTENCE")
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.complete_cases.persistence.',drug,'.rds')) %>% 
    mutate(sex = factor(sex),
           mother_tongue_fin_sv = factor(mother_tongue_fin_sv),
           urban = factor(urban),
           soc_prev_year = factor(soc_prev_year),
           persistent = factor(persistent)
    )
  
  
  if (drug != 'blood_pressure'){
    d <- d %>% mutate(secondary = factor(secondary))
  }
  
  nullmod <- glm(persistent ~ 1, d, family="binomial")
  
  if (drug == 'breast_cancer') {
    base <- 'persistent ~ age_first_purch + yob'
  } else {
    base <- 'persistent ~ age_first_purch + yob + sex'
  } 
  
  model <- glm(formula = base, data = d, family = "binomial")
  nullmod <- glm(persistent ~ 1, d, family="binomial")
  
  if (drug == 'breast_cancer') {
    base_mcfadden_r2 <- 1 - (logLik(model) - 2) /logLik(nullmod) 
  } else {
    base_mcfadden_r2 <- 1 - (logLik(model) - 3) /logLik(nullmod) 
  } 

  res <- NULL
  
  for (v in c("secondary", "CCI_score", "mother_tongue_fin_sv", "edu_years", "urban", "soc_prev_year")){
    
    if ( (drug %in% c('blood_pressure', 'doac', 'breast_cancer')) & (v == "secondary")){
      f <- base
    } else {
      f <- paste0(base, " + ", v)
    }
    
    model <- glm(formula = f, data = d, family = "binomial")
    res_v <- data.frame(drug = drug, predictor = v, base_mcfadden_r2 = base_mcfadden_r2, mcfadden_r2 = 1-logLik(model)/logLik(nullmod))
    res <- bind_rows(res, res_v)
    
  }
  
  fwrite(res, paste0('/data/projects/project_adherence/results/manuscript/predictors.continuous.model.persistence.r2.',drug,'.tsv'), sep = '\t', quote = F, na = NA)
}

my.cluster <- parallel::makeCluster(
  5, 
  type = "FORK")

registerDoParallel(cl = my.cluster)
foreach(i=1:length(drugs), .verbose = TRUE) %dopar% run_models(drugs[i])


tot <- NULL
for (drug in drugs){
  d <- fread(paste0('/data/projects/project_adherence/results/manuscript/predictors.continuous.model.adherence.r2.',drug,'.tsv'))
  tot <- bind_rows(tot, d)
}
tot <- tot %>% mutate(perc_change = ((adj_r2 - base_adj_r2) / base_adj_r2)*100)
fwrite(tot, '/data/projects/project_adherence/results/manuscript/predictors.continuous.model.adherence.r2.all.tsv', sep = '\t', quote = F, na = NA)


tot <- NULL
for (drug in drugs){
  d <- fread(paste0('/data/projects/project_adherence/results/manuscript/predictors.continuous.model.persistence.r2.',drug,'.tsv'))
  tot <- bind_rows(tot, d)
}
tot <- tot %>% mutate(perc_change = ((mcfadden_r2 - base_mcfadden_r2) / base_mcfadden_r2)*100)
fwrite(tot, '/data/projects/project_adherence/results/manuscript/predictors.continuous.model.persistence.r2.all.tsv', sep = '\t', quote = F, na = NA)

df <- fread('/data/projects/project_adherence/results/manuscript/predictors.continuous.model.adherence.r2.all.tsv')
df2 <- fread('/data/projects/project_adherence/results/manuscript/predictors.continuous.model.persistence.r2.all.tsv')
