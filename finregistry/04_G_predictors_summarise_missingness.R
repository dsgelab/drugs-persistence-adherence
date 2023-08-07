rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(foreach)
library(doParallel)

# Load files with the predictors we want to use for the model

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')
drug <- 'blood_pressure'
summarise_missingness <- function(demo, drug) {

  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.adherence.',drug,'.rds')) 
  if(drug == 'blood_pressure') {
    d <- d %>% mutate(secondary = 0)
  }
  
  prop_na <- function(x) round(sum(is.na(x))/length(x)*100, 2)
  
  na_summary <- d %>%
    summarise_all(list(na = prop_na)) %>% 
    mutate(drug = drug,
           tot_users = nrow(d),
           complete_rows = nrow(d %>% filter(complete.cases(.))),
           complete_rows_perc = round(nrow(d %>% filter(complete.cases(.)))/nrow(d)*100, 2),
           perc_in_soc_reg = round(sum(d$in_soc_reg)/nrow(d)*100, 2)) %>% 
    select(drug, tot_users, complete_rows, complete_rows_perc, everything())

  fwrite(na_summary, paste0('/data/projects/project_adherence/results/manuscript/predictors.adherence.summary.',drug,'.tsv'), sep = '\t', quote = F, na = NA)
  
  # Load all predictors for persistence
  d <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors_combined.persistence.',drug,'.rds')) 
  d <- d %>% filter(persistent == 0)
  if(drug == 'blood_pressure') {
    d <- d %>% mutate(secondary = 0)
  }
  
  prop_na <- function(x) round(sum(is.na(x))/length(x)*100, 2)
  
  na_summary <- d %>%
    summarise_all(list(na = prop_na)) %>% 
    mutate(drug = drug,
           tot_users = nrow(d),
           complete_rows = nrow(d %>% filter(complete.cases(.))),
           complete_rows_perc = round(nrow(d %>% filter(complete.cases(.)))/nrow(d)*100, 2),
           perc_in_soc_reg = round(sum(d$in_soc_reg)/nrow(d)*100, 2)) %>% 
    select(drug, tot_users, complete_rows, complete_rows_perc, everything())

  fwrite(na_summary, paste0('/data/projects/project_adherence/results/manuscript/predictors.persistence.summary.',drug,'.tsv'), sep = '\t', quote = F, na = NA)
}

my.cluster <- parallel::makeCluster(
  5, 
  type = "FORK")

registerDoParallel(cl = my.cluster)
foreach(i=1:length(drugs), .verbose = TRUE) %dopar% summarise_missingness(demo, drugs[i])

tot <- NULL
for (drug in drugs){
  d <- fread(paste0('/data/projects/project_adherence/results/manuscript/predictors.adherence.summary.',drug,'.tsv'))
  tot <- bind_rows(tot, d)
}
fwrite(tot, '/data/projects/project_adherence/results/manuscript/predictors.adherence.summary.tsv', sep = '\t', quote = F, na = NA)


tot <- NULL
for (drug in drugs){
  d <- fread(paste0('/data/projects/project_adherence/results/manuscript/predictors.persistence.summary.',drug,'.tsv'))
  tot <- bind_rows(tot, d)
}
fwrite(tot, '/data/projects/project_adherence/results/manuscript/predictors.persistence.summary.tsv', sep = '\t', quote = F, na = NA)


df <- fread('/data/projects/project_adherence/results/manuscript/predictors.adherence.summary.tsv') %>% 
  select(drug, tot_users, complete_rows_perc, everything())

df2 <- fread('/data/projects/project_adherence/results/manuscript/predictors.persistence.summary.tsv')%>% 
  select(drug, tot_users, complete_rows_perc, everything())

# # # # Plot
# 
# # "melt" data
# df_long_adh <- df %>% 
#   select(drug, age_first_purch_na:soc_prev_year_na) %>% 
#   pivot_longer(!drug, names_to = "predictor", values_to = "perc_na") %>% 
#   filter(perc_na != 0) %>% 
#   mutate(predictor = gsub("_na", "", predictor))
# 
# df_long_per <- df2 %>% 
#   select(drug, age_first_purch_na:soc_prev_year_na) %>% 
#   pivot_longer(!drug, names_to = "predictor", values_to = "perc_na") %>% 
#   filter(perc_na != 0) %>% 
#   mutate(predictor = gsub("_na", "", predictor))
# 
# library(ggplot2)
# options(bitmapType='cairo')
# 
# ggplot(df_long_adh, aes(x = predictor, y = perc_na, fill = drug)) +
#   geom_bar(stat='identity', position='dodge') +
#   geom_text(aes(label = perc_na, ), position = position_dodge(.9)) +
#   coord_flip() +
#   theme_minimal()
# 
# ggplot(df_long_per, aes(x = predictor, y = perc_na, fill = drug)) +
#   geom_bar(stat='identity', position='dodge') +
#   geom_text(aes(label = perc_na, ), position = position_dodge(.9)) +
#   coord_flip() +
#   theme_minimal()
