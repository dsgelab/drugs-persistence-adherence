rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)

# Filter breast_cancer and doac adherence and persistence files to keep only 2nd prevention individuals

for (drug in c("breast_cancer", "doac")) {
  
  id_adh_sec <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.adherence.prevention.',drug,'.rds')) %>% 
    filter(secondary == 1) %>% 
    pull(HETU)
  
  id_per_sec <- readRDS(paste0('/data/projects/project_adherence/data/predictors/predictors.persistence.prevention.',drug,'.rds')) %>% 
    filter(secondary == 1) %>% 
    pull(HETU)
  
  id_to_keep <- unique(c(id_adh_sec, id_per_sec))
  
  adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds'))
  saveRDS(adh, paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.ori.rds'))
  adh %>% 
    filter(HETU %in% id_to_keep) %>% 
    saveRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds'))
  
  per <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds'))
  saveRDS(per, paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'ori.rds'))
  per %>% 
    filter(HETU %in% id_to_keep) %>% 
    saveRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds'))
}