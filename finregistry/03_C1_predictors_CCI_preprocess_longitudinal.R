rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# Filter detailed longitudinal file to keep only drug users considered in the analysis and IDs and ICD codes
drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac', 'glaucoma')

# Pull IDs to extract from longitudinal file
IDs <- c()
for (drug in drugs){
  id_ad <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',drug,'.rds')) %>% 
    pull(HETU)
  id_pe <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.persistence.',drug,'.rds')) %>% 
    pull(HETU)
  
  IDs <- unique(c(IDs,id_ad,id_pe))
}

path <- "/data/processed_data/endpointer/supporting_files/main/longitudinal.txt.*"
files <- system(paste0("ls ", path), intern = T)

filter_long <- function(file, IDs){
  long <- fread(file, drop = c("PVM", "CODE2", "CODE3", "CODE4", "CATEGORY", "INDEX"))
  long <- long %>% 
    filter(FINREGISTRYID %in% IDs,
           SOURCE %in% c("INPAT", "OUTPAT", "OPER_IN", "OPER_OUT", "PRIM_OUT", "CANC"),
           EVENT_YRMNTH >= "1987-01")
  
  saveRDS(long, file = paste0('/data/projects/project_adherence/data/longitudinal/all.drugs.users.',basename(file),'.rds'))
  gc()
}

my.cluster <- parallel::makeCluster(
  length(files), 
  type = "FORK")

registerDoParallel(cl = my.cluster)

foreach(i=1:length(files), .verbose = TRUE) %dopar% filter_long(files[i], IDs)