rm(list=ls())

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

# run time: 64min with parallel, 200min user

vnr_map <- fread('/data/projects/project_adherence/data/vnr_info_combined.csv')

path <- "/data/processed_data/kela_purchase/"
files <- system(paste0("ls ", path), intern = T)

# Filter out files before 1998
files <- files[!grepl('1995|1996|1997', files)]

index_person <- fread('/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv') %>% 
  filter(index_person == 1) %>% 
  pull(FINREGISTRYID)

# Specify ATC codes fro each medication class
ATCs <- c('^C10AA', '^C0(2|3|8|9)', '^B01AC(04|30)', '^L02B(A01|G04|G06|G03)', '^B01AF(01|02|03|07)')
drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')          

extractATC <- function(path, file, to_keep, ATCs, drug) {
  print(file)
  purch <- fread(paste0(path, file))
  
  purch <- purch %>% 
    select(HETU, ATC, OSTOPV, PLKM, VNRO, ANJA) %>% 
    filter(HETU %in% to_keep,
           grepl(ATCs, ATC)) %>%
    mutate(VNRO = as.integer(VNRO)) %>% 
    left_join(vnr_map)
  
  saveRDS(purch, file = paste0('/data/projects/project_adherence/data/purch_per_year_filt/',file,'.',drug,'.index.rds'))
  gc()
}

my.cluster <- parallel::makeCluster(
  length(files), 
  type = "FORK")

registerDoParallel(cl = my.cluster)

for (j in 1:length(drugs)) {
  print(drugs[j])
  foreach(i=1:length(files)) %dopar% extractATC(path, files[i], index_person, ATCs[j], drugs[j])
}