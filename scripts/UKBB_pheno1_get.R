#### SET ####
rm(list=ls())

library(ggplot2)
library(data.table)
library(lubridate)
library(readxl)
library(dplyr)
library(Hmisc)
require(cowplot)
require(dplyr)
require(readr)
require(grid)
require(gridExtra)
library(ggpubr)
library(viridis)
library(lubridate)

source("~/code/R_rainclouds.R")
source('~/code/adherence_funs.R')

f <- function(x){
  list(mean(x),median(x),sd(x))
}
get_stat <- function(df){
  ds <- df %>% select(tot_pills,tot_days,tot_purch,age_first_purch,adherence)
  stat <- sapply(ds,f)
  rownames(stat) <- c("mean","median","std_dev")
  return(stat)
}

#covariates 
cov <- fread('~/data/export_for_map_analysis.tsv') %>% 
  select(eid = f.eid, sex = f.31.0.0, date_recruitment = f.53.0.0, age_recruitment = f.21022.0.0)
cov$eid <- as.character(cov$eid)

#pca covariates
pops <- fread("~/data/all_pops_non_eur_pruned_within_pop_pc_covs.tsv")

#samples in bgen
sam <- fread("~/genotype/bgen.samples.txt")

####  DOAC ####

doac_sum <- fread("~/data/doac_summary.tsv") %>%
  #filter for pca we have
  filter(eid %in% pops$s) %>%
  #filter for imputed bgen we have
  filter(eid %in% sam$V1) %>%
  #filter for overadherers
  filter(!adherence >= 1.1)  %>% 
  #create good adherers
  mutate(eid=as.character(eid),
         isgood = if_else(adherence >= 0.8, 1, 0)) %>%
  #merge with covariates info
  left_join(cov[,c("eid","sex")], by = "eid") 



# keep only individuls and purchases included in adherence calculation
doac_purch <- fread("~/data/doac_traj.tsv") %>%
  filter(eid %in% sam$V1) %>%
  filter(eid %in% pops$s) %>%
  filter(eid %in% doac_sum$eid,
         !is.na(days_next_purch))

# write pheno
fwrite(doac_sum,"~/data/doac_pheno.tsv", sep = "\t", quote = F)
fwrite(doac_purch,"~/data/doac_presc.tsv", sep = "\t", quote = F)

#### BLOODPRESSURE ####

bp_sum <- fread("~/data/bp_summary.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% sam$V1,
         !adherence >= 1.1) %>%
  mutate(eid=as.character(eid),
         isgood = if_else(adherence >= 0.8, 1, 0)) %>%
  left_join(cov[,c("eid","sex")], by = "eid") 

bp_purch <- fread("~/data/bp_traj.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% bp_sum$eid,
         !is.na(days_next_purch))



# write pheno
fwrite(bp_sum,"~/data/bp_pheno.tsv", sep = "\t", quote = F)
fwrite(bp_purch,"~/data/bp_presc.tsv", sep = "\t", quote = F)

#### BREST CANCER ####

bc_sum <- fread("~/data/bc_summary.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% sam$V1,
         !adherence >= 1.1 ) %>%
  mutate(eid=as.character(eid),
         isgood = if_else(adherence >= 0.8, 1, 0)) %>%
  left_join(cov[,c("eid","sex")], by = "eid") 

bc_purch <- fread("~/data/bc_traj.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% bc_sum$eid,
         !is.na(days_next_purch)) 



# write pheno
fwrite(bc_sum,"~/data/bc_pheno.tsv", sep = "\t", quote = F)
fwrite(bc_purch,"~/data/bc_presc.tsv", sep = "\t", quote = F)

#### ANTIPLATELET ####

ap_sum <- fread("~/data/ap_summary.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% sam$V1,
         !adherence >= 1.1 ) %>%
  mutate(eid=as.character(eid),
         isgood = if_else(adherence >= 0.8, 1, 0)) %>%
  left_join(cov[,c("eid","sex")], by = "eid") 

ap_purch <- fread("~/data/ap_traj.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% ap_sum$eid,
         !is.na(days_next_purch))

# write pheno
fwrite(ap_sum,"~/data/ap_pheno.tsv", sep = "\t", quote = F)
fwrite(ap_purch,"~/data/ap_presc.tsv", sep = "\t", quote = F)


#### GLAUCOMA ####

gla_sum <- fread("~/data/gla_summary.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% sam$V1,
         !adherence >= 1.1 ) %>%
  mutate(eid=as.character(eid),
         isgood = if_else(adherence >= 0.8, 1, 0)) %>%
  left_join(cov[,c("eid","sex")], by = "eid") 

gla_purch <- fread("~/data/gla_traj.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% gla_sum$eid,
         !is.na(days_next_purch))





# write pheno
fwrite(gla_sum,"~/data/gla_pheno.tsv", sep = "\t", quote = F)
fwrite(gla_purch,"~/data/gla_presc.tsv", sep = "\t", quote = F)


#### STATINS ####

stat_sum <- fread("~/data/st_summary.tsv") %>%
  filter(!adherence >= 1.1,
         eid %in% sam$V1,
         eid %in% pops$s) %>%
  mutate(eid=as.character(eid),
         isgood = if_else(adherence >= 0.8, 1, 0)) %>%
  left_join(cov[,c("eid","sex")], by = "eid") 


stat_purch <- fread("~/data/st_traj.tsv") %>%
  filter(eid %in% pops$s,
         eid %in% stat_sum$eid,
         !is.na(days_next_purch))



# write pheno
fwrite(stat_sum,"~/data/st_pheno.tsv", sep = "\t", quote = F)
fwrite(stat_purch,"~/data/st_presc.tsv", sep = "\t", quote = F)



