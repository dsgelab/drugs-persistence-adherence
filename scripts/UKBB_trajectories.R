####SET####

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

source("~/code/R_rainclouds.R")
source('~/code/adherence_funs.R')


# Exctract individual level purchase trajectories from ukbb data
getTrajectoriesUKB <- function(scripts,dose=1) {
  df <- scripts %>%
    group_by(eid) %>%
    # remove duplicated
    distinct() %>%
    # filter trajectories with at least 4 purchases, and with the same provider for all
    filter(n()>=4,
           n_distinct(data_provider) == 1,
           issue_date >= date_recruitment) %>%
    # keep only trajectories with quantities for all the events
    filter(all(!is.na(quantity_num))) %>%
    # for each individual, order by event_age
    arrange(issue_date, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = as.numeric(difftime(lead(issue_date), issue_date, units = "d")),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific dose
      n_pills = as.numeric(as.numeric(quantity_num)/dose)) %>%
    # tODO: deal with duplicated where days_next is 0. just remove them for now
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = drug_name != lag(drug_name),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      n_pills = case_when(lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      # exclude gap >=150 days without pills
      days_next_purch = case_when(days_next_purch > n_pills+150 ~ NA_real_,
                                  TRUE ~ days_next_purch),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                          TRUE ~ n_pills),
      pills_norm = n_pills/days_next_purch,
      days_norm = days_next_purch/n_pills)
  return(df)
}

summarizeTrajectoriesUKB <- function(traj, thr=3.5) {
  t <- traj %>%
    group_by(eid) %>%
    summarise(tot_pills = sum(n_pills, na.rm = T),
              tot_days = sum(days_next_purch, na.rm = T),
              mean_days = mean(days_next_purch, na.rm = T),
              sd_days = sd(days_next_purch, na.rm = T),
              mean_days_norm = mean(pills_norm, na.rm = T),
              sd_days_norm = sd(days_norm, na.rm = T),
              mean_pills_norm = mean(pills_norm, na.rm = T),
              sd_pills_norm = sd(pills_norm, na.rm = T),
              tot_purch = length(which(!is.na(n_pills))),
              age_first_purch = first(lubridate::time_length(difftime(first(issue_date), date_recruitment), "years") + age_recruitment) ) %>%
    mutate(adherence = tot_pills/tot_days) %>%
    # filter out trajectories which result to have NA for the previous metrics (e.g. if there is one purchase left)
    filter_all(all_vars(!is.na(.))) %>%
    # filter out outliers
    filter(adherence <= thr) %>%
    mutate(age_bin = cut(age_first_purch, 5)) %>%
    mutate(adherence_std = as.numeric(scale(adherence)))
}

getTrajectoriesUKBgla <- function(scripts,dose=1) {
  df <- scripts %>%
    group_by(eid) %>%
    # remove duplicated
    distinct() %>%
    # filter trajectories with at least 4 purchases, and with the same provider for all
    filter(n()>=4,
           n_distinct(data_provider) == 1,
           issue_date >= date_recruitment,
           !type=="ambiguous") %>%
    # keep only trajectories with quantities for all the events
    filter(all(!is.na(quantity_num))) %>%
    # for each individual, order by event_age
    arrange(issue_date, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = as.numeric(difftime(lead(issue_date), issue_date, units = "d")),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific dose
      n_pills = as.numeric(if_else(type=="ml",
                                   (as.numeric(quantity_num)/0.05)/dose,
                                   (as.numeric(quantity_num)/dose) ) ) ) %>%
    # tODO: deal with duplicated where days_next is 0. just remove them for now
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = drug_name != lag(drug_name),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      n_pills = case_when(lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      # exclude gap >=150 days without pills
      days_next_purch = case_when(days_next_purch > n_pills+150 ~ NA_real_,
                                  TRUE ~ days_next_purch),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                          TRUE ~ n_pills),
      pills_norm = n_pills/days_next_purch,
      days_norm = days_next_purch/n_pills)
  return(df)
}


# functions to transform quantities
tablet <- function(x){
  return(as.numeric(gsub(".*?([0-9]+).*", "\\1", x)))
}

pack_per_tablet <- Vectorize(function(x) {
  qts <- as.integer(regmatches(x, gregexpr('[0-9]+',x))[[1]])
  return(ifelse(length(qts) == 1, qts[1], qts[1]*qts[2]))
}, vectorize.args = "x")

bottle <- function(x){
  return(as.numeric(gsub(".*?([0-9]+(\\.[0-9]+)?).*", "\\1", x)))
}

pack_per_bottle <- Vectorize(function(x) {
  qts <- as.numeric(regmatches(x, gregexpr('[0-9]+(\\.[0-9]+)?',x))[[1]])
  return(ifelse(length(qts) == 1, qts[1], exp(sum(log(qts))) ) )
},vectorize.args = "x")


#### READ ####

#remove those who withdrawn consent
w <- fread('~/data/ukb31063.withdrawn_participants.20210809.csv')

# f.31.0.0: sex
# f.53.0.0: Date of attending assessment centre
# f.21022.0.0: age at recruitment
cov <- fread('~/data/export_for_map_analysis.tsv') %>% 
  select(eid = f.eid, sex = f.31.0.0, date_recruitment = f.53.0.0, age_recruitment = f.21022.0.0)
cov$eid <- as.character(cov$eid)

#### STATINS TRAJECTORIES ####

#read and exclude withdrawn 
st <- fread("~/data/st_scripts.tsv") %>%
  filter(!eid %in% w$V1)
st$eid <- as.character(st$eid)


#Merge data
st_all <- st %>%
  left_join(cov, by = "eid") 
st_all$drug_name <- tolower(st_all$drug_name)


#get traj
st_t <- getTrajectoriesUKB(st_all)

#filter low frequency drug names 
tab_st <- sort(table(st_t$drug_name))
# nothing strange so we do not eliminate any drugs_name

#get median of pills_norm
st_d <- st_t %>%
  group_by(drug_name) %>%
  summarise_at(vars(pills_norm), funs(median(., na.rm=TRUE)))

summary(st_d)
table(round(st_d$pills_norm))

# set dose variable (fiter at maximum 2 per day and minimum 1)
st_d <- st_d %>% filter(pills_norm < 3) %>%
  mutate(dose = ifelse(pills_norm < 1, 1, round(pills_norm)))

summary(st_d)
table(st_d$dose)

# merge
st_all <- st_all %>% 
  left_join(st_d, by = "drug_name") %>%
  filter(!is.na(dose))

#get dose-weighted traj 
st_t <- getTrajectoriesUKB(st_all, dose=st_all$dose)
dim(st_t)

#summarize
st_s <- summarizeTrajectoriesUKB(st_t) 
dim(st_s)

#plot
p <- ggplot(data = st_s, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  ggtitle('Overall for STAT') +
  theme_minimal()

p

#median adherence
median(st_s$adherence)

##write tsv
fwrite(st_t,"~/data/st_traj.tsv",
       quote=F,sep = "\t")
fwrite(st_s,"~/data/st_summary.tsv",
       quote=F, sep = "\t")




#### DOAC TRAJECTORIES ####

doac <- fread("~/data/doac_scripts.tsv") %>%
  filter(!eid %in% w$V1)
doac$eid <- as.character(doac$eid)
                         
#Merge data
doac_all <- doac %>%
  left_join(cov, by = "eid") 
doac_all$drug_name <- tolower(doac_all$drug_name)

#get traj
doac_t <- getTrajectoriesUKB(doac_all)

#filter low frequency drug names 
tab_doac <- sort(table(doac_t$drug_name))
# nothing strange so we do not eliminate any drugs_name

#get median of pills_norm
doac_d <- doac_t %>%
  group_by(drug_name) %>%
  summarise_at(vars(pills_norm), funs(median(., na.rm=TRUE)))

summary(doac_d)
table(round(doac_d$pills_norm))

# set dose variable
doac_d <- doac_d %>% filter(pills_norm < 3) %>%
  mutate(dose = ifelse(pills_norm < 1, 1, round(pills_norm)))

summary(doac_d)
table(doac_d$dose)

doac_all <- doac_all %>% 
  left_join(doac_d, by = "drug_name") %>%
  filter(!is.na(dose))  


#get dose-weighted traj 
doac_t <- getTrajectoriesUKB(doac_all, dose=doac_all$dose)
dim(doac_t)

#summarize
doac_s <- summarizeTrajectoriesUKB(doac_t) 
dim(doac_s)



p1 <- ggplot(data = doac_s, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  ggtitle('Overall for DOAC') +
  theme_minimal()

p1

#median adherence
median(doac_s$adherence)

##write tsv
fwrite(doac_t,"~/data/doac_traj.tsv",
       quote=F,sep = "\t")
fwrite(doac_s,"~/data/doac_summary.tsv",
       quote=F, sep = "\t")


#### BLOODPRESSURE TRAJECTORIES ####

bp <- fread("~/data/bp_scripts.tsv") %>%
  filter(!eid %in% w$V1)
bp$eid <- as.character(bp$eid)

#Merge data
bp_all <- bp %>% 
  left_join(cov, by = "eid")

bp_all$drug_name <- tolower(bp_all$drug_name)

#get traj
bp_t <- getTrajectoriesUKB(bp_all)

#filter low frequency drug names
tab_bp <- sort(table(bp_t$drug_name))


#get median pills_norm
bp_d <- bp_t %>%
  group_by(drug_name) %>%
  summarise_at(vars(pills_norm), funs(median(., na.rm=TRUE)))

summary(bp_d)
table(round(bp_d$pills_norm))

# set dose variable
bp_d <- bp_d %>% filter(pills_norm < 3) %>%
  mutate(dose = ifelse(pills_norm < 1, 1, round(pills_norm)))

summary(bp_d)
table(bp_d$dose)

bp_all <- bp_all %>% 
  left_join(bp_d, by = "drug_name")%>%
  filter(!is.na(dose))  



bp_t <- getTrajectoriesUKB(bp_all, dose=bp_all$dose) 
dim(bp_t)

bp_s <- summarizeTrajectoriesUKB(bp_t,thr=3)
dim(bp_s)


p2 <- ggplot(data = bp_s, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2,3,4,5,6)) +
  ggtitle('Overall for Blood Pressure') +
  theme_minimal()

p2

median(bp_s$adherence)

fwrite(bp_t,"~/data/bp_traj.tsv",
       quote=F,sep = "\t")
fwrite(bp_s,"~/data/bp_summary.tsv",
       quote=F,sep = "\t")

#### ANTIPLATELET TRAJECTORIES ####

ap <- fread("~/data/ap_scripts.tsv") %>%
  filter(!eid %in% w$V1)
ap$eid <- as.character(ap$eid)

#Merge data
ap_all <- ap %>% 
  left_join(cov, by = "eid") 

ap_all$drug_name <- tolower(ap_all$drug_name)


#get traj
ap_t <- getTrajectoriesUKB(ap_all)


#filter low frequency drug names
tab_ap <- sort(table(ap_t$drug_name))


ap_d <- ap_t %>%
  group_by(drug_name) %>%
  summarise_at(vars(pills_norm), funs(median(., na.rm=TRUE)))

table(round(ap_d$pills_norm))
summary(ap_d)

# set dose variable
ap_d <- ap_d %>% filter(pills_norm < 3) %>%
  mutate(dose = ifelse(pills_norm < 1, 1, round(pills_norm)))

table(ap_d$dose)
summary(ap_d)


ap_all <- ap_all %>% 
  left_join(ap_d, by = "drug_name")%>%
  filter(!is.na(dose))  


#get dose-weighted traj
ap_t <- getTrajectoriesUKB(ap_all, dose=ap_all$dose)
dim(ap_t)

ap_s <- summarizeTrajectoriesUKB(ap_t)
dim(ap_s)


p3 <- ggplot(data = ap_s, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  ggtitle('Overall for Antiplatelet') +
  theme_minimal()

p3

median(ap_s$adherence)

fwrite(ap_t,"~/data/ap_traj.tsv",
       quote=F, sep = "\t")
fwrite(ap_s,"~/data/ap_summary.tsv",
       quote=F, sep = "\t")

#### BREAST CANCER TRAJECTORIES ####

bc <- fread("~/data/bc_scripts.tsv") %>%
  filter(!eid %in% w$V1)
bc$eid <- as.character(bc$eid)

#Merge data
bc_all <- bc %>% 
  left_join(cov, by = "eid") 

bc_all$drug_name <- tolower(bc_all$drug_name)


#get first traj
bc_t <- getTrajectoriesUKB(bc_all)


#filter low frequency drug names
tab_bc <- sort(table(bc_t$drug_name))


#get median of pills_norm
bc_d <- bc_t %>%
  group_by(drug_name) %>%
  summarise_at(vars(pills_norm), funs(median(., na.rm=TRUE)))

summary(bc_d)
table(round(bc_d$pills_norm))
# set dose variable
bc_d <- bc_d %>% filter(pills_norm < 3) %>%
  mutate(dose = ifelse(pills_norm < 1, 1, round(pills_norm)))

summary(bc_d)
table(bc_d$dose)

bc_all <- bc_all %>% 
  left_join(bc_d, by = "drug_name") %>%
  filter(!is.na(dose))  


#get dose-weighted traj 
bc_t <- getTrajectoriesUKB(bc_all, dose=bc_all$dose)
dim(bc_t)

#summarize
bc_s <- summarizeTrajectoriesUKB(bc_t)
dim(bc_s)

# distribution of adherence
p4 <- ggplot(data = bc_s, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  ggtitle('Overall for Breast Cancer') +
  theme_minimal()

p4

#median adherence
median(bc_s$adherence)

fwrite(bc_t,"~/data/bc_traj.tsv",
       quote=F,sep = "\t")
fwrite(bc_s,"~/data/bc_summary.tsv",
       quote=F,sep = "\t")

#### GLAUCOMA TRAJECTORIES ####

#read
gla <- fread("~/data/gla_scripts.tsv") %>%
  filter(!eid %in% w$V1)
gla$eid <- as.character(gla$eid)

#Merge data
gla_all <- gla %>% 
  filter(grepl("^1106",bnf_code)) %>% 
  left_join(cov, by = "eid") 

gla_all$drug_name <- tolower(gla_all$drug_name)
gla_all$quantity <- tolower(gla_all$quantity)
gla_all$quantity <- gsub(".{1,5}%","",gla_all$quantity)
gla_all$quantity <- gsub(".{1,3}micrograms","",gla_all$quantity)
gla_all$quantity <- gsub(".{1,4}mg","",gla_all$quantity)



gla_all <- gla_all %>%
  mutate(quantity_num = case_when(
    grepl('ml|millilitre|dose', quantity) &
      !grepl('pack|bottle', quantity) &
      !grepl('\\*', quantity) &
      !grepl('x', quantity)
    ~ bottle(quantity),
    grepl('ml|milliliter|dose', quantity) &
      (grepl('pack|bottle', quantity) |
       grepl('\\*', quantity) |
       grepl('x', quantity))    ~ as.numeric(pack_per_bottle(quantity)),
    TRUE ~ as.numeric(quantity)),
    type = case_when(grepl("dose",quantity) ~ "dose",
                     grepl("ml|millilitre",quantity) ~ "ml",
                     TRUE ~ "ambiguous"))


table(gla_all$type)
table(is.na(gla_all$quantity_num))

#get traj
gla_t <- getTrajectoriesUKBgla(gla_all)


tab_gla_t <- sort(table(gla_t$drug_name))


#get median of pills_norm
gla_d <- gla_t %>%
#  filter(drug_name %in% drugs_to_keep) %>%
  group_by(drug_name) %>%
  summarise_at(vars(pills_norm), funs(median(., na.rm=TRUE)))

summary(gla_d$pills_norm)
table(round(gla_d$pills_norm))

# set dose variable
gla_d <- gla_d %>% filter(pills_norm < 5) %>%
  mutate(dose = ifelse(pills_norm < 1, 1, round(pills_norm)))

summary(gla_d)
table(gla_d$dose)

gla_all <- gla_all %>% 
  left_join(gla_d, by = "drug_name") %>%
  filter(!is.na(dose))   


#get dose-weighted traj 
gla_t <- getTrajectoriesUKBgla(gla_all, dose=gla_all$dose)

#summarize
gla_s <- summarizeTrajectoriesUKB(gla_t)


# distribution of adherence
p5 <- ggplot(data = gla_s, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  ggtitle('Overall for Glaucoma') +
  theme_minimal()

p5

#median adherence
median(gla_s$adherence)

#write data
fwrite(gla_t,"~/data/gla_traj.tsv",quote=F,sep="\t")
fwrite(gla_s,"~/data/gla_summary.tsv",quote=F,sep="\t")







