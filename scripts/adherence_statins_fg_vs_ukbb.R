rm(list=ls())

library(ggplot2)
library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)

source('/home/cordioli/drugs/adherence_funs.R')

# Finngen
fg_t <- fread('~/drugs/data/R7_statins_summarized_150.txt')
fg_p <- fread('~/drugs/data/R7_statins_gap_150.txt') %>% 
  filter(FINNGENID %in% fg_t$FINNGENID)

# UKB, all scripts
ukb_t <- fread('~/drugs/data/UKB_statins_summarized.txt')
ukb_p <- fread('~/drugs/data/UKB_statins_gap_150.txt') %>% 
  filter(eid %in% ukb_t$eid)


# UKB, filter adfter enrolment
ukb_t2 <- fread('~/drugs/data/UKB_statins_summarized_after_enrolment.txt')
ukb_p2 <- fread('~/drugs/data/UKB_statins_gap_150_after_enrolment.txt') %>% 
  filter(eid %in% ukb_t2$eid)

min(ukb_p$issue_date)
min(ukb_p2$issue_date)

# Tot users
nrow(fg_t)
nrow(ukb_t)
nrow(ukb_t2)

# Tot scripts
sum(fg_t$tot_purch)
sum(ukb_t$tot_purch)
sum(ukb_t2$tot_purch)

# Mean scripts per individuals
summary(fg_t$tot_purch)
summary(ukb_t$tot_purch)
summary(ukb_t2$tot_purch)


summary(fg_t$tot_pills)
summary(ukb_t$tot_pills)
summary(ukb_t2$tot_pills)


mean(fg_t$tot_days)/365.25
mean(ukb_t$tot_days)/365.25
mean(ukb_t2$tot_days)/365.25


median(fg_p$pkoko_num)
median(ukb_p$quantity_num)
median(ukb_p2$quantity_num)

# Correlations
# adherence - age
cor(fg_t$adherence, fg_t$age_first_purch)
cor(ukb_t$adherence, ukb_t$age_first_purch)
cor(ukb_t2$adherence, ukb_t2$age_first_purch)
# adherence - days
cor(fg_t$adherence, fg_t$tot_days)
cor(ukb_t$adherence, ukb_t$tot_days)
cor(ukb_t2$adherence, ukb_t2$tot_days)
# adherence - "regularity"
cor(fg_t$adherence, fg_t$sd_days)
cor(ukb_t$adherence, ukb_t$sd_days)
cor(ukb_t2$adherence, ukb_t2$sd_days)

# Mean adherence by sex
mean(fg_t$adherence)
sd(fg_t$adherence)

mean(fg_t$adherence[fg_t$SEX == "female"])
sd(fg_t$adherence[fg_t$SEX == "female"])
mean(fg_t$adherence[fg_t$SEX == "male"])
sd(fg_t$adherence[fg_t$SEX == "male"])

# annotate ukb_t with sex
ukb_sex <- ukb_p %>% 
  select(eid, sex) %>% 
  distinct()
ukb_sex2 <- ukb_p2 %>% 
  select(eid, sex) %>% 
  distinct()
ukb_t <- ukb_t %>% 
  left_join(ukb_sex)
ukb_t2 <- ukb_t2 %>% 
  left_join(ukb_sex2)

# ukb sex: 0 female, 1 male
mean(ukb_t$adherence)
sd(ukb_t$adherence)
mean(ukb_t$adherence[ukb_t$sex == 0], na.rm = T)
sd(ukb_t$adherence[ukb_t$sex == 0], na.rm = T)
mean(ukb_t$adherence[ukb_t$sex == 1], na.rm = T)
sd(ukb_t$adherence[ukb_t$sex == 1], na.rm = T)


mean(ukb_t2$adherence)
sd(ukb_t2$adherence)
mean(ukb_t2$adherence[ukb_t$sex == 0], na.rm = T)
sd(ukb_t2$adherence[ukb_t$sex == 0], na.rm = T)
mean(ukb_t2$adherence[ukb_t$sex == 1], na.rm = T)
sd(ukb_t2$adherence[ukb_t$sex == 1], na.rm = T)


# Deltas between purchases/scripts, need to re-calculate the intervals because some were set to NA (long gaps)
fg_p <- fg_p %>% 
  group_by(FINNGENID) %>% 
  mutate(delta = round( (lead(EVENT_AGE) - EVENT_AGE)*365.25),
         delta_no_pills = delta - pkoko_num)

ukb_p <- ukb_p %>% 
  group_by(eid) %>% 
  mutate(delta = as.numeric(difftime(lead(issue_date), issue_date, units = "d")),
         delta_no_pills = delta - quantity_num)

ukb_p2 <- ukb_p2 %>% 
  group_by(eid) %>% 
  mutate(delta = as.numeric(difftime(lead(issue_date), issue_date, units = "d")),
         delta_no_pills = delta - quantity_num)

# Days without pills
gaps <- df$days_next_purch - df$n_pills
gaps <- gaps[gaps>0 & gaps<601 & !is.na(gaps)]

hist(gaps, breaks=br, include.lowest=TRUE, main=file_name)


# plots
require(grid)
require(gridExtra)

p1 <- ggplot(fg_p, aes(x=delta)) +
  geom_histogram(bins = 10, fill = "white", color="steelblue4", alpha=0.9) +
  ggtitle("Finngen, Statins - Days between purchases") +
  theme_minimal()
  
p2 <- ggplot(ukb_p2, aes(x=delta)) +
  geom_histogram(bins = 10, fill = "white", color="steelblue4", alpha=0.9) +
  ggtitle("UKBB, Statins - Days between prescriptions") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)

# plot filtering up to 1000 days
fg_p_500 <- fg_p %>% 
  filter(delta <= 500)

ukb_p2_500 <- ukb_p2 %>% 
  filter(delta <= 500)

p1 <- ggplot(fg_p_500, aes(x=delta)) +
  geom_histogram(binwidth = 50, fill = "white", color="steelblue4", alpha=0.9) +
  ggtitle("Finngen, Statins - Days between purchases") +
  theme_minimal()

p2 <- ggplot(ukb_p2_500, aes(x=delta)) +
  geom_histogram(binwidth = 50, fill = "white", color="steelblue4", alpha=0.9) +
  ggtitle("UKBB, Statins - Days between prescriptions") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)


p1 <- ggplot(fg_p_500[fg_p_500$delta_no_pills >= 0,], aes(x=delta_no_pills)) +
  geom_histogram(binwidth = 50, fill = "white", color="steelblue4", alpha=0.9) +
  ggtitle("Finngen, Statins - Days wo/ pills") +
  theme_minimal()

p2 <- ggplot(ukb_p2_500[ukb_p2_500$delta_no_pills >= 0,], aes(x=delta_no_pills)) +
  geom_histogram(binwidth = 50, fill = "white", color="steelblue4", alpha=0.9) +
  ggtitle("UKBB, Statins - Days wo/ pills") +
  theme_minimal()

grid.arrange(p1, p2, ncol = 2)
