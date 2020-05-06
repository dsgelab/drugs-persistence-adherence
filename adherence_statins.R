rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

dat <- fread('/home/cordioli/drugs/data/R5_v3_purch_vnr_98.gz')
endpoints <- fread('/home/cordioli/R5_pheno/finngen_R5_v3_endpoint.gz')

# # # # #
# Define 'chronic' users: events that lead to a strong need of statins:
# I9_ASO
# I9_CHD
# I9_ATHSCLE
# I9_CEREBVASC	Cerebrovascular diseases	
# I9_INTRACRA	Nontraumatic intracranial haemmorrhage	
# I9_SAH	Subarachnoid haemmorrhage	
# I9_ICH	Intracerebral haemmorrhage	
# I9_OTHINTRACRA	Other intracranial haemorrhages	
# I9_STR	Stroke, excluding SAH	
# I9_STR_SAH	Stroke, including SAH	
# I9_STR_EXH	Ischaemic Stroke, excluding all haemorrhages	
# I9_STENOSIS	Occlusion and stenosis of arteries, not leading to stroke
ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
                'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')
ep_chronic <- intersect(ep_chronic, colnames(endpoints))
ep_chronic_age <- paste0(ep_chronic, '_AGE')

# Select people with onset of one of the endpoints and age of onset
chronics <- endpoints %>%
  select(FINNGENID, all_of(ep_chronic), all_of(ep_chronic_age)) %>%
  filter_at(vars(-FINNGENID), any_vars(.==1)) %>%
  select(FINNGENID, all_of(ep_chronic_age))

# For each row, the min accross columns is the age of onset of the first of the CVD events
min_age <- apply(chronics[,2:length(ep_chronic)], 1, min)
chronics <- chronics %>%
  mutate(age_cvd_ev = min_age) %>%
  select(FINNGENID, age_cvd_ev)

sex <- endpoints %>%
  select(FINNGENID, SEX)


# # # # #
# Select statins
# Remove duplicated rows (same event with a different INDEX value)
# Keep trajectories with at least 4 purchases
# Keep trajectories with VNR info for all purchases
st <- dat %>%
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
  filter(grepl('^C10AA', CODE1)) %>%
  group_by(FINNGENID) %>%
  select(-INDEX) %>%
  distinct() %>%
  filter(n()>=4) %>%
  filter(all(!is.na(pkoko_num))) %>%
  arrange(FINNGENID, EVENT_AGE)
length(unique(st$FINNGENID))

# Compute days to next purchase
# Compute n pills as pkoko_num * n pkg
# TODO: Combine multiple purchase on the same day summing number of pills (where everything but n pills is the same)
# st_new <- st %>%
#   group_by_at(setdiff(names(st), c('n_pills', 'days_next_purch'))) %>%
#   summarise(n_pills=sum(n_pills),
#             days_next_purch = last(days_next_purch))
# For the moment, just filter out all purch. with days_next == 0
# TODO: impute n_packages where 0 (30 cases). just set it to 1 now
st$CODE4[st$CODE4==0] <- 1
st <- st %>%
  mutate(n_pills = pkoko_num*CODE4,
         days_next_purch = round( (lead(EVENT_AGE) - EVENT_AGE)*365.25)) %>%
  filter(days_next_purch != 0 | is.na(days_next_purch))

# - Adjust pills for each purchase: if I change type of statins after e.g. 20 days, and I've bought 100 pills, I should count only 20 pills for that purch
# - TODO: Tag gap/discontinuation, define discontinuation better here: for looking at adherence better to have a wider gap, i wouldn't count that purch
# otherwise and it would bias the adherence
# - Compute pill/day
st <- st %>%
  mutate(change_type = CODE1 != lag(CODE1),
         days_next_purch = case_when(days_next_purch >= 365*1.5 ~ NA_real_,
                                     TRUE ~ days_next_purch),
         n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                             lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                             TRUE ~ n_pills),
         pills_norm = n_pills/days_next_purch)

# Summarise adherence per each trajectory as tot pills/tot days
t <- st %>%
  group_by(FINNGENID) %>%
  summarise(tot_pills = sum(n_pills, na.rm = T),
            tot_days = sum(days_next_purch, na.rm = T),
            mean_days = mean(days_next_purch, na.rm = T),
            sd_days = sd(days_next_purch, na.rm = T),
            mean_pills_norm = mean(pills_norm, na.rm = T),
            sd_pills_norm = sd(pills_norm, na.rm = T),
            tot_purch = n()-1,
            age_first_purch = first(EVENT_AGE)) %>%
  mutate(adherence = tot_pills/tot_days) %>%
  filter_all(all_vars(!is.na(.))) %>%
  left_join(chronics) %>%
  mutate(chronic = case_when(!is.na(age_cvd_ev) ~ 1,
                             TRUE ~ 0),
         after_cvd = case_when(age_cvd_ev <= age_first_purch ~ 1,
                               TRUE ~ 0),
         adherence_bin = case_when(adherence >= 0.8 ~ 1,
                                   TRUE ~ 0),
         adherence_cat = case_when(adherence < 0.8 ~ 1,
                                   between(adherence, 0.8, 1.2) ~ 2,
                                   TRUE ~ 3))

# Normal
t_norm <- t %>%
  filter(adherence <=1.2)

# Overbuyers
t_overb <- t %>%
  filter(adherence > 1.2)


# ppl taking half pill a day:
# - adherence 0.45 - 0.55
# - mean(pills_norm) 0.45 - 0.55
# - SD(pills_norm) < 0.25
t_half <- t %>%
  filter(between(adherence, .45, .55),
         between(mean_pills_norm, .45, .55),
         sd_pills_norm <= 0.25)


library(cowplot)
library(dplyr)
library(readr)
source("R/RainCloudPlots/tutorial_R/R_rainclouds.R")

ggplot(data = t, aes(y = adherence, x = factor(chronic), fill = factor(chronic))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  theme_minimal()

ggplot(data = t_norm, aes(y = adherence, x = factor(chronic), fill = factor(chronic))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  theme_minimal()


ggplot(t, aes(x=factor(chronic), y=adherence, group=chronic)) + 
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~SEX)

ggplot(t[t$adherence <= 1,], aes(x=factor(chronic), y=adherence, group=chronic)) + 
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~SEX)

ggplot(data = t[t$adherence <= 1,], aes(y = adherence, x = factor(chronic), fill = factor(chronic))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  theme_minimal() +
  facet_wrap(~SEX)


ggplot(t[t$adherence <= 1,], aes(x=age_bin, y=adherence, group=age_bin)) + 
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~factor(chronic))

table(t_norm$age_bin)

# # # Infer daily dose, examples
case_half <- st[st$FINNGENID == 'FGCRS4JHU2', ]
d <- density(case_half$days_norm, na.rm = T)
plot(d)
d <- density(case_half$pills_norm, na.rm = T)
plot(d)
case_one <- st[st$FINNGENID == 'FG334EMYMN', ]
d <- density(case_one$pills_norm, na.rm = T)
plot(d)
case_oneb <- st[st$FINNGENID == 'FGJN6YVWVN', ]
d <- density(case_oneb$pills_norm, na.rm = T)
plot(d)
d <- density(case_half$pills_norm, na.rm = T)
plot(d)


# # # 

# Merge phenotypes with covariates for GWAS
cov_pheno <- fread('/home/cordioli/R5_pheno/R5_cov_pheno_1.0.txt.gz')
cov_pheno <- cov_pheno[,1:grep("PC20", colnames(cov_pheno))]

t$SEX <- NULL
cov_pheno <- cov_pheno %>%
  left_join(t, by = 'FINNGENID')

fwrite(cov_pheno, '/home/cordioli/drugs/data/adherence_cov_pheno.txt', sep = '\t', quote = F)

phenolist_binary <- colnames(cov_pheno)[grepl('_bin$', colnames(cov_pheno))]
phenolist_categoric <- c('adherence', 'adherence_cat')

fwrite(list(phenolist_binary), '/home/cordioli/drugs/data/adherence_binary_phenolist.txt', col.names = F)
fwrite(list(phenolist_categoric), '/home/cordioli/drugs/data/adherence_quanti_phenolist.txt', col.names = F)
