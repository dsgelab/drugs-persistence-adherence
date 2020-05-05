rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

dat <- fread('/home/cordioli/drugs/R5_v3_purch_vnr_98.gz')
endpoints <- fread('/home/cordioli/R5_pheno/finngen_R5_v3_endpoint.gz')

# # Define 'chronic' users: events that lead to a strong need of statins:
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

# Extract purchases of ATC of interest
# Extract individuals with at least 4 purchases
atc <- '^C10AA'
df <- dat %>%
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
  filter(grepl(atc, CODE1)) %>%
  group_by(FINNGENID) %>%
  select(-INDEX) %>%
  distinct() %>%
  filter(n()>=4) %>%
  arrange(FINNGENID, EVENT_AGE)

# Compute days_next_purch
df <- df %>%
  mutate(days_next_purch = round( (lead(EVENT_AGE) - EVENT_AGE)*365.25)) %>%
  filter(days_next_purch != 0 | is.na(days_next_purch))

# Tag new_start: if lag(days_next) >= lag(n_pkg*max_pill + gap), new_start=1
# Try different values
df <- df %>%
  mutate(start_100_60 = case_when(is.na(days_next_purch) ~ 0,
                                  lag(days_next_purch) >= lag(CODE4)*100 + 60 ~ 1,
                                  TRUE ~ 0),
         start_100_90 = case_when(is.na(days_next_purch) ~ 0,
                                  lag(days_next_purch) >= lag(CODE4)*100 + 90 ~ 1,
                                  TRUE ~ 0),
         start_100_120 = case_when(is.na(days_next_purch) ~ 0,
                                  lag(days_next_purch) >= lag(CODE4)*100 + 120 ~ 1,
                                  TRUE ~ 0),
         start_150_60 = case_when(is.na(days_next_purch) ~ 0,
                                  lag(days_next_purch) >= lag(CODE4)*150 + 60 ~ 1,
                                  TRUE ~ 0),
         start_150_90 = case_when(is.na(days_next_purch) ~ 0,
                                  lag(days_next_purch) >= lag(CODE4)*150 + 90 ~ 1,
                                  TRUE ~ 0),
         start_200_30 = case_when(is.na(days_next_purch) ~ 0,
                                  lag(days_next_purch) >= lag(CODE4)*200 + 30 ~ 1,
                                  TRUE ~ 0))

# Summarize each trajectory
t <- df %>%
  group_by(FINNGENID) %>%
  summarise(tot_days = sum(days_next_purch, na.rm = T),
            mean_days = mean(days_next_purch, na.rm = T),
            sd_days = sd(days_next_purch, na.rm = T),
            tot_purch = n(),
            age_first_purch = first(EVENT_AGE),
            gaps_100_60 = sum(start_100_60),
            gaps_100_90 = sum(start_100_90),
            gaps_100_120 = sum(start_100_120),
            gaps_150_60 = sum(start_150_60),
            gaps_150_90 = sum(start_150_90),
            gaps_200_30 = sum(start_200_30)) %>%
  filter_all(all_vars(!is.na(.))) %>%  # if individual has just one purchase I'll have NA 
  left_join(chronics) %>%
  left_join(sex) %>%
  mutate(age_bin = cut(age_first_purch, 5),
         chronic = case_when(!is.na(age_cvd_ev) ~ 1,
                             TRUE ~ 0),
         start_after_cvd_ev = case_when(age_cvd_ev <= age_first_purch ~ 1,
                                        TRUE ~ 0)) %>%
  mutate_at(vars(starts_with("gaps_")), 
            .funs = list(binary = ~case_when(. > 0 ~ 1, TRUE ~ 0))) %>%
  mutate_at(vars(starts_with("gaps_") & !ends_with('_bin')), 
            .funs = list(categoric = ~case_when(. > 5 ~ 5, TRUE ~ .)))

# Some plots
p1 <- ggplot(t, aes(x=factor(binary_disc)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  theme_minimal()
p2 <- ggplot(t, aes(x=factor(multi_disc)))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme_minimal()
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)


ggplot(t, aes(x=factor(binary_disc)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  facet_grid(.~age_bin) +
  theme_minimal()


# Merge phenotypes with covariates for GWAS
covs <- fread('/home/cordioli/R5_pheno/R5_cov_pheno_1.0.txt.gz')

covs <- covs[,1:grep("PC20", colnames(covs))]

t$SEX <- NULL
covs <- covs %>%
  left_join(t, by = 'FINNGENID')

fwrite(covs, '/home/cordioli/drugs/discontinuation_cov_pheno.txt', sep = '\t', quote = F)

phenolist_binary <- colnames(covs)[grepl('_binary$', colnames(covs))]
phenolist_categoric <- colnames(covs)[grepl('_categoric$', colnames(covs))]

fwrite(list(phenolist_binary), '/home/cordioli/drugs/discontinuation_binary_phenolist.txt', col.names = F)
fwrite(list(phenolist_categoric), '/home/cordioli/drugs/discontinuation_categoric_phenolist.txt', col.names = F)