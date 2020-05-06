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

# -Adjust pills for each purchase: if I change type of statins after e.g. 20 days, and I've bought 100 pills, I should count only 20 pills for that purch
# -Compute pill/day
# -Tag gap/discontinuation: days_next_purchase > n_pills + 60/90/120
st <- st %>%
  mutate(change_type = CODE1 != lag(CODE1),
         n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                             lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                             TRUE ~ n_pills),
        pills_norm = n_pills/days_next_purch,
        start_60 = case_when(is.na(days_next_purch) ~ 0,
                             lag(days_next_purch) >= lag(n_pills) + 60 ~ 1,
                             TRUE ~ 0),
        start_90 = case_when(is.na(days_next_purch) ~ 0,
                             lag(days_next_purch) >= lag(n_pills) + 90 ~ 1,
                             TRUE ~ 0),
        start_120 = case_when(is.na(days_next_purch) ~ 0,
                             lag(days_next_purch) >= lag(n_pills) + 120 ~ 1,
                             TRUE ~ 0))

# Summarize each trajectory
t <- st %>%
  group_by(FINNGENID) %>%
  summarise(tot_days = sum(days_next_purch, na.rm = T),
            mean_days = mean(days_next_purch, na.rm = T),
            sd_days = sd(days_next_purch, na.rm = T),
            tot_purch = n(),
            age_first_purch = first(EVENT_AGE),
            gaps_60 = sum(start_60),
            gaps_90 = sum(start_90),
            gaps_120 = sum(start_120)) %>%
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
  mutate_at(vars(starts_with("gaps_") & !ends_with('_binary')), 
            .funs = list(categoric = ~case_when(. > 5 ~ 5, TRUE ~ .)))

# Some plots
p1 <- ggplot(t, aes(x=factor(gaps_60_binary)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  theme_minimal()
p2 <- ggplot(t, aes(x=factor(gaps_90_binary)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  theme_minimal()
p3 <- ggplot(t, aes(x=factor(gaps_120_binary)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  theme_minimal()
p4 <- ggplot(t, aes(x=factor(gaps_60_categoric)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  theme_minimal()
p5 <- ggplot(t, aes(x=factor(gaps_90_categoric)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  theme_minimal()
p6 <- ggplot(t, aes(x=factor(gaps_120_categoric)))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme_minimal()

library(gridExtra)
grid.arrange(p1, p2, p3,
             p4, p5, p6, nrow = 2, ncol = 3)

# distribution by age
p1 <- ggplot(t, aes(x=factor(gaps_60_binary)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  facet_grid(.~age_bin) +
  theme_minimal()
p2 <- ggplot(t, aes(x=factor(gaps_90_binary)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  facet_grid(.~age_bin) +
  theme_minimal()
p3 <- ggplot(t, aes(x=factor(gaps_120_binary)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  facet_grid(.~age_bin) +
  theme_minimal()
p4 <- ggplot(t, aes(x=factor(gaps_60_categoric)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  facet_grid(.~age_bin) +
  theme_minimal()
p5 <- ggplot(t, aes(x=factor(gaps_90_categoric)))+
  geom_bar(stat="count", width=0.7, fill="steelblue") +
  facet_grid(.~age_bin) +
  theme_minimal()
p6 <- ggplot(t, aes(x=factor(gaps_120_categoric)))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  facet_grid(.~age_bin) +
  theme_minimal()

grid.arrange(p1, p2, p3,
             p4, p5, p6, nrow = 2, ncol = 3)



# Merge phenotypes with covariates for GWAS
cov_pheno <- fread('/home/cordioli/R5_pheno/R5_cov_pheno_1.0.txt.gz')
cov_pheno <- cov_pheno[,1:grep("PC20", colnames(cov_pheno))]

t$SEX <- NULL
cov_pheno <- cov_pheno %>%
  left_join(t, by = 'FINNGENID')

fwrite(cov_pheno, '/home/cordioli/drugs/data/discontinuation_cov_pheno_VNR0605.txt', sep = '\t', quote = F)

phenolist_binary <- colnames(cov_pheno)[grepl('_binary$', colnames(cov_pheno))]
phenolist_categoric <- colnames(cov_pheno)[grepl('_categoric$', colnames(cov_pheno))]

fwrite(list(phenolist_binary), '/home/cordioli/drugs/data/discontinuation_binary_phenolist_VNR0605.txt', col.names = F)
fwrite(list(phenolist_categoric), '/home/cordioli/drugs/data/discontinuation_categoric_phenolist_VNR0605.txt', col.names = F)

# Summarise number case-control
N <- 66350
ncases <- cov_pheno %>%
  filter(!is.na(age_first_purch)) %>%
  select(ends_with('_binary')) %>%
  summarise_each(sum)

nctrl <- N - ncases
