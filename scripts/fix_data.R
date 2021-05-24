rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

source('/home/cordioli/drugs/adherence_funs.R')

purch <- fread('/home/cordioli/drugs/data/finngen_R6_v2_purch_vnr_98.gz')

ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
                'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Complete trajectories
st <- getTrajectories(purch,'^C10AA')
ATCs <- '^C10AA'

#purch <- dat
# TODO: impute n_packages where 0 (comes like that from the registries). Just set it to 1 now
#purch$CODE4[purch$CODE4==0] <- 1
df <- purch %>%
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
  filter(grepl(ATCs, CODE1)) %>%
  group_by(FINNGENID) %>%
  # remove rows duplicated because they have different INDEX values in the registries
  select(-INDEX) %>%
  distinct() %>%
  # filter trajectories with at least for purchases
  filter(n()>=4) %>%
  # filter trajectories with VNR info (package size) for all the events
  filter(all(!is.na(pkoko_num))) %>%
  # for each individual, order by event_age
  arrange(EVENT_AGE, .by_group = T)  %>%
  # calculate n pills for each purchse
  mutate(n_pills = pkoko_num*CODE4,
         # compute difference between purchases
         days_next_purch = round((lead(EVENT_AGE) - EVENT_AGE)*365.25),
         change_type = CODE1 != lag(CODE1),
         # exclude gap >=1.5 years
         # days_next_purch = case_when(days_next_purch >= 365*1.5 ~ NA_real_,
         #                            TRUE ~ days_next_purch),
         # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
         # set to NA last purchase and 'gap' purchases (where days_nest is.na)
         n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                             lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                             TRUE ~ n_pills),
         pills_norm = n_pills/days_next_purch,
         days_norm = days_next_purch/n_pills)


# Get rows where days_next_purch==0 and the following one
dup_idx <- which(df$days_next_purch == 0 & !is.na(df$days_next_purch))
length(dup_idx)/nrow(df)*100
# 0.24%
dup_idx <- sort(c(dup_idx, dup_idx+1))
dup <- df[dup_idx,]
dup$pid <- seq(1,nrow(dup))

# Case 1: same age_event, same approx_day, different VNRs and strengths.
# gap to the next purchase suggest patient is taking 2 pills
# 6116 purchases (0.23%)
dup1 <- dup %>%
  filter(EVENT_AGE == lead(EVENT_AGE) & APPROX_EVENT_DAY == lead(APPROX_EVENT_DAY) & vnr != lead(vnr))
plus.one <- dup1$pid+1
dup1 <- rbind(dup1, dup[dup$pid %in% plus.one,])
dup1 <- dup1[order(dup1$pid),]

# Case 2: same event_age, different approx_event_day (we cannot resolve gap<3.6 days having age with 2 decimals), same VNR
# 292 (0.01%)
dup2 <- dup %>%
  filter(EVENT_AGE == lead(EVENT_AGE) & APPROX_EVENT_DAY != lead(APPROX_EVENT_DAY) & vnr == lead(vnr))
plus.one <- dup2$pid+1
dup2 <- rbind(dup2, dup[dup$pid %in% plus.one,])

# Case 3: same event_age, different approx_event_day (we cannot resolve gap<3.6 days having age with 2 decimals), differnt VNR
# 367 (0.01%)
dup3 <- dup %>%
  filter(EVENT_AGE == lead(EVENT_AGE) & APPROX_EVENT_DAY != lead(APPROX_EVENT_DAY) & vnr != lead(vnr))
plus.one <- dup3$pid+1
dup3 <- rbind(dup3, dup[dup$pid %in% plus.one,])


# cutoff to separate trajectories
ddf <- df 

# Days without pills
gaps <- ddf$days_next_purch - ddf$n_pills
gaps <- gaps[gaps>0 & !is.na(gaps)]

min(gaps)
max(gaps)

br <- seq(0,7250,by=50)
ranges <- paste(head(br,-1), br[-1], sep=" - ")

freq <- hist(gaps, breaks=br, include.lowest=TRUE, plot=FALSE)

dd <- data.frame(range = ranges, frequency = freq$counts)

hist(gaps, breaks=br, include.lowest=TRUE)


gaps <- gaps[gaps<601 & !is.na(gaps)]
br <- seq(0,600,by=10)
ranges <- paste(head(br,-1), br[-1], sep=" - ")
hist(gaps, breaks=br, include.lowest=TRUE)
freq <- hist(gaps, breaks=br, include.lowest=TRUE, plot=FALSE)
dd <- data.frame(range = ranges, frequency = freq$counts)

gaps <- gaps[gaps<101 & !is.na(gaps)]
br <- seq(0,100,by=10)
hist(gaps, breaks=br, include.lowest=TRUE)


summary(ddf$n_pills)
n_pills <- ddf$n_pills[!is.na(ddf$n_pills)]
summary(n_pills)
hist(n_pills)
