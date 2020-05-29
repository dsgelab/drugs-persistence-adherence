rm(list=ls())

library(dplyr)
library(data.table)
library(ggplot2)

dd <- fread('/home/cordioli/drugs/statins_10_15.tsv')

no_pkg <- dd %>%
  filter(CODE4 == 0) %>%
  pull(FINNGENID)

dd <- dd %>%
  filter(!FINNGENID %in% no_pkg)


# missing VNR:
# length(which(is.na(d$pkoko_num)))
# 262078

# length(unique(d$CODE3))
# 279

vnr_missing <- dd %>%
  filter(is.na(pkoko_num)) %>%
  select(CODE3, CODE1, pkoko_num, NAME) %>%
  distinct()
  
n_ind <- length(
  dd %>%
  ungroup() %>%
  filter(!is.na(pkoko_num)) %>%
  distinct(FINNGENID) %>%
  pull(FINNGENID)
)

n_purch <- length(
  dd %>%
    ungroup() %>%
    filter(!is.na(pkoko_num)) %>%
    pull(FINNGENID)
)

# plot VNR
vnr <- dd %>%
  ungroup() %>%
  filter(!is.na(pkoko_num)) %>%
  mutate(n_pkg = CODE4) %>%
  count(CODE3, CODE1, valmiste, vahvuus, pkoko_num, n_pkg) %>%
  mutate(tot_pkg = n_pkg*n) %>%
  group_by(CODE3, CODE1, valmiste, vahvuus, pkoko_num) %>%
  summarise(tot = sum(tot_pkg)) %>%
  mutate(per_purch = tot/n_purch,
         per_year = tot/(5*n_purch),
         per_pers = tot/(5*n_ind))

vnr <- dd %>%
  ungroup() %>%
  filter(!is.na(pkoko_num)) %>%
  mutate(n_pkg = CODE4) %>%
  count(CODE1, NAME, vahvuus, pkoko_num, n_pkg) %>%
  mutate(tot_pkg = n_pkg*n) %>%
  group_by(CODE1, NAME, vahvuus, pkoko_num) %>%
  summarise(tot = sum(tot_pkg)) %>%
  mutate(per_purch = tot/n_purch,
         per_year = tot/(5*n_purch),
         per_pers = tot/(5*n_ind),
         lab = paste0(CODE1, ' - ', NAME, ' - ', vahvuus, ' - ', pkoko_num))


ggplot(vnr, aes(x = lab, y = tot)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
  

# # # Apply RFM only on data with VNR

ids <- dd %>%
  filter(is.na(pkoko_num)) %>%
  pull(FINNGENID)

d <- dd %>%
  filter(!FINNGENID %in% ids)

# diff day
d <- d %>%
  arrange(FINNGENID, EVENT_AGE) %>%
  group_by(FINNGENID) %>%
  mutate(diff_day = round( (EVENT_AGE - lag(EVENT_AGE))*365.25 ))

# Split trajectories:
# - diff_day > Npkg * Npills

# Npills for each purchase
d$n_pills <- d$CODE4*d$pkoko_num
summary(d$pkoko_num)
summary(d$n_pills)

d <- d %>%
  arrange(FINNGENID, APPROX_EVENT_DAY) %>%
  group_by(FINNGENID) %>%
  mutate(new_traj = as.integer(diff_day > lag(n_pills) + 30))

d$new_traj[which(is.na(d$new_traj))] <- 1

d$diff_day2 <- d$diff_day
d$diff_day2[which(d$new_traj == 1)] <- NA

t <- d %>%
  group_by(FINNGENID) %>%
  summarise(tot_pills = sum(n_pills),
            SD_days = sd(diff_day2, na.rm = T),
            N_traj = sum(new_traj),
            N_purch = n(),
            t_len = sum(diff_day2, na.rm = T)) %>%
  mutate(freq = N_purch/t_len) %>%
  filter(N_purch > 4,
         !is.na(SD_days))



tm <- melt(t[,c(2,3,4,7)])

ggplot(tm, aes(x=value)) + 
  geom_histogram(bins = 50) +
  facet_wrap(~variable, scales = 'free')

cor(t[,c(2,3,4,7)])

# Variables distribution and correlations
library(ggplot2)
library(GGally)
ggpairs(t, columns=c(2,3,4,7)) + 
  theme_minimal() +
  ggtitle("Variables distribution and correlation")

hist(d$n_pills)

# Categorize variables
t <- t %>%
  mutate(tot_pills_c = cut(tot_pills, 5, include.lowest=TRUE, labels=c(1,2,3,4,5)),
         SD_days_c = cut(SD_days, 5, include.lowest=TRUE, labels=c(1,2,3,4,5)),
         N_traj_c = cut(N_traj, 5, include.lowest=TRUE, labels=c(1,2,3,4,5)),
         freq_c = cut(freq, 5, include.lowest=TRUE, labels=c(1,2,3,4,5))
         )


# # # K-Means
# function to compute total within-cluster sum of square 
# Compute and plot wss for k = 1 to k = 15
library(tidyverse)
library(cluster)
library(factoextra)

fviz_nbclust(t[,c(2,3,4,7)], kmeans, method = "wss")
# k = 4 seems the optimal for wss
clusters <- kmeans(t[,c(2,3,4,7)], 4, nstart = 25)
t$cl <- clusters$cluster
sil <- silhouette(clusters$cluster, dist(t[,c(2,3,4,7)]))
plot(sil, col=1:4, border=NA)

ggpairs(t, columns=c(2,3,4,7), aes(colour = as.factor(cl), alpha = .5)) + 
  theme_minimal() +
  ggtitle("")


# clustering on categorical variables
fviz_nbclust(t[,9:12], kmeans, method = "wss")
# k = 4 seems the optimal for wss
clusters <- kmeans(t[,9:12], 4, nstart = 25)
t$cl_c <- clusters$cluster
sil <- silhouette(clusters$cluster, dist(t[,9:12]))
plot(sil, col=1:4, border=NA)

ggpairs(t, columns=c(2,3,4,7), aes(colour = as.factor(cl_c), alpha = .5)) + 
  theme_minimal() +
  ggtitle("")

cl <- t %>%
  select(cl, tot_pills, SD_days, N_traj, freq) %>%
  group_by(cl) %>%
  summarise_all("mean")

cl_c <- t %>%
  select(cl_c, tot_pills, SD_days, N_traj, freq) %>%
  group_by(cl_c) %>%
  summarise_all("mean")

cl
cl_c


# Normalize data
t <- t %>%
  mutate(tot_pills_n = scale(tot_pills),
         SD_days_n = scale(SD_days),
         N_traj_n = scale(N_traj),
         freq_n = scale(freq)
  )
# clustering on norm variables
fviz_nbclust(t[,14:17], kmeans, method = "wss")
# k = 4 seems the optimal for wss
clusters <- kmeans(t[,14:17], 4, nstart = 25)
t$cl_c <- clusters$cluster
sil <- silhouette(clusters$cluster, dist(t[,14:17]))
plot(sil, col=1:4, border=NA)

ggpairs(t, columns=c(2,3,4,7), aes(colour = as.factor(cl_c), alpha = .5)) + 
  theme_minimal() +
  ggtitle("")


# Viz variables per age at first purchase, sex
t_age <- d %>%
  group_by(FINNGENID) %>%
  summarise(age_first = first(EVENT_AGE),
            tot_pills = sum(n_pills),
            SD_days = sd(diff_day2, na.rm = T),
            N_traj = sum(new_traj),
            N_purch = n(),
            t_len = sum(diff_day2, na.rm = T)) %>%
  mutate(tot_pill_n = tot_pills/t_len) %>%
  filter(N_purch > 4,
         !is.na(SD_days)) %>%
  mutate(age_bin = cut(age_first, 5, include.lowest=TRUE))

tt_age <- t_age %>%
  ungroup() %>%
  group_by(age_bin) %>%
  summarise(m_n = mean(tot_pill_n),
            m_tot = mean(tot_pills))

ggplot(tt_age, aes(x=age_bin, y=m_n)) + 
  geom_bar(stat="identity") +
  theme_minimal() +
  ylab("TOT PILLS / TOT DAYS")


# Plot trajectories and coloured by cluster
clusters_id <- t %>%
  select(FINNGENID, cl, cl_c)

d2 <- merge(d, clusters_id)

df <- d2 %>% 
  filter(FINNGENID %in% sample(unique(FINNGENID), 30)) %>%
  group_by(FINNGENID) %>%
  arrange(EVENT_AGE, cl_c)

ggplot(df,aes(x=as.Date(APPROX_EVENT_DAY),y=factor(FINNGENID),colour=FINNGENID,group=FINNGENID)) + 
  geom_point(aes(size=n_pills)) +
  scale_size_continuous(name="N pills", range = c(.6,3)) +
  geom_line() +
  theme_minimal() +
  scale_colour_discrete(guide = FALSE)

+
  theme(legend.position = "none")


df <- d2 %>% 
  filter(FINNGENID %in% sample(unique(FINNGENID), 30)) %>%
  group_by(FINNGENID) %>%
  arrange(EVENT_AGE, cl)

ggplot(df,aes(x=as.Date(APPROX_EVENT_DAY),y=FINNGENID,colour=factor(cl),group=FINNGENID)) + 
  geom_point(aes(size=n_pills)) +
  scale_size_continuous(name="N pills", range = c(.5,4)) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(breaks = c("1","2","3","4"),  values=c("red", "blue", "green", "black"))
