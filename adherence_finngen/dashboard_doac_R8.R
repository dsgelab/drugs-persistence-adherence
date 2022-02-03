rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

setwd('/home/ivm/')

traj <- fread('drugs/data/R8_doac_summarized.txt')
traj.lt75 <- fread('drugs/data/R8_doac_summarized_lt75.txt')

# keep only individuals and purchases included in adherence calculation
purch <- fread('drugs/data/R8_doac.txt') %>%
  filter(FINNGENID %in% traj$FINNGENID,
         !is.na(days_next_purch))

# stats
nrow(traj)
table(traj$chronic)
table(traj$chronic)[2]/nrow(traj)

sum(traj$tot_purch)

nrow(traj)/338183

summary(traj$tot_pills)
summary(traj$tot_days)
summary(traj$tot_purch)
summary(purch$pkoko_num)

mean(traj$adherence)
sd(traj$adherence)
sum(traj$adherence>0.8)/nrow(traj)
 # 2nd prevention
traj_2 <- traj %>% 
  filter(chronic == 1)
mean(traj_2$adherence)
sd(traj_2$adherence)
sum(traj_2$adherence>0.8)/nrow(traj_2)

cor(traj$adherence,traj$tot_days)
cor(traj$adherence,traj$age_first_purch)


# # ATCs included
ATCs <- purch %>%
  select(CODE1, valmiste) %>%
  mutate(valmiste = sapply(strsplit(valmiste, " "), "[", 1)) %>%
  distinct() %>%
  arrange(CODE1)


# Adherence distribution, all vs truncated at 75y
source('/home/ivm/drugs/scripts/adherence_funs.R')

pdf('/home/ivm/drugs/plots/R8/doac_adherence.pdf', width = 4, height = 6)
plotAdherenceDensity(traj)
dev.off()

# plotAdherenceDensity(traj.lt75)

pdf('/home/ivm/drugs/plots/R8/doac_adherence_age.pdf', width = 8, height = 6)
plotAdherenceByAge(traj)
dev.off()

# plotAdherenceByAge(traj.lt75)


# # # Purchases over time: years/months
yy <- as.numeric(format(purch$APPROX_EVENT_DAY,'%Y'))
mm <- as.numeric(format(purch$APPROX_EVENT_DAY,'%m'))

# assign all purchases in Jan 2020 to Dec 2019 (this is due to date randomization, end of register is in 31.12.2019)
yy[yy==2020] <- 2019
mm[yy==2020] <- 12


purch_year <- data.frame(table(yy))

purch_month_year <- data.frame(table(mm, yy))

purch_month_year <- purch_month_year %>%
  mutate(Freq_mm_yy = Freq) %>%
  select(-Freq) %>%
  filter(yy != 2020) %>%
  left_join(purch_year)


pdf('/home/ivm/drugs/plots/R8/doac_purch_years.pdf', width = 7, height = 5)
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('/home/ivm/drugs/plots/R8/doac_purch_months.pdf', width = 7, height = 5)
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>%
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

pdf('/home/ivm/drugs/plots/R8/doac_purch_avg_months.pdf', width = 7, height = 5)
ggplot(data=purch_av_month_year, aes(x=mm, y=avg_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()


# # adherence over time (per month)
# keep only individuls and purchases included in adherence calculation
pp <- purch %>%
  select(FINNGENID, APPROX_EVENT_DAY, n_pills, days_next_purch) %>%
  mutate(year = as.numeric(format(purch$APPROX_EVENT_DAY,'%Y')),
         month = as.numeric(format(purch$APPROX_EVENT_DAY,'%m'))) %>%
  group_by(FINNGENID, year, month) %>%
  summarise(adh = sum(n_pills, na.rm = T) / sum(days_next_purch, na.rm = T))
  
adh_avg_month <- pp %>%
  group_by(month) %>%
  summarise(avg_adh_month = mean(adh))


# combine N purchases and adherence by month on the same plot
adh_avg_month$month <- factor(adh_avg_month$month)

adh_purch_avg_month <- adh_avg_month %>%
  left_join(purch_av_month_year, by = c("month" = "mm")) %>%
  select(month, mean_adherence = avg_adh_month, mean_purchases = avg_month)

melted <- melt(adh_purch_avg_month)

pdf('/home/ivm/drugs/plots/R8/doac_adherence_purchases_per_month.pdf', width = 4, height = 3)
ggplot(data=melted, aes(x=month, y=value)) +
  geom_line(aes(group = 1, colour = variable)) +
  geom_point() +
  theme_minimal() +
  facet_grid(variable~., scales = "free") +
  theme(legend.position = "none")
dev.off()

