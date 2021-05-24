rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

traj <- fread('/home/cordioli/drugs/data/R7_breast_cancer_summarized_150.txt')

# keep only individuls and purchases included in adherence calculation
purch <- fread('/home/cordioli/drugs/data/R7_breast_cancer_gap_150.txt') %>%
  filter(FINNGENID %in% traj$FINNGENID,
         !is.na(days_next_purch))

# ATCs included
ATCs <- purch %>%
  select(CODE1) %>%
  distinct() %>%
  arrange(CODE1)
  
table(traj$chronic)

# # # Purchases over time: years/months
yy <- as.numeric(format(purch$APPROX_EVENT_DAY,'%Y'))
mm <- as.numeric(format(purch$APPROX_EVENT_DAY,'%m'))

purch_year <- data.frame(table(yy))

purch_month_year <- data.frame(table(mm, yy))

purch_month_year <- purch_month_year %>%
  mutate(Freq_mm_yy = Freq) %>%
  select(-Freq) %>%
  filter(yy != 2020) %>%
  left_join(purch_year)


pdf('/home/cordioli/drugs/plots/R7/breast_cancer_purch_years.pdf', width = 7, height = 5)
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal()
dev.off()

pdf('/home/cordioli/drugs/plots/R7/breast_cancer_purch_months.pdf', width = 7, height = 5)
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>%
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

pdf('/home/cordioli/drugs/plots/R7/breast_cancer_purch_avg_months.pdf', width = 7, height = 5)
ggplot(data=purch_av_month_year, aes(x=mm, y=avg_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()


summary(traj$tot_pills)
summary(traj$tot_days)
summary(traj$tot_purch)
summary(purch$pkoko_num)

# # adherence over time (per month)
# keep only individuls and purchases included in adherence calculation
pp <- purch %>%
  select(FINNGENID, APPROX_EVENT_DAY, n_pills, days_next_purch) %>%
  mutate(year = as.numeric(format(purch$APPROX_EVENT_DAY,'%Y')),
         month = as.numeric(format(purch$APPROX_EVENT_DAY,'%m'))) %>%
  group_by(month) %>%
  summarise(adh = sum(n_pills, na.rm = T) / sum(days_next_purch, na.rm = T))
  
adh_avg_month <- pp %>%
  group_by(month) %>%
  summarise(avg_adh_month = mean(adh))

pdf('/home/cordioli/drugs/plots/R7/breast_cancer_adherence_month_avg_across_ind.pdf', width = 7, height = 5)
ggplot(data=adh_avg_month, aes(x=month, y=avg_adh_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()

pp <- purch %>%
  select(FINNGENID, APPROX_EVENT_DAY, n_pills, days_next_purch) %>%
  mutate(year = as.numeric(format(purch$APPROX_EVENT_DAY,'%Y')),
         month = as.numeric(format(purch$APPROX_EVENT_DAY,'%m'))) %>%
  group_by(FINNGENID, year, month) %>%
  summarise(adh = sum(n_pills, na.rm = T) / sum(days_next_purch, na.rm = T))

adh_avg_month <- pp %>%
  group_by(month) %>%
  summarise(avg_adh_month = mean(adh))

pdf('/home/cordioli/drugs/plots/R7/breast_cancer_adherence_month_avg_per_ind.pdf', width = 7, height = 5)
ggplot(data=adh_avg_month, aes(x=month, y=avg_adh_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()


# combine N purchases and adherence by month on the same plot
adh_avg_month$month <- factor(adh_avg_month$month)

adh_purch_avg_month <- adh_avg_month %>%
  left_join(purch_av_month_year, by = c("month" = "mm")) %>%
  select(month, mean_adherence = avg_adh_month, mean_purchases = avg_month)

melted <- melt(adh_purch_avg_month)

pdf('/home/cordioli/drugs/plots/R7/breast_cancer_adherence_purchases_per_month.pdf', width = 4, height = 3)
ggplot(data=melted, aes(x=month, y=value)) +
  geom_line(aes(group = 1, colour = variable)) +
  geom_point() +
  theme_minimal() +
  facet_grid(variable~., scales = "free") +
  theme(legend.position = "none")
dev.off()