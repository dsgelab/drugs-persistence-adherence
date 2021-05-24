rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

traj <- fread('/home/cordioli/drugs/data/UKB_statins_summarized.txt')

# keep only individuls and purchases included in adherence calculation

purch <- fread('/home/cordioli/drugs/data/UKB_statins_gap_150.txt') %>%
  filter(eid %in% traj$eid,
         !is.na(days_next_purch))

# # ATCs included
# ATCs <- purch %>%
#   select(CODE1, valmiste) %>%
#   mutate(valmiste = sapply(strsplit(valmiste, " "), "[", 1)) %>%
#   distinct() %>%
#   arrange(CODE1)
#   
# table(traj$chronic)

# # # Purchases over time: years/months
yy <- as.numeric(format(purch$issue_date,'%Y'))
mm <- as.numeric(format(purch$issue_date,'%m'))

purch_year <- data.frame(table(yy))

purch_month_year <- data.frame(table(mm, yy))

purch_month_year <- purch_month_year %>%
  mutate(Freq_mm_yy = Freq) %>%
  select(-Freq) %>%
  filter(yy != 2020) %>%
  left_join(purch_year)


pdf('/home/cordioli/drugs/plots/ukbb/UKB_statins_purch_years.pdf', width = 7, height = 5)
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

pdf('/home/cordioli/drugs/plots/ukbb/UKB_statins_purch_months.pdf', width = 7, height = 5)
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>%
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

pdf('/home/cordioli/drugs/plots/ukbb/UKB_statins_purch_avg_months.pdf', width = 7, height = 5)
ggplot(data=purch_av_month_year, aes(x=mm, y=avg_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()


summary(traj$tot_pills)
summary(traj$tot_days)
summary(traj$tot_purch)
summary(purch$quantity_num)

# # adherence over time (per month)
# keep only individuls and purchases included in adherence calculation
pp <- purch %>%
  select(eid, issue_date, n_pills, days_next_purch) %>%
  mutate(year = as.numeric(format(purch$issue_date,'%Y')),
         month = as.numeric(format(purch$issue_date,'%m'))) %>%
  group_by(eid, year, month) %>%
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

pdf('/home/cordioli/drugs/plots/ukbb/UKB_statins_adherence_purchases_per_month.pdf', width = 4, height = 3)
ggplot(data=melted, aes(x=month, y=value)) +
  geom_line(aes(group = 1, colour = variable)) +
  geom_point() +
  theme_minimal() +
  facet_grid(variable~., scales = "free") +
  theme(legend.position = "none")
dev.off()


# Adherence
source('~/drugs/adherence_funs.R')

traj$chronic <- 0
p <- plotAdherenceDensity(traj)

p2 <- plotAdherenceByAge(traj)
str(traj)
