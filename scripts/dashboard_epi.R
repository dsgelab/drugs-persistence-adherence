rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

traj <- fread('/home/cordioli/drugs/data/R6_epi_summarized_150.txt') %>%
  filter(chronic == 1)

# keep only individuls and purchases included in adherence calculation
purch <- fread('/home/cordioli/drugs/data/R6_epi_gap_150.txt') %>%
  filter(FINNGENID %in% traj$FINNGENID,
         !is.na(days_next_purch))
  
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


pdf('/home/cordioli/drugs/plots/epi_purch_years.pdf', width = 7, height = 5)
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

pdf('/home/cordioli/drugs/plots/epi_purch_months.pdf', width = 7, height = 5)
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>%
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

pdf('/home/cordioli/drugs/plots/epi_purch_avg_months.pdf', width = 7, height = 5)
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

pdf('/home/cordioli/drugs/plots/epi_adherence_month_avg_across_ind.pdf', width = 7, height = 5)
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

pdf('/home/cordioli/drugs/plots/epi_adherence_month_avg_per_ind.pdf', width = 7, height = 5)
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

pdf('/home/cordioli/drugs/plots/epi_adherence_purchases_per_month.pdf', width = 4, height = 3)
ggplot(data=melted, aes(x=month, y=value)) +
  geom_line(aes(group = 1, colour = variable)) +
  geom_point() +
  theme_minimal() +
  facet_grid(variable~., scales = "free") +
  theme(legend.position = "none")
dev.off()


# Adherence
require(cowplot)
require(dplyr)
require(readr)
require(grid)
require(gridExtra)
source("/home/cordioli/R/RainCloudPlots/tutorial_R/R_rainclouds.R")
palette <-  c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data = traj, aes(y = adherence, x = 0))+
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8, colour = "dodgerblue4", fill = "dodgerblue3") +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5, colour = "dodgerblue4", fill = "dodgerblue3") +
  coord_flip() +
  scale_x_discrete(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  scale_fill_manual(values=palette) +
  theme_minimal() +
  theme(legend.position = "none")

palette <-  c("#999999", "#E69F00", "#009E73","#0072B2", "#56B4E9", "#D55E00", "#CC79A7", "#F0E442")
ggplot(data = traj, aes(y = adherence, x = factor(age_bin), fill = factor(age_bin))) +
  geom_flat_violin(position = position_nudge(x = 0, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_discrete(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1, 2)) +
  scale_fill_manual(values=palette) +
  theme_minimal() +
  theme(legend.position = "none")
