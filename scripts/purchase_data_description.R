rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

purch <- fread('/home/cordioli/drugs/data/finngen_R7_purch_vnr_98.gz')

# system("gcsfuse --only-dir finngen_R6/phenotype_2.0/data/ --file-mode 444 finngen-production-library-red /home/cordioli/mount_fg/")

d <- fread('/home/cordioli/drugs/data/finngen_R7_cov.txt')

colnames(purch)

# Tot individuals
length(unique(purch$FINNGENID))
# 258,281 / 271,112: 95%

# R7:
# > length(unique(purch$FINNGENID))
# [1] 305428
#  
# > nrow(d)
# [1] 309154
#  
# > length(unique(purch$FINNGENID))/nrow(d)
# [1] 0.9879478


# Tot purchases
nrow(purch)
# 54,140,783
# R7: 63,084,035

# # # ATC counts
atc_map <- fread('drugs/data/ATC_translate_984.tsv') %>%
  select(Id, Longname)

ATCs_count <- data.frame(table(purch$CODE1)) %>%
  left_join(atc_map, by = c("Var1" = "Id")) %>%
  arrange(desc(Freq))

# subgroup: 4 digits, level 3, therapeutic/pharmacological subgroup
ATCs_subgroup <- substr(purch$CODE1,1,4)
ATCs_groups_count <- data.frame(table(ATCs_subgroup)) %>%
  left_join(atc_map, by = c("ATCs_subgroup" = "Id")) %>%
  arrange(desc(Freq))

# plot first 20 meds, first 20 groups
top30 <- ATCs_count[1:20,]
top30$Longname <- factor(top30$Longname, levels = top30$Longname)

pdf('/home/cordioli/drugs/plots/top20_ATCs.pdf', width = 6, height = 4)
ggplot(data=top30, aes(x=reorder(Longname, Freq), y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  coord_flip() +
  theme_minimal()
dev.off()

top10 <- ATCs_groups_count[1:20,]
top10$Longname <- factor(top10$Longname, levels = top10$Longname)

pdf('/home/cordioli/drugs/plots/top20_ATCgroups.pdf', width = 6, height = 4)
ggplot(data=top10, aes(x=reorder(Longname, Freq), y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  coord_flip() +
  theme_minimal()
dev.off()

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


pdf('/home/cordioli/drugs/plots/R7_purch_years.pdf', width = 7, height = 5)
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('/home/cordioli/drugs/plots/R7_purch_months.pdf', width = 7, height = 5)
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>%
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

pdf('/home/cordioli/drugs/plots/R7_purch_avg_months.pdf', width = 7, height = 5)
ggplot(data=purch_av_month_year, aes(x=mm, y=avg_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()
