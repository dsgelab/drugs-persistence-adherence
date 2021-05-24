rm(list=ls())

library(ggplot2)
library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)

# # drop individuals in withdraw list and write latest version
w <- fread('~/drugs/data/ukbb/w31063_20210201.csv')

scripts1 <- fread('~/drugs/data/ukbb/ukb31063.gp_scripts.20191008.txt') %>%
  filter(!eid %in% w$V1)
fwrite(scripts1,'~/drugs/data/ukbb/ukb31063.gp_scripts.20191008w.txt', sep = '\t')

scripts2 <- fread('~/drugs/data/ukbb/ukb31063.gp_scripts.20200706.txt') %>%
  filter(!eid %in% w$V1)
fwrite(scripts2,'~/drugs/data/ukbb/ukb31063.gp_scripts.20200706w.txt', sep = '\t')

# Differences between the two versions:
# length(unique(scripts1$eid))
# [1] 222105
# length(unique(scripts2$eid))
# [1] 222105
# length(intersect(unique(scripts1$eid),unique(scripts2$eid)))
# [1] 222105

# we can use the one with latest date

scripts <- scripts2

str(scripts)
# Classes ‘data.table’ and 'data.frame':	57706541 obs. of  8 variables:
# $ eid          : int  1000035 1000035 1000035 1000035 1000035 1000035 1000035 1000035 1000035 1000035 ...
# $ data_provider: int  3 3 3 3 3 3 3 3 3 3 ...
# $ issue_date   : chr  "16/01/2006" "26/04/2006" "04/06/2008" "08/12/2010" ...
# $ read_2       : chr  "" "" "" "" ...
# $ bnf_code     : chr  "14.04.05.00.00" "05.01.01.03.00" "13.04.02.00.00" "01.06.04.00.00" ...
# $ dmd_code     : integer64 NA NA NA NA NA NA NA NA ... 
# $ drug_name    : chr  "Revaxis vaccine suspension for injection 0.5ml pre-filled syringes (sanofi pasteur MSD Ltd)" "Amoxicillin 500mg capsules" "Clobetasone 0.05% cream" "Macrogol compound oral powder sachets NPF sugar free" ...
# $ quantity     : chr  "1 pre-filled syringe(s) - 0.5 ml pre-filled syringe" "21 capsule(s) - 500 mg" "1 pack of 30 gram(s)" "1 pack of 30 sachet(s)" ...
# - attr(*, ".internal.selfref")=<externalptr> 

nrow(scripts)
# 57,704,460

length(unique(scripts$eid))
# 222,105 indidivuals

# process BNF codes and dates
scripts <- scripts %>%
  mutate(bnf_code = gsub("\\.", "", bnf_code),
         # process dates with lubridates
         issue_date = lubridate::dmy(issue_date))

fwrite(scripts,'~/drugs/data/ukbb/ukb31063.gp_scripts.20200706w.processed.txt', sep = '\t')

# min - max dates
summary(scripts$issue_date)

# # # check how many entries per each code (read_2, BNF, dmd)
# N prescription annotated with BNF code: 43,770,600 (76%)
sort(unique(scripts$bnf_code))
n_bnf <- scripts %>%
  filter(!bnf_code %in% c("", "00000000")) %>%
  nrow()
n_bnf/nrow(scripts)

#Dates:
# where clinical event or prescription date precedes participant date of birth it has been altered to
# 01/01/1901.
# Where the date matches participant date of birth it has been altered to 02/02/1902.
# Where the date follows participant date of birth but is in the year of their birth it has been
# altered to 03/03/1903.
# Where the date was in the future this has been changed to 07/07/2037 as these are likely to
# have been entered as a place-holder or other system default.

scripts_bnf_dateok <- scripts %>%
  filter(!bnf_code %in% c("", "00000000")) %>% 
  filter(between(issue_date,as.Date("1903-03-04"),as.Date("2037-07-06")))%>%
  select(eid, data_provider, issue_date, bnf_code, drug_name, quantity)

fwrite(scripts_bnf_dateok,'~/drugs/data/ukbb/ukb31063.gp_scripts.20200706w.processed.bnf_ok.date_ok.txt', sep = '\t')

nrow(scripts_bnf)-nrow(scripts_bnf_dateok)
# 7275 records with date not valid

summary(scripts_bnf_dateok$issue_date)
# Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
# "1945-07-26" "2006-09-05" "2011-01-19" "2009-10-27" "2014-02-19" "2019-08-25" 


# # # Prescription over time: years/months
yy <- as.numeric(format(scripts_bnf_dateok$issue_date,'%Y'))
mm <- as.numeric(format(scripts_bnf_dateok$issue_date,'%m'))

purch_year <- data.frame(table(yy))

purch_month_year <- data.frame(table(mm, yy))

purch_month_year <- purch_month_year %>%
  mutate(Freq_mm_yy = Freq) %>%
  select(-Freq) %>%
  filter(yy != 2020) %>%
  left_join(purch_year)


pdf('/home/cordioli/drugs/plots/UKBB_purch_years.pdf', width = 10, height = 5)
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

# restrict to 2002 on to be comparable with finnene

purch_month_year <- purch_month_year %>% 
  mutate(year = as.integer(as.character(yy))) %>% 
  filter(year >= 2002)

pdf('/home/cordioli/drugs/plots/UKBB_purch_months_2002.pdf', width = 7, height = 5)
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>% 
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

pdf('/home/cordioli/drugs/plots/UKBB_purch_avg_months_2002.pdf', width = 7, height = 5)
ggplot(data=purch_av_month_year, aes(x=mm, y=avg_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  theme_minimal()
dev.off()



# N prescription annotated with read_2 and not BNF: 13,761,249
sort(unique(scripts$read_2))
n_read2 <- scripts %>%
  filter(bnf_code %in% c("", "00000000") & read_2 != "") %>%
  nrow()
n_read2 / nrow(scripts)


# N prescription annotated only with dmd: 123,444
n_dmd <- scripts %>%
  filter(bnf_code %in% c("", "00000000") & read_2 == "" & !is.na(dmd_code) & dmd_code != 0) %>%
  nrow()

n_dmd / nrow(scripts)


# N prescription annotated with a BNF code present in the mapping: 1,925,358
scripts %>%
  filter(bnf_code %in% bnf_map$`BNF Code`) %>%
  nrow()

# BNF mapping
bnf_map <- read_excel('drugs/data/ukbb/BNF Snomed Mapping data 20210218.xlsx') %>%
  as.data.frame()

length(unique(scripts$bnf_code))
# 5204
length(unique(bnf_map$`BNF Code`))
# 49171
length(intersect(unique(scripts$bnf_code), unique(bnf_map$`BNF Code`)))
# 2816

summary(nchar(scripts$bnf_code))
summary(nchar(bnf_map$`BNF Code`))

head(scripts)
head(bnf_map)

# BNF - ATC mapping
bnf_atc_map <- fread('drugs/data/bnf_atc_map.tsv')

# # # 
bfs <- c('^02120', '^020[2|4|5]|^020602','^060102') #bnf codes for 1) lipid lowering 2) blood pressure 3) T2D medicines

# STATINS: ATC ^C10AA
statins_ATCs <- c("C10AA01", "C10AA02", "C10AA03", "C10AA04", "C10AA05", "C10AA06", "C10AA07")

# check mapping
bnf_atc_map %>%
  filter(atc %in% statins_ATCs)

# 1: 2120000 C10AA05
# 2: 2120000 C10AA07
# 3: 2120000 C10AA01
# 4: 2120000 C10AA04
# 5: 2120000 C10AA03

statins <- scripts %>%
  filter(grepl('^02120', bnf_code))

statins2 <- statins %>%
  inner_join(bnf_map, by = c("bnf_code" = "BNF Code"))

sort(unique(statins$drug_name))

head(statins)

sort(unique(statins$quantity))


# BP
load('drugs/data/r6_BP_ATCs')
bp_ATCs <- ATCs

bnf_atc_map %>%
  filter(atc %in% bp_ATCs) %>%
  arrange(atc)

length(intersect(bnf_atc_map$atc, bp_ATCs))


# clopi - dipy
cldy_ATCs <- c("B01AC04", "B01AC30")

bnf_atc_map %>%
  filter(atc %in% cldy_ATCs) %>%
  arrange(atc)


# DOAC
doac_ATCs <- c("B01AF01", "B01AF02", "B01AF03", "B01AF07")

bnf_atc_map %>%
  filter(atc %in% doac_ATCs) %>%
  arrange(atc)


# Breast cancer

bc_ATCs <- c("L02BA01", "L02BG03", "L02BG04", "L02BG06")

bnf_atc_map %>%
  filter(atc %in% bc_ATCs) %>%
  arrange(atc)


# Glaucoma
# ^S01E(E|D)

bnf_atc_map %>%
  filter(grepl("^S01E(E|D)", atc)) %>%
  arrange(atc)
