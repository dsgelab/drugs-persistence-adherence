####SET####

rm(list=ls())

library(ggplot2)
library(data.table)
library(lubridate)
library(readxl)
library(dplyr)
library(Hmisc)
require(cowplot)
require(dplyr)
require(readr)
require(grid)
require(gridExtra)

source("~/code/R_rainclouds.R")
source('~/code/adherence_funs.R')

f <- function(x){
  list(mean(x),median(x),sd(x))
}
get_stat <- function(df){
  ds <- df %>% select(tot_pills,tot_days,tot_purch,age_first_purch,adherence)
  stat <- sapply(ds,f)
  rownames(stat) <- c("mean","median","std_dev")
  return(stat)
}



#### DOAC ####
doac <- fread("~/data/doac_pheno_final.tsv")
#### stats ####
length(doac$eid)

doac_stat <- get_stat(doac)

doac_sex <- group_by(doac,sex) %>% 
  summarise(mean = mean(adherence),
            median = median(adherence),
            std_dev = sd(adherence),
            n = n())

doac_good <- group_by(doac,isgood) %>% 
  summarise(med_age = median(age_first_purch),
            med_days = median(tot_days),
            med_purch = median(tot_purch),
            n = n())


##### correlations 
cor(doac_sum$adherence,doac_sum$age_first_purch)
cor(doac_sum$adherence,doac_sum$tot_days)

#### plots ####
# # # Purchases over time: years/months
yy <- as.numeric(format(doac_purch$issue_date,'%Y'))
mm <- as.numeric(format(doac_purch$issue_date,'%m'))

purch_year <- data.frame(table(yy))

purch_month_year <- data.frame(table(mm, yy))

purch_month_year <- purch_month_year %>%
  mutate(Freq_mm_yy = Freq) %>%
  select(-Freq) %>%
  filter(yy != 2020) %>%
  left_join(purch_year)


png('/home/andreacorbetta/plots/doac_purch_years.png')
ggplot(data=purch_year, aes(x=yy, y=Freq)) +
  geom_bar(stat="identity", fill = "dodgerblue3") +
  theme_minimal() +
  ggtitle('Purchases per years - DOAC') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

png('/home/andreacorbetta/plots/doac_purch_months.png')
ggplot(data=purch_month_year, aes(x=mm, y=Freq_mm_yy/Freq, colour = yy)) +
  geom_line(aes(group = yy)) +
  ggtitle('Purchases per month - DOAC') +
  theme_minimal()
dev.off()

purch_av_month_year <- purch_month_year %>%
  group_by(mm) %>%
  summarise(avg_month = mean(Freq_mm_yy),
            avg_prop_month = mean(Freq_mm_yy/Freq)) 

#plot avg month purchases
png('/home/andreacorbetta/plots/doac_purch_avg_months.png')
ggplot(data=purch_av_month_year, aes(x=mm, y=avg_month)) +
  geom_line(aes(group = 1)) +
  geom_point() +
  ggtitle('Average purchases per month - DOAC') +
  theme_minimal()
dev.off()





# # adherence over time (per month)
# keep only individuls and purchases included in adherence calculation
pp <- doac_purch %>%
  select(eid, issue_date, n_pills, days_next_purch) %>%
  mutate(year = as.numeric(format(doac_purch$issue_date,'%Y')),
         month = as.numeric(format(doac_purch$issue_date,'%m'))) %>%
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

melted <- melt(as.data.table(adh_purch_avg_month))

#plot purch and adherence"
png('/home/andreacorbetta/plots/doac_adherence_purchases_per_month.png')
ggplot(data=melted, aes(x=month, y=value)) +
  geom_line(aes(group = 1, colour = variable)) +
  geom_point() +
  ggtitle('Adherence and purchases - DOAC') +
  theme_minimal() +
  facet_grid(variable~., scales = "free") +
  theme(legend.position = "none")
dev.off()


# Adherence


#plot overall
doac_sum$chronic <- 0
p <- plotAdherenceDensity(doac_sum,title="DOAC")
png('/home/andreacorbetta/plots/doac_adherence_overall.png',height=240)
p
dev.off()


#plot summary
p2 <- plotAdherenceByAge(doac_sum,title="DOAC")
png("/home/andreacorbetta/plots/doac_adherence_byage.png")
p2
dev.off()

###plot oonly by age
p3 <- ggplot(doac_sum, aes(x=age_bin, y=adherence, fill=age_bin)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  ggtitle('Adherence by age DOAC') +
  theme_minimal() +
  theme(plot.title=element_text(size=10)) +
  theme(legend.position = "none")

png('/home/andreacorbetta/plots/doac_adherence_only_byage.png',width=240)
p3
dev.off()


## new stuff

da_doac <- ggplot(data=doac_sum, aes(x=tot_days, y=adherence)) +
  geom_point(aes(colour = factor(sex)),size=0.7) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "blue") +
  scale_y_continuous(limits=c(0,1.1),breaks = c(0.0,0.3,0.6,0.9,1.1)) +
  theme_minimal() +
  ggtitle('adh vs days - DOAC') +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))


age_doac <- ggplot(data=doac_sum, aes(x=age_first_purch, y=adherence)) +
  geom_point(aes(colour = factor(sex)),size=1) +
  scale_y_continuous(limits=c(0,1.1),breaks = c(0.0,0.3,0.6,0.9,1.1)) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "blue") +
  theme_minimal() +
  ggtitle('adh vs age - DOAC') +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))






#### AP ####
ap <- fread("~/data/antiplat_pheno_final.tsv")
bc <- fread("~/data/breastcanc_pheno_final.tsv")
bp <- fread("~/data/bloodpres_pheno_final.tsv")
gla <- fread("~/data/glaucoma_pheno_final.tsv")
stat <- fread("~/data/statins_pheno_final.tsv")


#### COMBINE ####


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(da_doac)

ggarrange(da_doac + ggtitle("") + rremove("legend"),
          da_bp + ggtitle("") + rremove("legend"),
          da_bc  + ggtitle("") + rremove("legend"),
          da_ap  + ggtitle("") + rremove("legend"), 
          da_gla  + ggtitle("") + rremove("legend"),
          mylegend,
          labels = c("DOAC", "BP", "BC", "AP","GLA"),
          ncol = 2, nrow = 3)


mylegend<-g_legend(age_doac)

ggarrange(age_doac + ggtitle("") + rremove("legend"),
          age_bp + ggtitle("") + rremove("legend"),
          age_bc  + ggtitle("") + rremove("legend"),
          age_ap  + ggtitle("") + rremove("legend"), 
          age_gla  + ggtitle("") + rremove("legend"),
          mylegend,
          labels = c("DOAC", "BP", "BC", "AP","GLA"),
          ncol = 2, nrow = 3)

