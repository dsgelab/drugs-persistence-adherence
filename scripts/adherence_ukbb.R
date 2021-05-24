rm(list=ls())

library(ggplot2)
library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)

source('/home/cordioli/drugs/adherence_funs.R')

write.flag <- T

# # # Pre-process
# read_excel('drugs/data/ukbb/BNF Snomed Mapping data 20210218.xlsx') %>%
#   as.data.frame() %>%
#   select(bnf_code = `BNF Code`, bnf_name = `BNF Name`, SNOMED_code = `SNOMED Code`, strength = Strength,
#          unit = `Unit Of Measure`, pkg_size = `Pack`) %>%
#   fwrite('drugs/data/ukbb/bnf_snomed_mapping_processed.tsv', sep = '\t')

# read in files
# scripts <- fread('~/drugs/data/ukbb/ukb31063.gp_scripts.20200706w.processed.bnf_ok.date_ok.txt')

# f.31.0.0: sex
# f.53.0.0: Date of attending assessment centre
# f.21022.0.0: age at recruitment

cov <- fread('~/drugs/data/ukbb/export_for_map_analysis.tsv') %>% 
  select(eid = f.eid, sex = f.31.0.0, date_recruitment = f.53.0.0, age_recruitment = f.21022.0.0)


# # # # # # # #
#   STATINS   #
# # # # # # # #

# Define 'chronic' users: events that lead to a strong need of statins:
# I9_ASO
# I9_CHD Coronary heart disease
# I9_ATHSCLE Atherosclerosis
# I9_CEREBVASC	Cerebrovascular diseases	
# I9_INTRACRA	Nontraumatic intracranial haemmorrhage	
# I9_SAH	Subarachnoid haemmorrhage	
# I9_ICH	Intracerebral haemmorrhage	
# I9_OTHINTRACRA	Other intracranial haemorrhages	
# I9_STR	Stroke, excluding SAH	
# I9_STR_SAH	Stroke, including SAH	
# I9_STR_EXH	Ischaemic Stroke, excluding all haemorrhages	
# I9_STENOSIS	Occlusion and stenosis of arteries, not leading to stroke
# ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
#                 'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')


# Complete trajectories
st_list <- fread('~/drugs/data/ukbb/ukbb_statins.tsv') %>% 
  pull(drug_name)

st_all <- fread('~/drugs/data/ukbb/UKB_all_statins.tsv')
length(unique(st_all$eid))

st_all <- st_all %>% 
  left_join(cov, by = "eid")

st <- getTrajectoriesUKB(st_all, st_list)
length(unique(st$eid))

if (write.flag == T) {
  fwrite(st, '/home/cordioli/drugs/data/UKB_statins_gap_150_after_enrolment.txt', sep = '\t', quote = F)
}

# Summarised trajectories
stt <- summarizeTrajectoriesUKB(st)
stt$chronic <- 0

if (write.flag == T) {
  fwrite(stt, '/home/cordioli/drugs/data/UKB_statins_summarized_after_enrolment.txt', sep = '\t', quote = F)
}

plotAdherenceDensity(stt)
plotAdherenceByAge(stt)


# # # # # # # # # # # # # # # # # #  
#   BLOOD PRESSURE MEDICATIONS    #
# # # # # # # # # # # # # # # # # #

# Define 'chronic' users: events that lead to a strong need of bp:

# TODO: add endpoints
# ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
#                 'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Individual trajectories
bp <- getTrajectories(purch,'^C0(2|3|8|9)')

if (write.flag == T) {
  fwrite(bp, '/home/cordioli/drugs/data/R6_bp_gap_150.txt', sep = '\t', quote = F)
}

# Summarised trajectories
bpp <- summarizeTrajectories(bp)

bpp <- bpp %>%
  left_join(covs)

# Difference between sexes
t.test(bpp$adherence[bpp$SEX_IMPUTED==0], bpp$adherence[bpp$SEX_IMPUTED==1])$p.value


cor(bpp$adherence,bpp$sd_days)
cor.test(bpp$adherence,bpp$sd_days)$p.value

summary(lm(adherence ~ sd_days, data = bpp))

# h2: statins:  0.0242, 2.19E-03. Bp:  0.0394, 1.49E-07
# Age first related event
# age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

# stt <- stt %>%
#   left_join(age_first) %>%
#   mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
#                              TRUE ~ 0))

if (write.flag == T) {
  fwrite(bpp, '/home/cordioli/drugs/data/R6_bp_summarized_150.txt', sep = '\t', quote = F)
}


# # # # # # # # # # # # # # # # # #  
#   CLOPIDOGREL & ASA+DIP   #
# # # # # # # # # # # # # # # # # #
# clopi '^B01AC04'
# dipydogrol '^B01AC30' -> adjust pills *0.5 (dose is 2pills/day)

# Define 'chronic' users: events that lead to a strong need of bp:

ep_chronic <- c('I9_STR_SAH', 'I9_TIA', 'I9_MI')

ATCs <- data.frame(atc = c("B01AC04","B01AC30"),
                   dose = c(1,2),
                   stringsAsFactors = F)

# Complete trajectories
cd <- getTrajectories2(purch,ATCs)

if (write.flag == T) {
  fwrite(cd, '/home/cordioli/drugs/data/R6_clopi_dipy_gap_150.txt', sep = '\t', quote = F)
}

# Summarised trajectories
ccdd <- summarizeTrajectories(cd)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

ccdd <- ccdd %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

ccdd <- ccdd %>%
  left_join(covs)

p <- plotAdherenceDensity(ccdd)


if (write.flag == T) {
  fwrite(ccdd, '/home/cordioli/drugs/data/R6_clopi_dipy_summarized_150.txt', sep = '\t', quote = F)
}



# # # # # # # # # # #
#   BREAST CANCER   #
# # # # # # # # # # #

# If initiated, the medications should continue to at least  ~5 years after the diagnosis date. Some cancers require a longer period, but that is the general minimum.
# They usually start within 1 year after the diagnosis (after chemotherapy etc is done), and a successful therapy could therefore be defined as at least 4 years of purchases. The C3_BREAST diagnosis in FinnGen endpoints is a good one.
# These four are ATC-kodes that are relevant, and they should be combined to one variable, because in many cases, one can be switched to another:
# L02BA01 Tamoxifen
# L02BG04 Letrozole
# L02BG06 Exemestane
# L02BG03 Anastrozole
# In brief, the timeline would be
# If a person has started of any of the 4 medications -> looking at whether the purchases discontinue before 4 years. Some discontinuation may of course be that the breast cancer has recurred and requires new aggressive treatments,
# how have you handled such cases for other diseases?
 
ep_chronic <- c('C3_BREAST')

# Complete trajectories
bc <- getTrajectories(purch,'^L02B(A01|G04|G06|G03)')

if (write.flag == T) {
  fwrite(bc, '/home/cordioli/drugs/data/R6_breast_cancer_gap_150.txt', sep = '\t', quote = F)
}

# Summarised trajectories
bcc <- summarizeTrajectories(bc)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, bcc$FINNGENID, ep_chronic)

bcc <- bcc %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

# # # 
# Exclude non-chronic and re-standardize adherence for GWAS
bcc <- bcc %>%
  filter(chronic == 1) %>%
  mutate(adherence_std = as.numeric(scale(adherence)))

if (write.flag == T) {
  fwrite(bcc, '/home/cordioli/drugs/data/R6_breast_cancer_summarized_150.txt', sep = '\t', quote = F)
}


# # # # # # # # # #
#   ASTHMA MEDS   #
# # # # # # # # # #

ep_chronic <- c('J10_ASTHMA')

# Complete trajectories
as <- getTrajectories(purch,'^R0(3AK|3AC|3BA|1AD)')

# Summarised trajectories
aas <- summarizeTrajectories(as)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, aas$FINNGENID, ep_chronic)

aas <- aas %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

aas <- aas %>%
  left_join(covs)

# Visualize pkoko distribution
pkoko <- as %>%
  ungroup() %>%
  count(pkoko)

ggplot(pkoko, aes(x = reorder(pkoko, n), y = n)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal()


# Difference between sexes
t.test(aas$adherence[aas$SEX_IMPUTED==0], aas$adherence[aas$SEX_IMPUTED==1], alternative = 'greater')$p.value

# Primary VS secondary prevention
t.test(aas$adherence[aas$chronic==0], aas$adherence[aas$chronic==1], paired = F)

# Regularity
cor(aas$adherence,aas$sd_days)
cor.test(aas$adherence,aas$sd_days)


if (write.flag == T) {
  fwrite(aas, '/home/cordioli/drugs/data/R6_asthma_summarized.txt', sep = '\t', quote = F)
}

# # # 
# Merge phenotypes with covariates for GWAS
# ad <- aas %>%
#   mutate(statins = adherence_std,
#          age_first_statins = age_first_purch) %>%
#   select(FINNGENID,statins,age_first_statins)
# 
# covs <- covs %>%
#   left_join(ad)
# 
# if (write.flag == T) {
#   fwrite(covs, '/home/cordioli/drugs/data/R6_cov_pheno_adherence.txt', sep = '\t', quote = F)
# }



# # # # # # # #
#   DOAC      #
# # # # # # # #

# Define 'chronic' users: events that lead to a strong need of statins:

ep_chronic <- c('I9_DVTANDPULM','I9_VTE','I9_AF')

ATCs <- data.frame(atc = c("B01AF01","B01AF02","B01AF03","B01AE07"),
                   dose = c(1,2,1,2),
                   stringsAsFactors = F)

# Complete trajectories
do <- getTrajectories2(purch,ATCs)

if (write.flag == T) {
  fwrite(do, '/home/cordioli/drugs/data/R6_doac_gap_150.txt', sep = '\t', quote = F)
}

# Summarised trajectories
doo <- summarizeTrajectories(do, 1.1)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, doo$FINNGENID, ep_chronic)

doo <- doo %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

doo <- doo %>%
  left_join(covs)

p <- plotAdherenceDensity(doo)

# Difference between sexes
t.test(doo$adherence[doo$SEX_IMPUTED==0], doo$adherence[doo$SEX_IMPUTED==1], alternative = 'greater')$p.value

# Primary VS secondary prevention
t.test(doo$adherence[doo$chronic==0], doo$adherence[doo$chronic==1], paired = F)

# Regularity
cor(doo$adherence,doo$sd_days)
cor.test(doo$adherence,doo$sd_days)


if (write.flag == T) {
  fwrite(doo, '/home/cordioli/drugs/data/R6_doac_summarized_150.txt', sep = '\t', quote = F)
}



# # # # 
# # Merge phenotypes with covariates for GWAS
# ad <- doo %>%
#   mutate(statins = adherence_std,
#          age_first_statins = age_first_purch) %>%
#   select(FINNGENID,statins,age_first_statins)
# 
# covs <- fread('/home/cordioli/drugs/data/R6_cov_pheno_adherence.txt')
# 
# covs <- covs %>%
#   left_join(ad)
# 
# if (write.flag == T) {
#   fwrite(covs, '/home/cordioli/drugs/data/R6_cov_pheno_adherence.txt', sep = '\t', quote = F)
# }


# # # # # # # #
#   GLAUCOMA  #
# # # # # # # #

# # # Some pkokoNum need manual check becasue they are wrong
# gl <- purch %>%
#   mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
#   filter(grepl('^S01E(E|D)', CODE1))
# 
# pkoko_list <- c("1X2,5 ML", "3 x 3 ml", "3X5ML", "30 x 0.4 g", "90 x 0.3 ml", "30 x 0.2 ml", "90 x 0.2 ml", "90 (18 x 5) x 0.2 ml", "3 x 2.5 ml")
# 
# getPkokoNum <- function(pkoko) {
#   b <- strsplit(gsub(" |ml|g","",tolower(gsub(",", "\\.", pkoko))), "x")
#   return( unlist(lapply(b, function(i) as.numeric(i[1])*as.numeric(i[2]))) )
# }
# 
# gl$pkoko_num[gl$pkoko %in% pkoko_list] <- getPkokoNum(gl$pkoko[gl$pkoko %in% pkoko_list])
# 
# pkoko <- gl %>%
#   ungroup() %>%
#   select(pkoko, pkoko_num) %>%
#   distinct
# 
# fwrite(gl, '/home/cordioli/drugs/data/R6_glaucoma_all_pkoko_ok.txt', sep = '\t', quote = F)

# Define 'chronic' users: events that lead to a strong need of the medication:

ep_chronic <- c('H7_GLAUCPRIMOPEN')

# Complete trajectories
g <- fread('/home/cordioli/drugs/data/R6_glaucoma_all_pkoko_ok.txt')

gl <- getTrajectories(g,'^S01E(E|D)',0.1)

plot(density(gl$pills_norm, na.rm = T))

getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(gl$pills_norm)

# Summarised trajectories
gll <- summarizeTrajectories(gl, 1.1)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, gll$FINNGENID, ep_chronic)

gll <- gll %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

gll <- gll %>%
  left_join(covs)

p <- plotAdherenceDensity(gll)

# Difference between sexes
t.test(gll$adherence[gll$SEX_IMPUTED==0], gll$adherence[gll$SEX_IMPUTED==1], alternative = 'greater')$p.value

# Primary VS secondary prevention
t.test(gll$adherence[gll$chronic==0], gll$adherence[gll$chronic==1], paired = F)

# Regularity
cor(gll$adherence,gll$sd_days)
cor.test(gll$adherence,gll$sd_days)

if (write.flag == T) {
  fwrite(gl, '/home/cordioli/drugs/data/R6_glauc_gap_150.txt', sep = '\t', quote = F)
}

if (write.flag == T) {
  fwrite(gll, '/home/cordioli/drugs/data/R6_glauc_summarized_150.txt', sep = '\t', quote = F)
}


# # # # # # # # # # #
#   EPILEPSY        #
# # # # # # # # # # #

m <- fread('/home/cordioli/drugs/data/epilepsy_meds.txt')
meds <- m$drug_code[m$perc_epilepsy>=40]
# meds <- substring(paste(sub('N03A','',meds), collapse = "|"),2)

meds_regex <- "A02|B02|B52|D01|F02|F03|F04|G04|G06|X14|X15|X17|X18|X22|X23"

ep_chronic <- c('G6_EPLEPSY')

# Complete trajectories
# '^C0(2|3|8|9)')
epi <- getTrajectories(purch, paste0('^N03A(',meds_regex,')'))

# Summarised trajectories
epii <- summarizeTrajectories(epi, 10)
# Age first related event
age_first <- getAgeFirstEndpoint(ep, epii$FINNGENID, ep_chronic)

epii <- epii %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

plotAdherenceDensity(epii)

# Peaks at 1, 2, 4
unique(epi$CODE1[between(epi$pills_norm,0.9,1.1)])
unique(epi$CODE1[between(epi$pills_norm,1.9,2.1)])
unique(epi$CODE1[between(epi$pills_norm,3.9,4.1)])


# ppl disgnosed with G6_EPLEPSY
epi_diagnosed <- sum(ep$G6_EPLEPSY, na.rm = T)

# ppl with G6_EPLEPSY in the trajectories
length(intersect(epi_diagnosed, epii$FINNGENID))


epi <- epi %>%
  filter(FINNGENID %in% epii$FINNGENID)


if (write.flag == T) {
  fwrite(epi, '/home/cordioli/drugs/data/R6_epi_gap_150.txt', sep = '\t', quote = F)
}

if (write.flag == T) {
  fwrite(epii, '/home/cordioli/drugs/data/R6_epi_summarized_150.txt', sep = '\t', quote = F)
}