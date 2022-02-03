rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)


setwd('/home/ivm/')

source('drugs/scripts/adherence_funs.R')

purch <- fread('drugs/data/finngen_R8_purch_vnr_98.gz')
ep <- fread('/finngen/library-red/finngen_R8/phenotype_3.0/data/finngen_R8_endpoint.txt.gz')

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
ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
                'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Complete trajectories
st <- getTrajectoriesDates(purch,'^C10AA')

fwrite(st, 'drugs/data/R8_statins.txt', sep = '\t', quote = F)

st <- fread('drugs/data/R8_statins.txt')

# Summarise trajectories
stt <- summarizeTrajectories(st, 18)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

stt <- stt %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

fwrite(stt, 'drugs/data/R8_statins_summarized.txt', sep = '\t', quote = F)

# Truncate follow up at 75y
d.lt75 <- fread('drugs/data/R8_statins.txt') %>% 
  filter(EVENT_AGE <= 75)
dd.lt75 <- summarizeTrajectories(d.lt75, 18)

age_first <- getAgeFirstEndpoint(ep, dd.lt75$FINNGENID, ep_chronic)

dd.lt75 <- dd.lt75 %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))
fwrite(dd.lt75, 'drugs/data/R8_statins_summarized_lt75.txt', sep = '\t', quote = F)


# # # # # # # # # # # # # # # # # #  
#   BLOOD PRESSURE MEDICATIONS    #
# # # # # # # # # # # # # # # # # #

# Individual trajectories
bp <- getTrajectoriesDates(purch,'^C0(2|3|8|9)')

fwrite(bp, 'drugs/data/R8_blood_pressure.txt', sep = '\t', quote = F)


# Summarised trajectories
bpp <- summarizeTrajectories(bp, 18)

fwrite(bpp, 'drugs/data/R8_blood_pressure_summarized.txt', sep = '\t', quote = F)


# Truncate follow up at 75y
d.lt75 <- fread('drugs/data/R8_blood_pressure.txt') %>% 
  filter(EVENT_AGE <= 75)
dd.lt75 <- summarizeTrajectories(d.lt75, 18)

fwrite(dd.lt75, 'drugs/data/R8_blood_pressure_summarized_lt75.txt', sep = '\t', quote = F)



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
cd <- getTrajectoriesDates2(purch,ATCs)
fwrite(cd, 'drugs/data/R8_clopi_dipy.txt', sep = '\t', quote = F)

# Summarised trajectories
ccdd <- summarizeTrajectories(cd, 18)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

ccdd <- ccdd %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

fwrite(ccdd, 'drugs/data/R8_clopi_dipy_summarized.txt', sep = '\t', quote = F)


# Truncate follow up at 75y
d.lt75 <- fread('drugs/data/R8_clopi_dipy.txt') %>% 
  filter(EVENT_AGE <= 75)
dd.lt75 <- summarizeTrajectories(d.lt75, 18)

age_first <- getAgeFirstEndpoint(ep, dd.lt75$FINNGENID, ep_chronic)

dd.lt75 <- dd.lt75 %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))
fwrite(dd.lt75, 'drugs/data/R8_clopi_dipy_summarized_lt75.txt', sep = '\t', quote = F)



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
bc <- getTrajectoriesDates(purch,'^L02B(A01|G04|G06|G03)')


fwrite(bc, 'drugs/data/R8_breast_cancer.txt', sep = '\t', quote = F)


# Summarised trajectories
bcc <- summarizeTrajectories(bc)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, bcc$FINNGENID, ep_chronic)

bcc <- bcc %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0)) %>% 
  # exclude non-chronic, re-standardize adherence for GWAS and recalculate age bins
  filter(chronic == 1) %>%
  mutate(adherence_std = as.numeric(scale(adherence)),
         age_bin = cut(age_first_purch, 5))


fwrite(bcc, 'drugs/data/R8_breast_cancer_summarized.txt', sep = '\t', quote = F)


# Truncate follow up at 75y
d.lt75 <- fread('drugs/data/R8_breast_cancer.txt') %>% 
  filter(EVENT_AGE <= 75)
dd.lt75 <- summarizeTrajectories(d.lt75)

age_first <- getAgeFirstEndpoint(ep, dd.lt75$FINNGENID, ep_chronic)

dd.lt75 <- dd.lt75 %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0)) %>% 
  # exclude non-chronic, re-standardize adherence for GWAS and recalculate age bins
  filter(chronic == 1) %>%
  mutate(adherence_std = as.numeric(scale(adherence)),
         age_bin = cut(age_first_purch, 5))


fwrite(dd.lt75, 'drugs/data/R8_breast_cancer_summarized_lt75.txt', sep = '\t', quote = F)


# # # # # # # #
#   DOAC      #
# # # # # # # #

# Define 'chronic' users: events that lead to a strong need of medication:

ep_chronic <- c('I9_DVTANDPULM','I9_VTE','I9_AF')

ATCs <- data.frame(atc = c("B01AF01","B01AF02","B01AF03","B01AE07"),
                   dose = c(1,2,1,2),
                   stringsAsFactors = F)

# Complete trajectories
do <- getTrajectoriesDates2(purch,ATCs)


fwrite(do, 'drugs/data/R8_doac.txt', sep = '\t', quote = F)


# Summarised trajectories
doo <- summarizeTrajectories(do, 18)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, doo$FINNGENID, ep_chronic)

doo <- doo %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

fwrite(doo, 'drugs/data/R8_doac_summarized.txt', sep = '\t', quote = F)

# Truncate follow up at 75y
d.lt75 <- fread('drugs/data/R8_doac.txt') %>% 
  filter(EVENT_AGE <= 75)
dd.lt75 <- summarizeTrajectories(d.lt75, 18)

age_first <- getAgeFirstEndpoint(ep, dd.lt75$FINNGENID, ep_chronic)

dd.lt75 <- dd.lt75 %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

fwrite(dd.lt75, 'drugs/data/R8_doac_summarized_lt75.txt', sep = '\t', quote = F)


# # # # # # # #
#   GLAUCOMA  #
# # # # # # # #

# # # Some pkokoNum need manual check becasue they are wrong
gl <- purch %>%
  mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
  filter(grepl('^S01E(E|D)', CODE1))

pkoko_list <- c("1X2,5 ML", "3 x 3 ml", "3X5ML", "30 x 0.4 g", "90 x 0.3 ml", "30 x 0.2 ml", "90 x 0.2 ml", "90 (18 x 5) x 0.2 ml", "3 x 2.5 ml")

getPkokoNum <- function(pkoko) {
  b <- strsplit(gsub(" |ml|g","",tolower(gsub(",", "\\.", pkoko))), "x")
  return( unlist(lapply(b, function(i) as.numeric(i[1])*as.numeric(i[2]))) )
}

gl$pkoko_num[gl$pkoko %in% pkoko_list] <- getPkokoNum(gl$pkoko[gl$pkoko %in% pkoko_list])

pkoko <- gl %>%
  ungroup() %>%
  select(pkoko, pkoko_num) %>%
  distinct

fwrite(gl, 'drugs/data/R8_glaucoma_all_pkoko_ok.txt', sep = '\t', quote = F)

# Define 'chronic' users: events that lead to a strong need of the medication:

ep_chronic <- c('H7_GLAUCPRIMOPEN')

# Complete trajectories
g <- fread('drugs/data/R8_glaucoma_all_pkoko_ok.txt')

gl <- getTrajectoriesDates(g,'^S01E(E|D)', 0.1)

plot(density(gl$pills_norm, na.rm = T))

getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(gl$pills_norm)

fwrite(gl, 'drugs/data/R8_glauc.txt', sep = '\t', quote = F)

gl <- fread('drugs/data/R8_glauc.txt')

# Summarised trajectories
gll <- summarizeTrajectories(gl, 18)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, gll$FINNGENID, ep_chronic)

gll <- gll %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

fwrite(gll, 'drugs/data/R8_glauc_summarized.txt', sep = '\t', quote = F)


# Truncate follow up at 75y
d.lt75 <- fread('drugs/data/R8_glauc.txt') %>% 
  filter(EVENT_AGE <= 75)
dd.lt75 <- summarizeTrajectories(d.lt75, 18)

age_first <- getAgeFirstEndpoint(ep, dd.lt75$FINNGENID, ep_chronic)

dd.lt75 <- dd.lt75 %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

fwrite(dd.lt75, 'drugs/data/R8_glauc_summarized_lt75.txt', sep = '\t', quote = F)