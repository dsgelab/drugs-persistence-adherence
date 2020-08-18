rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

source('/home/cordioli/drugs/adherence_funs.R')

write.flag <- F

purch <- fread('/home/cordioli/drugs/data/R5_v3_purch_vnr_98.gz')
ep <- fread('/home/cordioli/R5_pheno/finngen_R5_v3_endpoint.gz')
covs <- fread('/home/cordioli/drugs/data/R5_cov.txt')


# # # # # # # #
#   STATINS   #
# # # # # # # #

# Define 'chronic' users: events that lead to a strong need of statins:
# I9_ASO
# I9_CHD
# I9_ATHSCLE
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
st <- getTrajectories(purch,'^C10AA')

# Summarised trajectories
stt <- summarizeTrajectories(st)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

stt <- stt %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

stt <- stt %>%
  left_join(covs)

# Difference between sexes
t.test(stt$adherence[stt$SEX_IMPUTED==0], stt$adherence[stt$SEX_IMPUTED==1], alternative = 'greater')$p.value

# Primary VS secondary prevention
t.test(stt$adherence[stt$chronic==0], stt$adherence[stt$chronic==1], paired = F)

# Regularity
cor(stt$adherence,stt$sd_days)
cor.test(stt$adherence,stt$sd_days)


if (write.flag == T) {
  fwrite(stt, '/home/cordioli/drugs/data/statins_summarized.txt', sep = '\t', quote = F)
}

# # # 
# Merge phenotypes with covariates for GWAS
ad <- stt %>%
  mutate(statins = adherence_std,
         age_first_statins = age_first_purch) %>%
  select(FINNGENID,statins,age_first_statins)

covs <- covs %>%
  left_join(ad)

if (write.flag == T) {
  fwrite(covs, '/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt', sep = '\t', quote = F)
}



# # # # # # # # # # # # # # # # # #  
#   BLOOD PRESSURE MEDICATIONS    #
# # # # # # # # # # # # # # # # # #

# Define 'chronic' users: events that lead to a strong need of bp:

# TODO: add endpoints
# ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
#                 'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Individual trajectories
bp <- getTrajectories(purch,'^C0(2|3|8|9)')

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
  fwrite(bpp, '/home/cordioli/drugs/data/bp_summarized.txt', sep = '\t', quote = F)
}

# # # 
# Merge phenotypes with covariates for GWAS
ad <- bpp %>%
  mutate(blood_pressure = adherence_std,
         age_first_blood_pressure = age_first_purch) %>%
  select(FINNGENID,blood_pressure,age_first_blood_pressure)

covs <- fread('/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt')

covs <- covs %>%
  left_join(ad)

if (write.flag == T) {
  fwrite(covs, '/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt', sep = '\t', quote = F)
}

# # # # # # # # # # # # # # # # # #  
#   CLOPIDOGREL & ASA+DIP   #
# # # # # # # # # # # # # # # # # #
# clopi '^B01AC04'
# dipydogrol '^B01AC30' -> adjust pills *0.5 (dose is 2pills/day)

# Define 'chronic' users: events that lead to a strong need of bp:

# TODO: add endpoints
# ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
#                 'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Individual trajectories
clopi <- getTrajectories(purch,'^B01AC04')
dipy <- getTrajectories(purch,'^B01AC30')

# TODO: refactor this in the function
# Dipy: adjust pills *0.5 (dose is 2pills/day)
dipy$n_pills <- dipy$n_pills*0.5

# Summarised trajectories
cl <- summarizeTrajectories(clopi)
di <- summarizeTrajectories(dipy)

cldi <- cl %>%
  bind_rows(di) %>%
  distinct(FINNGENID, .keep_all = T)

# Age first related event
# age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

# stt <- stt %>%
#   left_join(age_first) %>%
#   mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
#                              TRUE ~ 0))
if (write.flag == T) {
  fwrite(cldi, '/home/cordioli/drugs/data/clopi_dipy_summarized.txt', sep = '\t', quote = F)
}

# # # 
# Merge phenotypes with covariates for GWAS
ad <- cldi %>%
  mutate(clopi_dipy = as.numeric(scale(adherence)),
         age_first_clopi_dipy = age_first_purch) %>%
  select(FINNGENID,clopi_dipy,age_first_clopi_dipy)

covs <- fread('/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt')

covs <- covs %>%
  left_join(ad)
if (write.flag == T) {
  fwrite(covs, '/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt', sep = '\t', quote = F)
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

# Summarised trajectories
bcc <- summarizeTrajectories(bc)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, bcc$FINNGENID, ep_chronic)

bcc <- bcc %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))


if (write.flag == T) {
  fwrite(bcc, '/home/cordioli/drugs/data/breast_cancer_summarized.txt', sep = '\t', quote = F)
}


# # # 
# Exclude non-chronic and te-standardize adherence for GWAS
bcc <- bcc %>%
  filter(chronic == 1) %>%
  mutate(adherence_std = as.numeric(scale(adherence)))

# Merge phenotypes with covariates for GWAS
covs <- fread('/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt')

ad <- bcc %>%
  mutate(breast_canc = adherence_std,
         age_first_breast_canc = age_first_purch) %>%
  select(FINNGENID,breast_canc,age_first_breast_canc)

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt', sep = '\t', quote = F)



# # # # # # # # # # #
#   EPILEPSY        #
# # # # # # # # # # #

# 

# m <- fread('/home/cordioli/drugs/data/epilepsy_meds.txt')
# m <- m[m$perc_epilepsy>=40,]
# sort(m$drug_code)

# ep_chronic <- c('G6_EPLEPSY')
# 
# # Complete trajectories
# epi <- getTrajectories(purch,'^N03A|B02|B52|D01|F02|F03|F04|G04|G06|X14|X15|X17|X18|X22|X23')
# 
# # Summarised trajectories
# epii <- summarizeTrajectories(df, 10)
# 
# # Age first related event
# age_first <- getAgeFirstEndpoint(ep, epii$FINNGENID, ep_chronic)
# 
# epii <- epii %>%
#   left_join(age_first) %>%
#   mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
#                              TRUE ~ 0))
# # if (write.flag == T) {
# #   fwrite(epii, '/home/cordioli/drugs/data/breast_cancer_summarized.txt', sep = '\t', quote = F)
# # }
# 


phenolist <- c('statins','blood_pressure','clopi_dipy','breast_canc')
fwrite(list(phenolist), '/home/cordioli/drugs/data/adherence_phenolist.txt', col.names = F)

system('gsutil cp /home/cordioli/drugs/data/R5_cov_pheno_adherence.txt gs://mattia/drugs/')