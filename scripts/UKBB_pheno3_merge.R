#### SET ####
library(BBmisc)
library(stringr)
library(data.table)
library(dplyr)


d <- fread("~/data/doac_pheno_final.tsv")
a <- fread("~/data/antiplat_pheno_final.tsv")
bc <- fread("~/data/breastcanc_pheno_final.tsv")
bp <- fread("~/data/bloodpres_pheno_final.tsv")
g <- fread("~/data/glaucoma_pheno_final.tsv")
s <- fread("~/data/statins_pheno_final.tsv")

View(table(d$pop))
View(table(a$pop))
View(table(bc$pop))
View(table(bp$pop))
View(table(g$pop))
View(table(s$pop))


pops <- unique(s$pop)


#### DOAC ####
eur <- d %>%
  filter(pop=="EUR") %>%
  mutate(doac_eur = as.numeric(scale(adherence)),
          age_first_doac_eur = age_first_purch)

afr <- d %>%
  filter(pop=="AFR") %>%
  mutate(doac_afr = as.numeric(scale(adherence)),
         age_first_doac_afr = age_first_purch)

csa <- d %>%
  filter(pop=="CSA") %>%
  mutate(doac_csa = as.numeric(scale(adherence)),
         age_first_doac_csa = age_first_purch)

mid <- d %>%
  filter(pop=="MID") %>%
  mutate(doac_mid = as.numeric(scale(adherence)),
         age_first_doac_mid = age_first_purch)

amr <- d %>%
  filter(pop=="AMR") %>%
  mutate(doac_amr= as.numeric(scale(adherence)),
         age_first_doac_amr = age_first_purch)

eas <- d %>%
  filter(pop=="EAS") %>%
  mutate(doac_eas = as.numeric(scale(adherence)),
         age_first_doac_eas = age_first_purch)


doac_all <- merge(merge(merge(merge(merge(
  eur, 
  afr, all = TRUE),
  csa, all = TRUE),
  mid, all = TRUE),
  amr, all = TRUE),
  eas, all = TRUE)

#### ANTIPLATELET ####

eur <- a %>%
  filter(pop=="EUR") %>%
  mutate(antiplat_eur = as.numeric(scale(adherence)),
         age_first_antiplat_eur = age_first_purch)

afr <- a %>%
  filter(pop=="AFR") %>%
  mutate(antiplat_afr = as.numeric(scale(adherence)),
         age_first_antiplat_afr = age_first_purch)

csa <- a %>%
  filter(pop=="CSA") %>%
  mutate(antiplat_csa = as.numeric(scale(adherence)),
         age_first_antiplat_csa = age_first_purch)

mid <- a %>%
  filter(pop=="MID") %>%
  mutate(antiplat_mid = as.numeric(scale(adherence)),
         age_first_antiplat_mid = age_first_purch)

amr <- a %>%
  filter(pop=="AMR") %>%
  mutate(antiplat_amr= as.numeric(scale(adherence)),
         age_first_antiplat_amr = age_first_purch)

eas <- a %>%
  filter(pop=="EAS") %>%
  mutate(antiplat_eas = as.numeric(scale(adherence)),
         age_first_antiplat_eas = age_first_purch)


ap_all <- merge(merge(merge(merge(merge(
  eur, 
  afr, all = TRUE),
  csa, all = TRUE),
  mid, all = TRUE),
  amr, all = TRUE),
  eas, all = TRUE)





#### BREAST CANCER ####

eur <- bc %>%
  filter(pop=="EUR") %>%
  mutate(breastcanc_eur = as.numeric(scale(adherence)),
         age_first_breastcanc_eur = age_first_purch)

afr <- bc %>%
  filter(pop=="AFR") %>%
  mutate(breastcanc_afr = as.numeric(scale(adherence)),
         age_first_breastcanc_afr = age_first_purch)

csa <- bc %>%
  filter(pop=="CSA") %>%
  mutate(breastcanc_csa = as.numeric(scale(adherence)),
         age_first_breastcanc_csa = age_first_purch)

mid <- bc %>%
  filter(pop=="MID") %>%
  mutate(breastcanc_mid = as.numeric(scale(adherence)),
         age_first_breastcanc_mid = age_first_purch)

amr <- bc %>%
  filter(pop=="AMR") %>%
  mutate(breastcanc_amr= as.numeric(scale(adherence)),
         age_first_breastcanc_amr = age_first_purch)

eas <- bc %>%
  filter(pop=="EAS") %>%
  mutate(breastcanc_eas = as.numeric(scale(adherence)),
         age_first_breastcanc_eas = age_first_purch)


bc_all <- merge(merge(merge(merge(merge(
  eur, 
  afr, all = TRUE),
  csa, all = TRUE),
  mid, all = TRUE),
  amr, all = TRUE),
  eas, all = TRUE)

#### BLOOD PRESSURE ####

eur <- bp %>%
  filter(pop=="EUR") %>%
  mutate(bloodpres_eur = as.numeric(scale(adherence)),
         age_first_bloodpres_eur = age_first_purch)

afr <- bp %>%
  filter(pop=="AFR") %>%
  mutate(bloodpres_afr = as.numeric(scale(adherence)),
         age_first_bloodpres_afr = age_first_purch)

csa <- bp %>%
  filter(pop=="CSA") %>%
  mutate(bloodpres_csa = as.numeric(scale(adherence)),
         age_first_bloodpres_csa = age_first_purch)

mid <- bp %>%
  filter(pop=="MID") %>%
  mutate(bloodpres_mid = as.numeric(scale(adherence)),
         age_first_bloodpres_mid = age_first_purch)

amr <- bp %>%
  filter(pop=="AMR") %>%
  mutate(bloodpres_amr= as.numeric(scale(adherence)),
         age_first_bloodpres_amr = age_first_purch)

eas <- bp %>%
  filter(pop=="EAS") %>%
  mutate(bloodpres_eas = as.numeric(scale(adherence)),
         age_first_bloodpres_eas = age_first_purch)


bp_all <- merge(merge(merge(merge(merge(
  eur, 
  afr, all = TRUE),
  csa, all = TRUE),
  mid, all = TRUE),
  amr, all = TRUE),
  eas, all = TRUE)


#### GLAUCOMA ####
eur <- g %>%
  filter(pop=="EUR") %>%
  mutate(glaucoma_eur = as.numeric(scale(adherence)),
         age_first_glaucoma_eur = age_first_purch)

afr <- g %>%
  filter(pop=="AFR") %>%
  mutate(glaucoma_afr = as.numeric(scale(adherence)),
         age_first_glaucoma_afr = age_first_purch)

csa <- g %>%
  filter(pop=="CSA") %>%
  mutate(glaucoma_csa = as.numeric(scale(adherence)),
         age_first_glaucoma_csa = age_first_purch)

mid <- g %>%
  filter(pop=="MID") %>%
  mutate(glaucoma_mid = as.numeric(scale(adherence)),
         age_first_glaucoma_mid = age_first_purch)

amr <- g %>%
  filter(pop=="AMR") %>%
  mutate(glaucoma_amr= as.numeric(scale(adherence)),
         age_first_glaucoma_amr = age_first_purch)

eas <- g %>%
  filter(pop=="EAS") %>%
  mutate(glaucoma_eas = as.numeric(scale(adherence)),
         age_first_glaucoma_eas = age_first_purch)


g_all <- merge(merge(merge(merge(merge(
  eur, 
  afr, all = TRUE),
  csa, all = TRUE),
  mid, all = TRUE),
  amr, all = TRUE),
  eas, all = TRUE)




#### STATINS ####
eur <- s %>%
  filter(pop=="EUR") %>%
  mutate(statins_eur = as.numeric(scale(adherence)),
         age_first_statins_eur = age_first_purch)

afr <- s %>%
  filter(pop=="AFR") %>%
  mutate(statins_afr = as.numeric(scale(adherence)),
         age_first_statins_afr = age_first_purch)

csa <- s %>%
  filter(pop=="CSA") %>%
  mutate(statins_csa = as.numeric(scale(adherence)),
         age_first_statins_csa = age_first_purch)

mid <- s %>%
  filter(pop=="MID") %>%
  mutate(statins_mid = as.numeric(scale(adherence)),
         age_first_statins_mid = age_first_purch)

amr <- s %>%
  filter(pop=="AMR") %>%
  mutate(statins_amr= as.numeric(scale(adherence)),
         age_first_statins_amr = age_first_purch)

eas <- s %>%
  filter(pop=="EAS") %>%
  mutate(statins_eas = as.numeric(scale(adherence)),
         age_first_statins_eas = age_first_purch)


s_all <- merge(merge(merge(merge(merge(
  eur, 
  afr, all = TRUE),
  csa, all = TRUE),
  mid, all = TRUE),
  amr, all = TRUE),
  eas, all = TRUE)

#### MERGE ####
TUTTI <- merge(merge(merge(merge(merge(
  ap_all, 
  bc_all, all = TRUE),
  bp_all, all = TRUE),
  doac_all, all = TRUE),
  g_all, all = TRUE),
  s_all, all = TRUE)




fwrite(TUTTI,"~/data/ukbb_pheno_TUTTI.tsv")

colnames(TUTTI)
