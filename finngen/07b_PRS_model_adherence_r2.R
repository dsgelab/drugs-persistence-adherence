rm(list=ls())

library("data.table")
library("dplyr")
library("broom")
library(gridExtra)
library(ggplot2)
library(forcats)
library(gtools)

setwd('/home/ivm/')

# #### merge datasets ####
# adh <- fread("/home/ivm/drugs/data/R10_cov_pheno_adherence.txt") %>% mutate(FINNGENID = IID)
# prs <- fread("/home/ivm/drugs/data/traits_for_pgs.csv")
# prs$pheno <- gsub(" ","_", prs$pheno)
# prs$pheno <- gsub("-","_", prs$pheno)
# 
# for (file in prs$loc){
#   score <- fread(file)
#   score$score.std <- scale(score$SCORE1_AVG)
#   name <- paste0(prs$pheno[prs$loc == file],"_score")
#   colnames(score) <- c("#FID",	"IID",	"NMISS_ALLELE_CT",	"NAMED_ALLELE_DOSAGE_SUM",	"SCORE1_AVG", name)
#   adh <- merge(x = adh, y = score[ , c(2,6)], by.x= "FINNGENID", by.y = "IID", all.x=TRUE)
# }
# 
# fwrite(adh,"/home/ivm/drugs/data/R10_cov_pheno_adherence_prs.txt",col.names = T, row.names = F,
#        quote = F, sep="\t", na="NA")
# 
# adhprs <- fread("/home/ivm/drugs/data/R10_cov_pheno_adherence_prs.txt")
# 
# dd <- fread("/home/ivm/drugs/data/R10_blood_pressure_summarized.txt") %>% mutate(blood_pressure_adh = adherence) %>%
#   select(FINNGENID,blood_pressure_adh)
# adhprs <- merge(x = adhprs, y = dd , by = "FINNGENID", all.x=TRUE)
# 
# dd <- fread("/home/ivm/drugs/data/R10_breast_cancer_summarized.txt") %>% mutate(breast_canc_adh = adherence) %>%
#   select(FINNGENID,breast_canc_adh)
# adhprs <- merge(x = adhprs, y = dd , by = "FINNGENID", all.x=TRUE)
# 
# dd <- fread("/home/ivm/drugs/data/R10_clopi_dipy_summarized.txt") %>% mutate(clopi_dipy_adh = adherence) %>%
#   select(FINNGENID,clopi_dipy_adh)
# adhprs <- merge(x = adhprs, y = dd , by = "FINNGENID", all.x=TRUE)
# 
# dd <- fread("/home/ivm/drugs/data/R10_doac_summarized.txt") %>% mutate(doac_adh = adherence) %>%
#   select(FINNGENID,doac_adh)
# adhprs <- merge(x = adhprs, y = dd , by = "FINNGENID", all.x=TRUE)
# 
# dd <- fread("/home/ivm/drugs/data/R10_statins_summarized.txt") %>% mutate(statins_adh = adherence) %>%
#   select(FINNGENID,statins_adh)
# adhprs <- merge(x = adhprs, y = dd , by = "FINNGENID", all.x=TRUE)
# 
# 
# colnames(adhprs)
# 
# 
# fwrite(adhprs,"/home/ivm/drugs/data/R10_cov_pheno_adherence_prs.txt",col.names = T, row.names = F,
#        quote = F, sep="\t", na="NA")

adhprs <- fread("/home/ivm/drugs/data/R10_cov_pheno_adherence_prs.txt")


####  PRS MODEL ####

prs_n <- colnames(adhprs)[grepl("_score", colnames(adhprs))]

pheno <- c("statins",
           "blood_pressure",
           "clopi_dipy",
           "doac",
           "breast_canc")


adhprs <- adhprs %>% mutate(age_first_statins = scale(age_first_statins),                                
                            age_first_blood_pressure = scale(age_first_blood_pressure),                  
                            age_first_clopi_dipy = scale(age_first_clopi_dipy),                         
                            age_first_breast_canc = scale(age_first_breast_canc),                           
                            age_first_doac = scale(age_first_doac),
                            PC1 = scale(PC1),
                            PC2 = scale(PC2),
                            PC3 = scale(PC3),
                            PC4 = scale(PC4),
                            PC5 = scale(PC5),
                            PC6 = scale(PC6),
                            PC7 = scale(PC7),
                            PC8 = scale(PC8),
                            PC9 = scale(PC9),
                            PC10 = scale(PC10))

param <- data.frame()
#med <- "statins"

res <- NULL
for (med in pheno) {
  
  betas <- data.frame()
  age <- paste0("age_first_",med)
  medi <- paste0(med,"_adh")

  #i <- "Participation_ukb_FFQ_score"
  for (i in prs_n){
    print(i)
    
    # Baseline R2, with and withouth PCs
    if (med == "breast_canc"){
      base <- paste0(medi," ~ ", age)
      base_pc <- paste0(base, " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    } else {
      base <- paste0(medi," ~ factor(SEX_IMPUTED) + ", age)
      base_pc <- paste0(base, " + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    }
    
    model_base <- lm(formula = base, data = adhprs)
    model_base_pc <- lm(formula = base_pc, data = adhprs)
    base_r2 <- summary(model_base)$adj.r.squared
    base_pc_r2 <- summary(model_base_pc)$adj.r.squared
    
    # PGS r2, with and without PCs
    if (med == "breast_canc"){
      subdata <-adhprs %>% select(all_of(medi), all_of(age),
                                  "PC1", "PC2", "PC3", "PC4", "PC5",
                                  "PC6", "PC7", "PC8", "PC9", "PC10",
                                  all_of(i))
      model_pc <- paste0(medi," ~ .")
    } else {
      subdata <-adhprs %>% select(all_of(medi), all_of(age), "SEX_IMPUTED",
                                "PC1", "PC2", "PC3", "PC4", "PC5",
                                "PC6", "PC7", "PC8", "PC9", "PC10",
                                all_of(i))
      model_pc <- paste0(medi," ~ factor(SEX_IMPUTED) + .")
    }
    
    m_pcs <- lm(model_pc, subdata)
    pgs_pc_r2 <- summary(m_pcs)$adj.r.squared
    
    if (med == "breast_canc"){
      subdata <-adhprs %>% select(all_of(medi), all_of(age), all_of(i))
      model <- paste0(medi," ~ .")
    } else {
      subdata <-adhprs %>% select(all_of(medi), all_of(age), "SEX_IMPUTED", all_of(i))
      model <- paste0(medi," ~ factor(SEX_IMPUTED) + .")
    }
    
    m <- lm(model, subdata)
    pgs_r2 <- summary(m)$adj.r.squared
    
    res_i <- data.frame(drug = med, predictor = i, base_r2 = base_r2, pgs_r2 = pgs_r2, base_pc_r2 = base_pc_r2, pgs_pc_r2 = pgs_pc_r2)
    res <- bind_rows(res, res_i)
    

  }
}

fwrite(res,"/home/ivm/drugs/results/manuscript/PGS_adh_r2_2.tsv",sep = "\t",
       col.names = T,row.names = F, quote = F)