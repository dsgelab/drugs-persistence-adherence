rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

covs <- fread('/home/cordioli/drugs/data/R5_cov.txt')

st <- fread('drugs/data/statins_summarized.txt')
st <- merge(st, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

bp <- fread('drugs/data/bp_summarized.txt')
bp <- merge(bp, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

cd <- fread('drugs/data/clopi_dipy_summarized.txt')
cd <- merge(cd, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

br <- fread('drugs/data/breast_cancer_summarized.txt')
br <- merge(br, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

# Mean adherence for each drug
res <- data.frame(drug = c("statins", "bp_meds", "clopidogrel_dipyridamole", "breast_cancer_meds"))

res$N <- c(nrow(st), nrow(bp), nrow(cd), nrow(br))

res$N_secondary_prev <- c(sum(st$chronic), NA, NA, sum(br$chronic))

res$mean_age_first <- c(mean(st$age_first_purch), mean(bp$age_first_purch), mean(cd$age_first_purch), mean(br$age_first_purch))

res$mean_adherence <- c(mean(st$adherence), mean(bp$adherence), mean(cd$adherence), mean(br$adherence))

res$sd_adherence <- c(sd(st$adherence), sd(bp$adherence), sd(cd$adherence), sd(br$adherence))

res$adherence_gt_0.8 <- c(length(which(st$adherence > 0.8))/nrow(st),
                          length(which(bp$adherence > 0.8))/nrow(bp),
                          length(which(cd$adherence > 0.8))/nrow(cd),
                          length(which(br$adherence > 0.8))/nrow(br))

# 44 males in breast cancer!
res$mean_adherence_M <- c(mean(st$adherence[st$SEX == "male"]), mean(bp$adherence[bp$SEX == "male"]), mean(cd$adherence[cd$SEX == "male"]), mean(br$adherence[br$SEX == "male"]))
res$mean_adherence_F <- c(mean(st$adherence[st$SEX == "female"]), mean(bp$adherence[bp$SEX == "female"]), mean(cd$adherence[cd$SEX == "female"]), mean(br$adherence[br$SEX == "female"]))

res$mean_adherence_1st_prev <- c(mean(st$adherence[st$chronic == 0]), NA, NA, mean(br$adherence[br$chronic == 0]))
res$mean_adherence_2nd_prev <- c(mean(st$adherence[st$chronic == 1]), NA, NA, mean(br$adherence[br$chronic == 1]))

# correlation trajectory length
res$mean_treatment_years <- c(mean(st$tot_days), mean(bp$tot_days), mean(cd$tot_days), mean(br$tot_days))
res$mean_treatment_years <- res$mean_treatment_years/365.25
res$cor_length <- c(round(cor(st$adherence, st$tot_days), 4), round(cor(bp$adherence, bp$tot_days), 4), round(cor(cd$adherence, cd$tot_days), 4), round(cor(br$adherence, br$tot_days), 4))

# correlation regularity
res$mean_sd_days <- c(mean(st$sd_days_norm), mean(bp$sd_days_norm), mean(cd$sd_days_norm), mean(br$sd_days_norm))
res$cor_sd_days <- c(round(cor(st$adherence, st$sd_days_norm), 4), round(cor(bp$adherence, bp$sd_days_norm), 4), round(cor(cd$adherence, cd$sd_days_norm), 4), round(cor(br$adherence, br$sd_days_norm), 4))

resT <- data.frame(t(res)[-1,], stringsAsFactors = F)
colnames(resT) <- t(res)[1,]