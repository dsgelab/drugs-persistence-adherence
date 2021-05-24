rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

covs <- fread('/home/cordioli/drugs/data/finngen_R6_cov.txt')

st <- fread('drugs/data/R6_statins_summarized_150.txt')
#st <- merge(st, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

bp <- fread('drugs/data/R6_bp_summarized_150.txt')
#bp <- merge(bp, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

cd <- fread('drugs/data/R6_clopi_dipy_summarized_150.txt')
#cd <- merge(cd, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

br <- fread('drugs/data/R6_breast_cancer_summarized_150.txt')
br <- merge(br, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

as <- fread('drugs/data/R6_asthma_summarized.txt')
#as <- merge(as, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

do <- fread('drugs/data/R6_doac_summarized_150.txt')
#do <- merge(do, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")

gl <- fread('drugs/data/R6_glauc_summarized.txt')
#gl <- merge(gl, covs[, c("FINNGENID", "SEX")], by = "FINNGENID")


# Mean adherence for each drug
res <- data.frame(drug = c("statins", "bp_meds", "clopidogrel_dipyridamole", "breast_cancer_meds", "asthma_meds", "doac", "glaucoma_meds"))

res$N <- c(nrow(st), nrow(bp), nrow(cd), nrow(br), nrow(as), nrow(do), nrow(gl))

res$N_secondary_prev <- c(sum(st$chronic), NA, sum(cd$chronic), sum(br$chronic), sum(as$chronic), sum(do$chronic), sum(gl$chronic))

res$tot_purch <- c(sum(st$tot_purch), sum(bp$tot_purch), sum(cd$tot_purch), sum(br$tot_purch), sum(as$tot_purch), sum(do$tot_purch), sum(gl$tot_purch))

res$mean_age_first <- c(mean(st$age_first_purch), mean(bp$age_first_purch), mean(cd$age_first_purch), mean(br$age_first_purch), mean(as$age_first_purch), mean(do$age_first_purch), mean(gl$age_first_purch))

res$mean_adherence <- c(mean(st$adherence), mean(bp$adherence), mean(cd$adherence), mean(br$adherence), mean(as$adherence), mean(do$adherence), mean(gl$adherence))

res$sd_adherence <- c(sd(st$adherence), sd(bp$adherence), sd(cd$adherence), sd(br$adherence), sd(as$adherence), sd(do$adherence), sd(gl$adherence))

res$adherence_gt_0.8 <- c(length(which(st$adherence > 0.8))/nrow(st),
                          length(which(bp$adherence > 0.8))/nrow(bp),
                          length(which(cd$adherence > 0.8))/nrow(cd),
                          length(which(br$adherence > 0.8))/nrow(br),
                          length(which(as$adherence > 0.8))/nrow(as),
                          length(which(do$adherence > 0.8))/nrow(do),
                          length(which(gl$adherence > 0.8))/nrow(gl))

# 44 males in breast cancer!
res$mean_adherence_M <- c(mean(st$adherence[st$SEX == "male"]), mean(bp$adherence[bp$SEX == "male"]), mean(cd$adherence[cd$SEX == "male"]), mean(br$adherence[br$SEX == "male"]),
                          mean(as$adherence[as$SEX == "male"]), mean(do$adherence[do$SEX == "male"]), mean(gl$adherence[gl$SEX == "male"]))

res$mean_adherence_F <- c(mean(st$adherence[st$SEX == "female"]), mean(bp$adherence[bp$SEX == "female"]), mean(cd$adherence[cd$SEX == "female"]), mean(br$adherence[br$SEX == "female"]),
                          mean(as$adherence[as$SEX == "female"]), mean(do$adherence[do$SEX == "female"]), mean(gl$adherence[gl$SEX == "female"]))

res$sd_adherence_M <- c(sd(st$adherence[st$SEX == "male"]), sd(bp$adherence[bp$SEX == "male"]), sd(cd$adherence[cd$SEX == "male"]), sd(br$adherence[br$SEX == "male"]),
                        sd(as$adherence[as$SEX == "male"]), sd(do$adherence[do$SEX == "male"]), sd(gl$adherence[gl$SEX == "male"]))

res$sd_adherence_F <- c(sd(st$adherence[st$SEX == "female"]), sd(bp$adherence[bp$SEX == "female"]), sd(cd$adherence[cd$SEX == "female"]), sd(br$adherence[br$SEX == "female"]),
                        sd(as$adherence[as$SEX == "female"]), sd(do$adherence[do$SEX == "female"]), sd(gl$adherence[gl$SEX == "female"]))

res$mean_adherence_1st_prev <- c(mean(st$adherence[st$chronic == 0]), NA, mean(cd$adherence[cd$chronic == 0]), mean(br$adherence[br$chronic == 0]),
                                 mean(as$adherence[as$chronic == 0]), mean(do$adherence[do$chronic == 0]), mean(gl$adherence[gl$chronic == 0]))

res$mean_adherence_2nd_prev <- c(mean(st$adherence[st$chronic == 1]), NA, mean(cd$adherence[cd$chronic == 1]), mean(br$adherence[br$chronic == 1]),
                                 mean(as$adherence[as$chronic == 1]), mean(do$adherence[do$chronic == 1]), mean(gl$adherence[gl$chronic == 1]))

res$sd_adherence_1st_prev <- c(sd(st$adherence[st$chronic == 0]), NA, sd(cd$adherence[cd$chronic == 0]), sd(br$adherence[br$chronic == 0]),
                               sd(as$adherence[as$chronic == 0]), sd(do$adherence[do$chronic == 0]), sd(gl$adherence[gl$chronic == 0]))

res$sd_adherence_2nd_prev <- c(sd(st$adherence[st$chronic == 1]), NA, sd(cd$adherence[cd$chronic == 1]), sd(br$adherence[br$chronic == 1]),
                               sd(as$adherence[as$chronic == 1]), sd(do$adherence[do$chronic == 1]), sd(gl$adherence[gl$chronic == 1]))


# correlation trajectory length
res$mean_treatment_years <- c(mean(st$tot_days), mean(bp$tot_days), mean(cd$tot_days), mean(br$tot_days), mean(as$tot_days), mean(do$tot_days), mean(gl$tot_days))
res$mean_treatment_years <- res$mean_treatment_years/365.25
res$cor_length <- c(round(cor(st$adherence, st$tot_days), 4), round(cor(bp$adherence, bp$tot_days), 4), round(cor(cd$adherence, cd$tot_days), 4), round(cor(br$adherence, br$tot_days), 4),
                    round(cor(as$adherence, as$tot_days), 4), round(cor(do$adherence, do$tot_days), 4), round(cor(gl$adherence, gl$tot_days), 4))


# correlation regularity
res$mean_sd_days <- c(mean(st$sd_days_norm), mean(bp$sd_days_norm), mean(cd$sd_days_norm), mean(br$sd_days_norm),
                      mean(as$sd_days_norm), mean(do$sd_days_norm), mean(gl$sd_days_norm))
res$cor_sd_days <- c(round(cor(st$adherence, st$sd_days_norm), 4), round(cor(bp$adherence, bp$sd_days_norm), 4), round(cor(cd$adherence, cd$sd_days_norm), 4), round(cor(br$adherence, br$sd_days_norm), 4),
                     round(cor(as$adherence, as$sd_days_norm), 4), round(cor(do$adherence, do$sd_days_norm), 4), round(cor(gl$adherence, gl$sd_days_norm), 4))



resT <- data.frame(t(res)[-1,], stringsAsFactors = F)
colnames(resT) <- t(res)[1,]