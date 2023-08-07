rm(list = ls())

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
'%notlike%'= Negate('%like%')
'%notin%'= Negate('%in%')

setwd("~")

adhprs <- fread("/home/ivm/drugs/data/R10_cov_pheno_adherence_prs.txt")

drug_names <- data.frame(Medication = c("blood_pressure", "breast_canc", "clopi_dipy", "doac", "statins"),
                         name = c("Blood pressure medications", "Breast cancer medications", "Antiplatelets", "Anticoagulants", "Statins"))

d <- adhprs %>% 
  select(all_of(ends_with("_adh")), Breast_cancer_score, Systolic_blood_pressure_score, Low_density_lipoproteins_score, Low_density_lipoproteins_adjusted_score, Stroke_score, all_of(starts_with("Participation_"))) %>% 
  pivot_longer(cols = c("blood_pressure_adh", "breast_canc_adh", "clopi_dipy_adh", "doac_adh", "statins_adh"),
               names_to = "Medication",
               values_to = "Adherence") %>% 
  mutate(Medication = gsub("_adh", "", Medication)) %>% 
  left_join(drug_names) %>% 
  select(-Medication) %>% 
  mutate(Medication = factor(name)) %>% 
  select(-name) %>% 
  mutate(across(Breast_cancer_score:Participation_ukb_WRB_score, ~ cut(., breaks = quantile(., probs = seq(0, 1, by = 0.1)), include.lowest = TRUE, labels = FALSE))) %>% 
  pivot_longer(cols = Breast_cancer_score:Participation_ukb_WRB_score,
               names_to = "Trait",
               values_to = "Decile") %>% 
  mutate(Trait = gsub("_score", "", Trait)) %>% 
  filter(!is.na(Adherence)) %>% 
  group_by(Medication, Trait, Decile) %>% 
  summarise(Mean_adh = mean(Adherence), SD_adh = sd(Adherence))

meds <- c("Blood pressure medications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer medications")
d$Medication <- factor(d$Medication, levels = meds)
d$Trait <- factor(d$Trait, levels = sort(unique(d$Trait)))
palette <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00")

         
png("/home/ivm/drugs/results/manuscript/pgs/plots_pgs-dec_adherence", width=2000, height=1500, res = 300)
ggplot(d, aes(x = Decile, y = Mean_adh, color = Medication)) +
  facet_wrap(Medication~Trait, ncol = 5) +
  geom_pointrange(aes(ymin = Mean_adh - SD_adh, ymax = Mean_adh + SD_adh),
                  position = position_dodge(width = 0.5)) +
  theme_minimal() +
  scale_color_manual(values = palette)
dev.off()


d2 <- adhprs %>% 
  select(all_of(ends_with("_adh")), Breast_cancer_score, Systolic_blood_pressure_score, Low_density_lipoproteins_score, Low_density_lipoproteins_adjusted_score, Stroke_score, all_of(starts_with("Participation_"))) %>% 
  pivot_longer(cols = c("blood_pressure_adh", "breast_canc_adh", "clopi_dipy_adh", "doac_adh", "statins_adh"),
               names_to = "Medication",
               values_to = "Adherence") %>% 
  mutate(Medication = gsub("_adh", "", Medication)) %>% 
  left_join(drug_names) %>% 
  select(-Medication) %>% 
  mutate(Medication = factor(name)) %>% 
  select(-name) %>% 
  mutate(across(Breast_cancer_score:Participation_ukb_WRB_score, ~ as.numeric(scale(.)))) %>% 
  pivot_longer(cols = Breast_cancer_score:Participation_ukb_WRB_score,
               names_to = "Trait",
               values_to = "Scaled_score") %>% 
  mutate(Trait = gsub("_score", "", Trait)) %>% 
  filter(!is.na(Adherence)) %>%
  mutate(Good_adh = ifelse(Adherence > 0.8, 1, 0)) %>% 
  group_by(Medication, Trait, Good_adh) %>% 
  summarise(Mean_score = mean(Scaled_score), SD_score = sd(Scaled_score))


meds <- c("Blood pressure medications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer medications")
d2$Medication <- factor(d2$Medication, levels = meds)
d2$Trait <- factor(d2$Trait, levels = sort(unique(d2$Trait)))
palette <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00")


png("/home/ivm/drugs/results/manuscript/pgs/plots_pgs-dec_adherence", width=2000, height=1500, res = 300)
ggplot(d2, aes(x = Good_adh, y = Mean_score, color = Medication)) +
  facet_wrap(Medication~Trait, ncol = 5) +
  geom_pointrange(aes(ymin = Mean_score - SD_score, ymax = Mean_score + SD_score),
                  position = position_dodge(width = 0.5)) +
  theme_minimal() +
  scale_color_manual(values = palette)
dev.off()

