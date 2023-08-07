rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(forcats)
library(cowplot)

traits <- c("ADHD",
            "Age smoking initiation",
            "Alzheimer's disease",
            "Asthma",
            "Bipolar disorder",
            "BMI",
            "Breast cancer",
            "Chronic pain",
            "Cigarettes per day",
            "Cognitive performance",
            "Coronary artery disease",
            "Depression",
            "Diabetes",
            "Diastolic blood pressure",
            "Drinks per week",
            "Educational attainment",
            "Glycated hemoglobin",
            "Insomnia",
            "Lifespan",
            "Low-density lipoproteins",
            "Low-density lipoproteins adjusted",
            "Neuroticism",
            "Prostate cancer",
            "Rheumatoid arthritis",            
            "Risk tolerance",            
            "Schizophrenia",
            "Stroke",
            "Systolic blood pressure",
            "Triglycerides",
            "Participation in UKB aide-memoire",
            "Participation in UKB food questionnaire",
            "Completion of UKB mental health questionnaire",
            "Participation in UKB wearable")

drugs <- c("Statins",
           "Blood pressure\nmedications",
           "Antiplatelets",
           "Anticoagulants",
           "Breast cancer\nmedications")

alpha <- 0.05
multiplier <- qnorm(1 - alpha / 2)

# Read in adherence / stop file
pgs <- fread("/Users/cordioli/Projects/drugs/manuscript/pgs/pgs_results/R10_betas_prs_adh.tsv")

traits_mapping <- data.frame(term = unique(pgs$term),
                             trait = traits)

drugs_mapping <- data.frame(medicament = unique(pgs$medicament),
                            drug = drugs)

pgs <- pgs %>% 
  left_join(traits_mapping) %>% 
  left_join(drugs_mapping) %>% 
  mutate(lower = estimate - (multiplier * std.error),
         upper = estimate + (multiplier * std.error),
         perc_change = estimate/1.1,
         perc_lower = lower/1.1,
         perc_upper = upper/1.1)

#### Assign class to each trait #### 
# Psychiatric
pgs$class[pgs$trait=="Depression"] <- "Psychiatric"
pgs$class[pgs$trait=="Neuroticism"] <- "Psychiatric"
pgs$class[pgs$trait=="Schizophrenia"] <- "Psychiatric"
pgs$class[pgs$trait=="Bipolar disorder"] <- "Psychiatric"
pgs$class[pgs$trait=="ADHD"] <- "Psychiatric"

# Biomarker
pgs$class[pgs$trait=="BMI"] <- "Biomarker"
pgs$class[pgs$trait=="Systolic blood pressure"] <- "Biomarker"
pgs$class[pgs$trait=="Diastolic blood pressure"] <- "Biomarker"
pgs$class[pgs$trait=="Glycated hemoglobin"] <- "Biomarker"
pgs$class[pgs$trait=="Triglycerides"] <- "Biomarker"
pgs$class[pgs$trait=="Low-density lipoproteins"] <- "Biomarker"
pgs$class[pgs$trait=="Low-density lipoproteins adjusted"] <- "Biomarker"

# Disease liability
pgs$class[pgs$trait=="Coronary artery disease"] <- "Disease liability"
pgs$class[pgs$trait=="Stroke"] <- "Disease liability"
pgs$class[pgs$trait=="Atrial Fibrillation"] <- "Disease liability"
pgs$class[pgs$trait=="Alzheimer's disease"] <- "Disease liability"
pgs$class[pgs$trait=="Rheumatoid arthritis"] <- "Disease liability"
pgs$class[pgs$trait=="Asthma"] <- "Disease liability"
pgs$class[pgs$trait=="Prostate cancer"] <- "Disease liability"
pgs$class[pgs$trait=="Breast cancer"] <- "Disease liability"
pgs$class[pgs$trait=="Irritable bowel syndrome"] <- "Disease liability"
pgs$class[pgs$trait=="Glaucoma"] <- "Disease liability"
pgs$class[pgs$trait=="Diabetes"] <- "Disease liability"
pgs$class[pgs$trait=="Chronic pain"] <- "Disease liability"
pgs$class[pgs$trait=="Lifespan"] <- "Disease liability"

# Behavioural / Psychological
pgs$class[pgs$trait=="Educational attainment"] <- "Behavioural/\nPsychological"
pgs$class[pgs$trait=="Insomnia"] <- "Behavioural/\nPsychological"
pgs$class[pgs$trait=="Cognitive performance"] <- "Behavioural/\nPsychological"
pgs$class[pgs$trait=="Age smoking initiation"] <- "Behavioural/\nPsychological"
pgs$class[pgs$trait=="Cigarettes per day"] <- "Behavioural/\nPsychological"
pgs$class[pgs$trait=="Drinks per week"] <- "Behavioural/\nPsychological"
pgs$class[pgs$trait=="Risk tolerance"] <- "Behavioural/\nPsychological"


#participatoion
pgs$class[pgs$trait=="Participation in UKB aide-memoire"] <- "Participation"
pgs$class[pgs$trait=="Participation in UKB food questionnaire"] <- "Participation"
pgs$class[pgs$trait=="Completion of UKB mental health questionnaire"] <- "Participation"
pgs$class[pgs$trait=="Participation in UKB wearable"] <- "Participation"


#### plot ####
pgs$trait <- factor(pgs$trait, levels = sort(unique(pgs$trait)))
pgs$drug <- factor(pgs$drug, levels = c("Blood pressure\nmedications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer\nmedications"))
pgs$class <- factor(pgs$class, levels = sort(unique(pgs$class)))

fwrite(pgs, 'Projects/drugs/manuscript/pgs/pgs_adh_final.tsv', sep = "\t")

pgs <- pgs %>%
  group_by(drug) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"),
         significant = ifelse(p.value < 0.05, 1, 0),
         fdr_significant = ifelse(fdr < 0.05, 1, 0))

palette <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00")

text_size = 7
theme.size = 1.25
astrix.size = 2
line.size = 0.25

png("/Users/cordioli/Projects/drugs/manuscript/figures/SFx_b_PGS_adh.png", width = 12, height = 8, units = "in", res = 300)
ggplot(pgs) +
  facet_grid(factor(class) ~ factor(drug), scales = "free_y", space = "free_y") +
  geom_pointrange(aes(y=trait,x=perc_change,xmin=perc_lower, xmax=perc_upper, color = factor(fdr_significant)), size = .12) +
  theme_minimal() +
  scale_x_continuous(labels = scales::percent) +
  scale_colour_manual(values=c("gray80", "dodgerblue4"), name = "Significance", labels = c("Not significant", "Significant after FDR correction")) +
  geom_vline(aes(xintercept = 0)) +
  ylab("Trait PGS") +
  xlab("Percentage change in adherence (95% CI)\nper SD increase in PGS")
dev.off()