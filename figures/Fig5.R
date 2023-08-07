rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggplot2)
library(forcats)
library(cowplot)

traits <- c("ADHD",
            "Alzheimer's disease",
            "Cognitive performance",
            "Stroke",
            "LDL cholesterol",
            "Depression",
            "Diabetes",
            "Asthma",
            "Bipolar disorder",
            "BMI",
            "Breast cancer",
            "Coronary artery disease",
            "Diastolic blood pressure",
            "Drinks per week",
            "Educational attainment",
            "Glycated hemoglobin",
            "Insomnia",
            "LDL cholesterol adjusted*",
            "Lifespan",
            "Chronic pain",
            "Neuroticism",
            "Participation in UKB** aide-memoire",
            "Participation in UKB** food questionnaire",
            "Completion of UKB** mental health questionnaire",
            "Participation in UKB** wearable",
            "Prostate cancer",
            "Rheumatoid arthritis",            
            "Risk tolerance",            
            "Schizophrenia",
            "Age smoking initiation",
            "Cigarettes per day",
            "Systolic blood pressure",
            "Triglycerides"
)

drugs <- c("Statins",
           "Blood pressure medications",
           "Breast cancer medications",
           "Antiplatelets",
           "Anticoagulants")

alpha <- 0.05
multiplier <- qnorm(1 - alpha / 2)

# Read in adherence / stop file
rgs_adh <- fread("/Users/cordioli/Projects/drugs/manuscript/rgs/rgs_adherence.tsv") %>% 
  mutate(pheno = "Adherence",
         p1 = ifelse(p1 == "breast_canc", "breast_cancer", p1))

rgs_per <- fread("/Users/cordioli/Projects/drugs/manuscript/rgs/rgs_persistence.tsv") %>% 
  mutate(pheno = "Persistence")

rgs <- bind_rows(rgs_adh, rgs_per)

traits_mapping <- data.frame(p2 = unique(rgs$p2),
                         trait = traits)

drugs_mapping <- data.frame(p1 = unique(rgs$p1),
                            drug = drugs)

rgs <- rgs %>% 
  left_join(traits_mapping) %>% 
  left_join(drugs_mapping) %>% 
  mutate(lower = rg - (multiplier * se),
         upper = rg + (multiplier * se),
         trait = factor(trait, levels=sort(traits)),
         drug = factor(drug),
         fdr = p.adjust(p, method = "fdr"))

#### Assign class to each trait #### 
# Psychiatric
rgs$class[rgs$trait=="Depression"] <- "Psychiatric"
rgs$class[rgs$trait=="Neuroticism"] <- "Psychiatric"
rgs$class[rgs$trait=="Schizophrenia"] <- "Psychiatric"
rgs$class[rgs$trait=="Bipolar disorder"] <- "Psychiatric"
rgs$class[rgs$trait=="ADHD"] <- "Psychiatric"

# Biomarker
rgs$class[rgs$trait=="BMI"] <- "Biomarker"
rgs$class[rgs$trait=="Systolic blood pressure"] <- "Biomarker"
rgs$class[rgs$trait=="Diastolic blood pressure"] <- "Biomarker"
rgs$class[rgs$trait=="Glycated hemoglobin"] <- "Biomarker"
rgs$class[rgs$trait=="Triglycerides"] <- "Biomarker"
rgs$class[rgs$trait=="LDL cholesterol"] <- "Biomarker"
rgs$class[rgs$trait=="LDL cholesterol adjusted*"] <- "Biomarker"

# Disease liability
rgs$class[rgs$trait=="Coronary artery disease"] <- "Disease liability"
rgs$class[rgs$trait=="Stroke"] <- "Disease liability"
rgs$class[rgs$trait=="Atrial Fibrillation"] <- "Disease liability"
rgs$class[rgs$trait=="Alzheimer's disease"] <- "Disease liability"
rgs$class[rgs$trait=="Rheumatoid arthritis"] <- "Disease liability"
rgs$class[rgs$trait=="Asthma"] <- "Disease liability"
rgs$class[rgs$trait=="Prostate cancer"] <- "Disease liability"
rgs$class[rgs$trait=="Breast cancer"] <- "Disease liability"
rgs$class[rgs$trait=="Irritable bowel syndrome"] <- "Disease liability"
rgs$class[rgs$trait=="Glaucoma"] <- "Disease liability"
rgs$class[rgs$trait=="Diabetes"] <- "Disease liability"
rgs$class[rgs$trait=="Chronic pain"] <- "Disease liability"
rgs$class[rgs$trait=="Lifespan"] <- "Disease liability"

# Behavioural / Psychological
rgs$class[rgs$trait=="Educational attainment"] <- "Behavioural / Psychological"
rgs$class[rgs$trait=="Insomnia"] <- "Behavioural / Psychological"
rgs$class[rgs$trait=="Cognitive performance"] <- "Behavioural / Psychological"
rgs$class[rgs$trait=="Age smoking initiation"] <- "Behavioural / Psychological"
rgs$class[rgs$trait=="Cigarettes per day"] <- "Behavioural / Psychological"
rgs$class[rgs$trait=="Drinks per week"] <- "Behavioural / Psychological"
rgs$class[rgs$trait=="Risk tolerance"] <- "Behavioural / Psychological"


#participatoion
rgs$class[rgs$trait=="Participation in UKB** aide-memoire"] <- "Participation"
rgs$class[rgs$trait=="Participation in UKB** food questionnaire"] <- "Participation"
rgs$class[rgs$trait=="Completion of UKB** mental health questionnaire"] <- "Participation"
rgs$class[rgs$trait=="Participation in UKB** wearable"] <- "Participation"


#### plot ####
rgs$trait <- factor(rgs$trait, levels = sort(unique(rgs$trait)))
rgs$drug <- factor(rgs$drug, levels = rev(c("Blood pressure medications","Statins","Antiplatelets","Anticoagulants","Breast cancer medications")))
rgs$class <- factor(rgs$class, levels = sort(unique(rgs$class)))

fwrite(rgs, 'Projects/drugs/manuscript/rgs/rgs_both_final.tsv', sep = "\t")

rgs_all <- rgs %>%
  mutate(drug = factor(drug, levels = rev(c("Blood pressure medications","Statins","Antiplatelets","Anticoagulants","Breast cancer medications")))) %>% 
  group_by(drug, pheno) %>% 
  mutate(fdr = p.adjust(p, method = "fdr"),
         significant = ifelse(p < 0.05, 1, 0),
         fdr_significant = ifelse(fdr < 0.05, 1, 0))

rgs_filter <- rgs %>%
  filter(drug %in% c("Blood pressure medications", "Statins", "Antiplatelets"))  %>% 
  mutate(drug = factor(drug, levels = rev(c("Blood pressure medications", "Statins", "Antiplatelets"))),
         significant = ifelse(p < 0.05, 1, 0),
         tile_size = ifelse(p > 0.85, 0.15, 1-p)) %>% 
  group_by(drug, pheno) %>% 
  mutate(fdr = p.adjust(p, method = "fdr"),
         fdr_significant = ifelse(fdr < 0.05, 1, 0))

rgs_filter$pheno <- factor(rgs_filter$pheno, levels = rev(c("Adherence", "Persistence")))

text_size = 7
theme.size = 1.25
astrix.size = 2
line.size = 0.25


rgs_filter <- rgs_filter %>% filter(!(drug == "Antiplatelets" & pheno == "Persistence"))

p <- ggplot(rgs_filter) +
  facet_grid(factor(pheno) ~ factor(class), scales = "free", space = "free") +
  # geom_raster(data = filter(dat.long, !is.na(rg)), aes(x = exposure.name, y = outcome.name, fill = rg)) +
  geom_tile(data = rgs_filter,
            aes(x = trait, y = drug, fill = rg, height = tile_size, width = tile_size)) +
  # geom_text(data = rgs_filter,
  #           aes(label = fdr_label, x = trait, y = drug), vjust = 0.75) +
  #scale_fill_distiller(palette = rev("RdBu"), limits = c(-1, 1)) +
  # scale_fill_gradientn(
  #   colours = RColorBrewer::brewer.pal(11, "RdBu"),
  #   breaks = seq(-1, 1, length.out = 11)
  # ) +
  scale_fill_gradient2(low="#67001F",
                       high="#053061",
                       mid="grey95",
                       name = "Genetic correlation",
                       limits = c(-1, 1)) +
  geom_point(data = subset(rgs_filter, significant == 1), 
             aes(x = trait, y = drug),
             size = 0.08) +
  geom_point(data = subset(rgs_filter, fdr_significant == 1),
             aes(x = trait, y = drug),
             size = 1.5,
             shape = 1) +
  theme_classic() +
  geom_vline(xintercept=seq(0.5, 33.5, 1),color="grey90",size=line.size) +
  geom_hline(yintercept=seq(0.5, 3.5, 1),color="grey90",size=line.size) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 35, hjust = 0),
        #aspect.ratio=1,
        legend.text = element_text(hjust = 1.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=margin(10,10,10,10, 'pt'),
        text = element_text(size=text_size, family="Arial"),
        strip.text.y = element_text(size = text_size-1, family="Arial"),
        axis.ticks = element_line(size = line.size),
        axis.line = element_line(size = line.size),
        strip.background = element_blank(),
        # panel.border = element_rect(colour = "black", fill = NA)
        # legend.key.height = unit(0.75, "line"),
        # legend.spacing.y = unit(-0.75, "line"),
        # plot.margin=grid::unit(c(0,0,0,0), "mm")
  ) +
  scale_x_discrete(position = "top") +
  ylim(rev(levels(rgs_filter$drug)))

png("/Users/cordioli/Projects/drugs/manuscript/figures/Figure5_panels.png", width = 7, height = 2.6, units = "in", res = 600)
p
dev.off()

png("/Users/cordioli/Projects/drugs/manuscript/figures/Figure5_legend.png", width = 10, height = 3, units = "in", res = 600)
p + theme(legend.position = 'right')
dev.off()

png("/Users/cordioli/Projects/drugs/manuscript/figures/Figure5_sign_legend.png", width = 2, height = 1.5, units = "in", res = 600)
df.legend <- data.frame(x = c(0, 0),
                        y = c(0.5,1.5),
                        significant = c(1,1),
                        fdr_significant = c(1,0),
                        labels = c("Significant\nafter FDR correction", "Significant\nat P < 0.05"))

ggplot(df.legend, aes(x = x, y = y, label = labels)) +
  geom_point(data = subset(df.legend, significant == 1), 
             size = 0.5) +
  geom_point(data = subset(df.legend, fdr_significant == 1),
             shape = 1,
             size = 2.5) +
  geom_text(aes(x=0.2), hjust = 0) +
  ylim(c(0,2)) +
  xlim(c(-.1,3)) +
  theme_void() +
  theme(legend.position = 'none',
        text = element_text(family="Arial"),
        strip.background = element_blank())
dev.off()