rm(list = ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)
library(ggforestplot)

source('/Users/cordioli/Projects/drugs/manuscript/predictors/geom_stripes.R')

df <- fread('/Users/cordioli/Projects/drugs/manuscript/predictors/predictors.binary.model.adherence.summary.tsv')
df2 <- fread('/Users/cordioli/Projects/drugs/manuscript/predictors/predictors.binary.model.persistence.summary.tsv')

df <- df %>%
  filter(!term %in% c("(Intercept)", "yob")) %>% 
  mutate(lower = estimate - 1.96*std.error,
         upper = estimate + 1.96*std.error,
         per_adh = estimate/1.1,
         per_adh_lower = lower/1.1,
         per_adh_upper = upper/1.1)

mother_tongue <- df %>% 
  filter(term == "mother_tongue_fin_sv1")
mean(mother_tongue$per_adh)
min(mother_tongue$per_adh)
max(mother_tongue$per_adh)

soc <- df %>% 
  filter(term == "soc_prev_year1")
mean(soc$per_adh)
min(soc$per_adh)
max(soc$per_adh)


df2 <- df2 %>%
  filter(!term %in% c("(Intercept)", "yob")) %>% 
  mutate(OR = exp(estimate),
         lower = exp(estimate - 1.96*std.error),
         upper = exp(estimate + 1.96*std.error))

sex <- df2 %>% 
  filter(term == "sex1")
mean(sex$OR)
min(sex$OR)
max(sex$OR)

df$trait[df$term == "age_first_purch_601"] <- "Age at initiation > 60"
df$trait[df$term == "sex1"] <- "Sex (female) *"
df$trait[df$term == "secondary1"] <- "Secondary prevention **"
df$trait[df$term == "CCI_score_severe_51"] <- "Charlson comorbidity index > 5\n(severe)"
df$trait[df$term == "mother_tongue_fin_sv1"] <- "Mother tongue (Fin/Swe)"
df$trait[df$term == "edu_years_bachelor_161"] <- "Education (University degree)"
df$trait[df$term == "urban1"] <- "Living in urban areas"
df$trait[df$term == "soc_prev_year1"] <- "Social assistance benefits"

df$Medication[df$drug == "statins"] <- "Statins"
df$Medication[df$drug == "doac"] <- "Anticoagulants"
df$Medication[df$drug == "breast_cancer"] <- "Breast cancer\nmedications"
df$Medication[df$drug == "blood_pressure"] <- "Blood pressure\nmedications"
df$Medication[df$drug == "antiplatelet"] <- "Antiplatelets"


df2$trait[df2$term == "age_first_purch_601"] <- "Age at initiation > 60"
df2$trait[df2$term == "sex1"] <- "Sex (female) *"
df2$trait[df2$term == "secondary1"] <- "Secondary prevention **"
df2$trait[df2$term == "CCI_score_severe_51"] <- "Charlson comorbidity index > 5\n(severe)"
df2$trait[df2$term == "mother_tongue_fin_sv1"] <- "Mother tongue (Fin/Swe)"
df2$trait[df2$term == "edu_years_bachelor_161"] <- "Education (University degree)"
df2$trait[df2$term == "urban1"] <- "Living in urban areas"
df2$trait[df2$term == "soc_prev_year1"] <- "Social assistance benefits"

df2$Medication[df2$drug == "statins"] <- "Statins"
df2$Medication[df2$drug == "doac"] <- "Anticoagulants"
df2$Medication[df2$drug == "breast_cancer"] <- "Breast cancer\nmedications"
df2$Medication[df2$drug == "blood_pressure"] <- "Blood pressure\nmedications"
df2$Medication[df2$drug == "antiplatelet"] <- "Antiplatelets"

# set order
drugs = c("Blood pressure\nmedications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer\nmedications")
df$Medication <- factor(df$Medication, levels = drugs)
df2$Medication <- factor(df2$Medication, levels = drugs)

df$Significance <- factor(ifelse(df$p.value < 0.05, "Significant at P < 0.05", "Not significant"), levels = c("Significant at P < 0.05", "Not significant"))


palette <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00")


png('/Users/cordioli/Projects/drugs/manuscript/figures/Figure3b.png', width = 7, height = 5, res = 300, units = "in")
ggforestplot::forestplot(
  df = df,
  estimate = per_adh,
  se = std.error/1.1,
  name = trait,
  pvalue = p.value,
  psignif = 0.05,
  colour = Medication,
  xlab = "Percentage change in adherence (95% CI)") +
  scale_color_manual(values=palette) +
  scale_x_continuous(labels = scales::percent) +
  ggtitle("b.") +
  theme(plot.title.position = "plot",
        legend.position = "none")
dev.off()

# remove antiplatelets
df2 <- df2 %>% filter(Medication != "Antiplatelets")
palette <- c("#999999", "#56B4E9", "#E69F00", "#D55E00")

png('/Users/cordioli/Projects/drugs/manuscript/figures/Figure3a.png', width = 7, height = 5, res = 300, units = "in")
ggforestplot::forestplot(
  df = df2,
  estimate = estimate,
  se = std.error,
  name = trait,
  pvalue = p.value,
  psignif = 0.05,
  colour = Medication,
  logodds = TRUE,
  xlab = "Odds of persistence (95% CI)") +
  scale_color_manual(values=palette) +
  ggtitle("a.") +
  theme(plot.title.position = "plot",
        legend.position = "none")
dev.off()


# print extra one for legend
png('/Users/cordioli/Projects/drugs/manuscript/figures/Figure3_legend.png', width = 7, height = 5, res = 300, units = "in")
ggforestplot::forestplot(
  df = df,
  estimate = per_adh,
  se = std.error/1.1,
  name = trait,
  pvalue = p.value,
  shape = Significance,
  psignif = 0.05,
  colour = Medication,
  xlab = "Percentage change in adherence (95% CI)") +
  scale_color_manual(values=palette) +
  scale_x_continuous(labels = scales::percent) +
  ggtitle("a.") +
  theme(plot.title.position = "plot") +
  ggplot2::scale_shape_manual(
    values = c(21L, 1L)) +
  guides(shape = guide_legend(override.aes = list(size=2)),
         color = guide_legend(override.aes = list(size=2)))

dev.off()