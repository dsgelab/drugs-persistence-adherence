rm(list=ls())

library(data.table)
library(dplyr)
library(broom)
library(gridExtra)
library(ggplot2)
library(forcats)
library(gtools)

res_adh <- fread("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence_all.tsv") %>% 
  mutate(pheno = "persistence",
         low = estimate-1.96*std.error,
         up = estimate+1.96*std.error,
         label = case_when(drug == "breast_canc" ~ "CYP2D6\nTamoxifen",
                           drug == "clopi_dipy" ~ "CYP2C19\nClopidogrel",
                           drug == "statins" ~ "SLCO1B1\nStatins (all)",
                           drug == "fluvastatins" ~ "CYP2C9\nFluvastatins",
                           drug == "rosuvastatins" ~ "ABCG2\nRosuvastatins"),
         Phenotype = case_when(term == "poor_function" ~ "Poor function",
                               term == "decreased_function" ~ "Decreased function",
                               term == "increased_function" ~ "Increased function",
                               term == "intermediate_metabolizer" ~ "Intermediate metabolizer",
                               term == "poor_metabolizer" ~ "Poor metabolizer",
                               term == "ultrarapid_metabolizer" ~ "Ultrarapid metabolizer"),
         Phenotype = factor(Phenotype, levels = c("Poor function", "Decreased function", "Increased function",
                                                  "Poor metabolizer", "Intermediate metabolizer", "Ultrarapid metabolizer")),
         label = factor(label, levels = rev(c("CYP2D6\nTamoxifen",
                                              "CYP2C19\nClopidogrel",
                                              "ABCG2\nRosuvastatins",
                                              "CYP2C9\nFluvastatins",
                                              "SLCO1B1\nStatins (all)"))),
         sign = ifelse(p.value < 0.05, 1, 0))

res_per <- fread("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence_all.tsv") %>% 
  mutate(pheno = "persistence",
         low = exp(estimate-1.96*std.error),
         up = exp(estimate+1.96*std.error),
         OR = exp(estimate),
         label = case_when(drug == "breast_canc" ~ "CYP2D6\nTamoxifen",
                           drug == "clopi_dipy" ~ "CYP2C19\nClopidogrel",
                           drug == "statins" ~ "SLCO1B1\nStatins (all)",
                           drug == "fluvastatins" ~ "CYP2C9\nFluvastatins",
                           drug == "rosuvastatins" ~ "ABCG2\nRosuvastatins"),
         Phenotype = case_when(term == "poor_function" ~ "Poor function",
                               term == "decreased_function" ~ "Decreased function",
                               term == "increased_function" ~ "Increased function",
                               term == "intermediate_metabolizer" ~ "Intermediate metabolizer",
                               term == "poor_metabolizer" ~ "Poor metabolizer",
                               term == "ultrarapid_metabolizer" ~ "Ultrarapid metabolizer"),
         Phenotype = factor(Phenotype, levels = c("Poor function", "Decreased function", "Increased function",
                                                  "Poor metabolizer", "Intermediate metabolizer", "Ultrarapid metabolizer")),
         label = factor(label, levels = rev(c("CYP2D6\nTamoxifen",
                                              "CYP2C19\nClopidogrel",
                                              "ABCG2\nRosuvastatins",
                                              "CYP2C9\nFluvastatins",
                                              "SLCO1B1\nStatins (all)"))),
         sign = ifelse(p.value < 0.05, 1, 0)) %>% 
  filter(!grepl("Clopidogrel", label))

# Persistence
# pdf('/home/ivm/drugs/results/manuscript/pharmacogenes/pgx_forest_plot_persistence.pdf', height = 5, width = 5)
# ggplot(res_per, aes(x = OR, y = Phenotype, xmin = 0, xmax = 2, color = factor(sign))) +
#   geom_point(size = 1) +
#   geom_errorbarh(aes(xmin = low, xmax = up), height = 0) +
#   facet_wrap(~label,  ncol=1, scales = "free") +
#   theme_minimal() +
#   labs(
#     x = "Odds of persistence (95% CI)",
#     y = "Phenotype") +
#   scale_colour_manual(values=c("gray80", "dodgerblue4"), name = "Significance",
#                       labels = c("Not significant", "Significant at P < 0.05"))
# dev.off()


pdf('/home/ivm/drugs/results/manuscript/pharmacogenes/pgx_forest_plot_persistence.pdf', height = 5.5, width = 5.5)
ggplot(res_per, aes(x = OR, y = Phenotype, xmin = 0, xmax = 2)) +
  geom_vline(xintercept = 1) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = low, xmax = up), height = 0) +
  facet_grid(label~., scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Odds of persistence (95% CI)",
    y = "Phenotype") +
  theme(text=element_text(size=14))
dev.off()

# Adherence
pdf('/home/ivm/drugs/results/manuscript/pharmacogenes/pgx_forest_plot_adherence.pdf', height = 5.5, width = 5.5)
ggplot(res_adh, aes(x = estimate, y = Phenotype, xmin = -0.2, xmax = 0.2)) +
  geom_vline(xintercept = 0) +
  geom_point(size = 1) +
  geom_errorbarh(aes(xmin = low, xmax = up), height = 0) +
  facet_grid(label~., scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Percentage change in adherence (95% CI)",
    y = "Phenotype") +
  theme(text=element_text(size=14))
dev.off()


res_per2 <- res_per %>% select(label, N_tot, Phenotype, N_dis, N_per, OR, low, up, p.value) 
