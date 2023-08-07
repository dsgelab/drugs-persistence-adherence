rm (list=ls())

library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)

options(bitmapType='cairo')

drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

all <- NULL

for (d in drugs) {
  adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',d,'.rds')) %>%
    mutate(drug = case_when (d == "statins" ~ "Statins",
                             d == "blood_pressure" ~ "Blood pressure medications",
                             d == "doac" ~ "Anticoagulants",
                             d == "antiplatelet" ~ "Antiplatelets",
                             d == "breast_cancer" ~ "Breast cancer medications")) %>% 
              select(drug, adherence)
  all <- bind_rows (all, adh)
}  

palette <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00")
meds <- c("Blood pressure\nmedications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer\nmedications")

# # # Panel a - proportion of persistence
df_per <- fread('/data/projects/project_adherence/results/manuscript/table1_FR.tsv') %>%
  mutate(drug = case_when (drug == "statins" ~ "Statins",
                           drug == "blood_pressure" ~ "Blood pressure\nmedications",
                           drug == "doac" ~ "Anticoagulants",
                           drug == "antiplatelet" ~ "Antiplatelets",
                           drug == "breast_cancer" ~ "Breast cancer\nmedications")) %>% 
  select(drug, prop = prop_persistent)

df_per$drug <- factor(df_per$drug, levels = meds)

panel_a <- ggplot (df_per, aes(y=prop, x=drug, fill=drug, label = scales::percent(prop, accuracy = 0.01))) +
  geom_bar(stat="identity", alpha = 0.9, width = .5) +
  ylab ("Proportion of persistent individuals\n(1+ year treatment)") +
  xlab ("Medication") +
  geom_text(position = position_dodge(width = .9), 
            vjust = -0.6,
            size = 4) + 
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  ggtitle("a.") +
  theme(panel.grid.minor.x = element_blank (),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title.position = "plot") +
  scale_fill_manual(values = palette)


# # # Panel b - proportion of good adherence
df_good <- all %>%
  group_by(drug) %>% 
  summarise(prop = sum(adherence >= 0.8)/n()) %>% 
  as.data.frame()

df_good$drug <- factor(df_good$drug, levels = meds)

panel_b <- ggplot (df_good, aes(y=prop, x=drug, fill=drug, label = scales::percent(prop, accuracy = 0.01))) +
  geom_bar(stat="identity", alpha = 0.9, width = .5) +
  ylab ("Proportion of good adherers\n(adherence > 0.8)") +
  xlab ("Medication") +
  geom_text(position = position_dodge(width = .9), 
            vjust = -0.6,
            size = 4) + 
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  ggtitle("b.") +
  theme(panel.grid.minor.x = element_blank (),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title.position = "plot") +
  scale_fill_manual(values = palette)


# order
meds <- c("Blood pressure medications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer medications")
all$drug <- factor(all$drug, levels = meds)

# DF for mean adherence
meandf <- all %>% 
  group_by(drug) %>% 
  summarise(mean_adh = mean(adherence))

# # # Panel c - adherence distribution
panel_c_free <- ggplot (all, aes (x = adherence, fill = drug)) +
  geom_histogram (bins = 200, alpha=0.9) +
  xlab ("Adherence") +
  ylab ("Number of individuals") +
  scale_x_continuous (breaks = c(0, 0.5, 0.8, 1, 1.1)) +
  geom_hline (yintercept = 0, alpha = 1, size = .1) +
  geom_vline(xintercept = 1, alpha = 1, size = .3) +
  geom_vline(data=meandf, aes(xintercept=mean_adh), alpha = .9, size = .3, linetype = "longdash") +
  ylim(c(0,NA)) +
  theme_minimal(base_size = 14) +
  ggtitle("c.") +
  theme(panel.grid.minor.x = element_blank (),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        plot.title.position = "plot",
        strip.background = element_rect(fill = alpha("steelblue", 0.2), color = NA)) +
  facet_wrap(~drug, scales = "free", ncol = 1, strip.position="top") +
  scale_fill_manual(values=palette)

png('/data/projects/project_adherence/results/manuscript/Figure2c.png', width = 5, height = 10, unit = "in", res = 300)
panel_c_free
dev.off()


# # # combine in biorender

# png('/data/projects/project_adherence/results/manuscript/Figure2a_scales_free.png', width = 11, height = 10, unit = "in", res = 300)
# grid.arrange(
#   panel_a, panel_b, panel_c_free,
#   widths = c(1.1, 0.9),
#   layout_matrix = rbind(c(1, 3),
#                         c(2, 3))
# )
# dev.off()
