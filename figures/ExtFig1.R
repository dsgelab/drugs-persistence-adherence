rm (list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)


drugs <- c('statins', 'blood_pressure', 'antiplatelet', 'breast_cancer', 'doac')

all <- NULL

for (d in drugs) {
  adh <- readRDS(paste0('/data/projects/project_adherence/data/trajectories/trajectories.summarized.',d,'.rds')) %>%
    mutate(drug = case_when (d == "statins" ~ "Statins",
                             d == "blood_pressure" ~ "Blood pressure medications",
                             d == "doac" ~ "Anticoagulants",
                             d == "antiplatelet" ~ "Antiplatelets",
                             d == "breast_cancer" ~ "Breast cancer medications")) %>% 
    select(drug, adherence, tot_days, tot_pills)
  all <- bind_rows (all, adh)
}  

# order
meds <- c("Blood pressure medications", "Statins", "Antiplatelets", "Anticoagulants", "Breast cancer medications")
all$drug <- factor(all$drug, levels = meds)
palette <- c("#999999", "#56B4E9", "#009E73", "#E69F00", "#D55E00")

scale_min <- round(min(all$adherence), digits = 1)

hex_dens_plot <- function(df){
  p <- ggplot(df) +
    stat_summary_hex(aes(x=tot_days, y=tot_pills, z=adherence), alpha = .75,
                     fun = function(x) if(abs(sum(x)) > 5) {mean(x)} else {NA},
                     bins = 100) +
    geom_abline(color = "gray50") +
    geom_density_2d(aes(x=tot_days, y=tot_pills), size = .15, color = "dodgerblue4") +
    scale_fill_continuous(type = "viridis",
                          name = "Adherence",
                          limits = c(0, 1.1)) +
    coord_equal() +
    theme_minimal() +
    xlab("Expected (N days)") +
    ylab("Observed (N pills)") +
    ggtitle(unique(df$drug))

}

p1 <- hex_dens_plot(all %>% filter(drug == "Blood pressure medications")) + theme(legend.position = "none")
p2 <- hex_dens_plot(all %>% filter(drug == "Statins")) + theme(legend.position = "none")
p3 <- hex_dens_plot(all %>% filter(drug == "Antiplatelets")) + theme(legend.position = "none")
p4 <- hex_dens_plot(all %>% filter(drug == "Anticoagulants")) + theme(legend.position = "none")
p5 <- hex_dens_plot(all %>% filter(drug == "Breast cancer medications")) + theme(legend.position = "none")

p6 <- hex_dens_plot(all %>% filter(drug == "Blood pressure medications"))

legend <- get_legend(p6)  

png('/data/projects/project_adherence/results/manuscript/Figure2d_2rows.png', width = 10, height = 5, unit = "in", res = 300)
grid.arrange(p1,p2,p3,p4,p5,legend, nrow = 2)
dev.off()
