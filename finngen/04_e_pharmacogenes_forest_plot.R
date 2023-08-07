rm(list=ls())

library(data.table)
library(dplyr)
library(gridExtra)
library(ggplot2)

res_adh_files <- system("ls /home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence*", intern = T)
res_per_files <- system("ls /home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence*", intern = T)

res_adh <- rbindlist(lapply(res_adh_files, function(name) fread(name))) %>% 
  mutate(pheno = "adherence",
         low = estimate-1.96*std.error,
         up = estimate+1.96*std.error,
         beta = estimate,
         fdr_significant = ifelse(fdr < 0.05, 1, 0))

res_per <- rbindlist(lapply(res_per_files, function(name) fread(name))) %>% 
  mutate(pheno = "persistence") %>% 
  mutate(low = exp(estimate-1.96*std.error),
         up = exp(estimate+1.96*std.error),
         OR = exp(estimate),
         fdr_significant = ifelse(fdr < 0.05, 1, 0))


# # # Drug-gene specific results for adherence
res_adh <- res_adh %>% 
  filter(complete.cases(.))

res_adh_statins <- res_adh %>%
  filter(drug == "statins",
         gene %in% c("slco1b1", "cyp2c9", "abcg2")) %>% 
  mutate(star_allele = factor(star_allele, levels = rev(sort(unique(star_allele)))))

res_adh_clopidogrel <- res_adh %>%
  filter(drug == "clopi_dipy",
         gene == "cyp2c19") %>% 
  mutate(star_allele = factor(star_allele, levels = rev(sort(unique(star_allele)))))

res_adh_tamoxifen <- res_adh %>%
  filter(drug == "breast_canc",
         gene == "cyp2d6") %>% 
  mutate(star_allele = factor(star_allele, levels = rev(sort(unique(star_allele)))))

res_per <- res_per %>% 
  filter(complete.cases(.))

res_per_statins <- res_per %>%
  filter(drug == "statins",
         gene %in% c("slco1b1", "cyp2c9", "abcg2")) %>% 
  mutate(star_allele = factor(star_allele, levels = rev(sort(unique(star_allele)))))

res_per_clopidogrel <- res_per %>%
  filter(drug == "clopi_dipy",
         gene == "cyp2c19") %>% 
  mutate(star_allele = factor(star_allele, levels = rev(sort(unique(star_allele)))))

res_per_tamoxifen <- res_per %>%
  filter(drug == "breast_canc",
         gene == "cyp2d6") %>% 
  mutate(star_allele = factor(star_allele, levels = rev(sort(unique(star_allele)))))


# # # Adherence forest plots
statins <- ggplot(res_adh_statins, aes(y = star_allele, color = factor(fdr_significant))) +
  geom_point(aes(x=beta), shape=15, size=2) +
  geom_linerange(aes(xmin=low, xmax=up)) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(gene)~., scales = "free", space = "free") +
  xlab("logOR (95% CI)") +
  ylab("") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("Not significant\nafter FDR correction", "Significant after FDR correction")) +
  ggtitle("a. statins") +
  theme(legend.position = "none") +
  xlim(-.15,.15)


png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_adherence_statins_forest_plot.png', height = 7, width = 5, units = "in", res = 300)
print(statins)
dev.off()

clopidogrel <- ggplot(res_adh_clopidogrel, aes(y = star_allele, color = factor(fdr_significant))) +
  geom_point(aes(x=beta), shape=15, size=2) +
  geom_linerange(aes(xmin=low, xmax=up)) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(gene)~., scales = "free", space = "free") +
  xlab("logOR (95% CI)") +
  ylab("") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("Not significant\nafter FDR correction", "Significant after FDR correction")) +
  ggtitle("b. clopidogrel") +
  theme(legend.position = "none") +
  xlim(-.1,.1)

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_adherence_clopidogrel_forest_plot.png', height = 5, width = 5, units = "in", res = 300)
print(clopidogrel)
dev.off()

tamoxifen <- ggplot(res_adh_tamoxifen, aes(y = star_allele, color = factor(fdr_significant))) +
  geom_point(aes(x=beta), shape=15, size=2) +
  geom_linerange(aes(xmin=low, xmax=up)) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(gene)~., scales = "free", space = "free") +
  xlab("logOR (95% CI)") +
  ylab("") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("Not significant\nafter FDR correction", "Significant after FDR correction")) +
  ggtitle("c. tamoxifen") +
  theme(legend.position = "none") +
  xlim(-.1,.1)

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_adherence_tamoxifen_forest_plot.png', height = 6, width = 5, units = "in", res = 300)
print(tamoxifen)
dev.off()

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_adherence_all_forest_plot.png', height = 7, width = 5, units = "in", res = 300)
grid.arrange(statins, clopidogrel, tamoxifen,
  layout_matrix = rbind(c(1, 2),
                        c(1, 3)))
dev.off()


# # # Persistence forest plots
statins <- ggplot(res_per_statins, aes(y = star_allele, color = factor(fdr_significant))) +
  geom_point(aes(x=OR), shape=15, size=2) +
  geom_linerange(aes(xmin=low, xmax=up)) +
  geom_vline(aes(xintercept = 1), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(gene)~., scales = "free", space = "free") +
  xlab("OR (95% CI)") +
  ylab("") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "",
                      labels = c("Not significant\nafter FDR correction", "Significant after FDR correction")) +
  ggtitle("a. statins") +
  theme(legend.position = "none") +
  xlim(0,3.5)

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_persistence_statins_forest_plot.png', height = 7, width = 6, units = "in", res = 300)
print(statins)
dev.off()

clopidogrel <- ggplot(res_per_clopidogrel, aes(y = star_allele, color = factor(fdr_significant))) +
  geom_point(aes(x=OR), shape=15, size=2) +
  geom_linerange(aes(xmin=low, xmax=up)) +
  geom_vline(aes(xintercept = 1), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(gene)~., scales = "free", space = "free") +
  xlab("OR (95% CI)") +
  ylab("") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "",
                      labels = c("Not significant\nafter FDR correction", "Significant after FDR correction")) +
  ggtitle("b. clopidogrel") +
  theme(legend.position = "none") +
  xlim(.1,2.6)

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_persistence_clopidogrel_forest_plot.png', height = 5, width = 5, units = "in", res = 300)
print(clopidogrel)
dev.off()

tamoxifen <- ggplot(res_per_tamoxifen, aes(y = star_allele, color = factor(fdr_significant))) +
  geom_point(aes(x=OR), shape=15, size=2) +
  geom_linerange(aes(xmin=low, xmax=up)) +
  geom_vline(aes(xintercept = 1), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(gene)~., scales = "free", space = "free") +
  xlab("OR (95% CI)") +
  ylab("") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "",
                      labels = c("Not significant\nafter FDR correction", "Significant after FDR correction")) +
  ggtitle("c. tamoxifen") +
  theme(legend.position = "none") +
  xlim(.1,2.6)

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_persistence_tamoxifen_forest_plot.png', height = 6, width = 8, units = "in", res = 300)
print(tamoxifen)
dev.off()

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_persistence_all_forest_plot.png', height = 7, width = 5, units = "in", res = 300)
grid.arrange(statins, clopidogrel, tamoxifen,
             layout_matrix = rbind(c(1, 2),
                                   c(1, 3)))
dev.off()


# # # QQ plots
adh_tot <- bind_rows(res_adh_statins, res_adh_clopidogrel, res_adh_tamoxifen)

plotd <- adh_tot %>% 
  mutate(observed_log10 = -log10(fdr),
         observed = fdr,
         dir = ifelse(beta > 0, "pos", "neg")) %>% 
  arrange(observed)

plotd$expected <- ppoints(nrow(plotd))
plotd$expected_log10 <- -log10(ppoints(nrow(plotd)))

