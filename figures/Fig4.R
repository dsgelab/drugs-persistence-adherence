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

freq <- fread('/home/ivm/drugs/results/manuscript/pharmacogenes/diplotype_freq.tsv')

res_adh_2 <- res_adh %>% left_join(freq, by = c("gene" = "gene", "star_allele" = "diplotype")) %>% filter(complete.cases(.))

res_per_2 <- res_per %>% left_join(freq, by = c("gene" = "gene", "star_allele" = "diplotype")) %>% filter(complete.cases(.))

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


# # # QQ plots

# Statins
plot_adh_statins <- res_adh_statins %>% 
  mutate(observed_log10 = -log10(p.value),
         observed = p.value) %>% 
  arrange(observed)
plot_adh_statins$expected <- ppoints(nrow(plot_adh_statins))
plot_adh_statins$expected_log10 <- -log10(ppoints(nrow(plot_adh_statins)))

# Clopidogrel
plot_adh_clopidogrel <- res_adh_clopidogrel %>% 
  mutate(observed_log10 = -log10(p.value),
         observed = p.value) %>% 
  arrange(observed)
plot_adh_clopidogrel$expected <- ppoints(nrow(plot_adh_clopidogrel))
plot_adh_clopidogrel$expected_log10 <- -log10(ppoints(nrow(plot_adh_clopidogrel)))

# Tamoxifen
plot_adh_tamoxifen <- res_adh_tamoxifen %>% 
  mutate(observed_log10 = -log10(p.value),
         observed = p.value) %>% 
  arrange(observed)
plot_adh_tamoxifen$expected <- ppoints(nrow(plot_adh_tamoxifen))
plot_adh_tamoxifen$expected_log10 <- -log10(ppoints(nrow(plot_adh_tamoxifen)))


plot_adh <- bind_rows(plot_adh_statins, plot_adh_clopidogrel, plot_adh_tamoxifen)

plot_adh$drug[plot_adh$drug == "statins"] <- "Statins"
plot_adh$drug[plot_adh$drug == "breast_canc"] <- "Tamoxifen"
plot_adh$drug[plot_adh$drug == "clopi_dipy"] <- "Clopidogrel"

plot_adh$gene <- toupper(plot_adh$gene)

png('/home/ivm/drugs/results/manuscript/pharmacogenes/pharmacogenes_adherence_qqplot_log10.png', height = 4, width = 9, units = "in", res = 300)
ggplot(plot_adh, aes(x=expected_log10, y=observed_log10, shape = factor(gene), col = factor(drug))) +
  geom_point(size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linewidth = .5) +
  coord_equal() +
  theme_minimal() +
  xlab( expression( paste("Expected ", -log[10](italic(P)))) ) +
  ylab( expression( paste("Observed ", -log[10](italic(P)))) ) +
  labs(col = "Drug", shape = "Gene") +
  facet_wrap(~factor(drug), nrow = 1) +
  scale_color_manual(values = c("#009E73", "#56B4E9", "#D55E00")) +
  scale_shape_manual(values = c(15,16,18,3,4)) +
  theme(strip.background = element_rect(fill = alpha("steelblue", 0.2), color = NA))
dev.off()



