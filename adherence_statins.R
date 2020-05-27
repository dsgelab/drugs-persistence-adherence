rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

source('/home/cordioli/drugs/adherence_funs.R')

purch <- fread('/home/cordioli/drugs/data/R5_v3_purch_vnr_98.gz')
ep <- fread('/home/cordioli/R5_pheno/finngen_R5_v3_endpoint.gz')
covs <- fread('/home/cordioli/drugs/data/R5_cov.txt')

# # # # # # # #
#   STATINS   #
# # # # # # # #

# Define 'chronic' users: events that lead to a strong need of statins:
# I9_ASO
# I9_CHD
# I9_ATHSCLE
# I9_CEREBVASC	Cerebrovascular diseases	
# I9_INTRACRA	Nontraumatic intracranial haemmorrhage	
# I9_SAH	Subarachnoid haemmorrhage	
# I9_ICH	Intracerebral haemmorrhage	
# I9_OTHINTRACRA	Other intracranial haemorrhages	
# I9_STR	Stroke, excluding SAH	
# I9_STR_SAH	Stroke, including SAH	
# I9_STR_EXH	Ischaemic Stroke, excluding all haemorrhages	
# I9_STENOSIS	Occlusion and stenosis of arteries, not leading to stroke
ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
                'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Individual trajectories
st <- getTrajectories(purch,'^C10AA')

# Summarised trajectories
stt <- summarizeTrajectories(st)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

stt <- stt %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

# Plots
library(cowplot)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
source("R/RainCloudPlots/tutorial_R/R_rainclouds.R")

p1 <- ggplot(data = stt, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ggtitle('Overall') +
  theme_minimal()

p2 <- ggplot(data = stt, aes(y = adherence, x = factor(chronic), fill = factor(chronic))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_discrete(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ggtitle('Stratified - starting after CVD event') +
  theme_minimal()
grid.arrange(p1, p2, nrow = 2)

palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00")
p3 <- ggplot(stt, aes(x=age_bin, y=adherence, fill=age_bin)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  ggtitle('Adherence') +
  theme_minimal() +
  scale_fill_manual(values=palette) +
  theme(legend.position = "none")

p4 <- ggplot(stt, aes(x=age_bin, fill=age_bin)) +
  geom_bar(stat = 'count') +
  ggtitle('N individuals') +
  labs(fill = "Age at first purch.") +
  scale_x_discrete(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values=palette) +
  theme_minimal()
plot_grid(p4, p3, align = "h", ncol = 2, rel_widths = c(2/3, 1/3))

# Correlation adherence - age 1st purchase
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(stt$adherence, stt$age_first_purch), 4) ), 
                          x = 0.1, y = 0.97, hjust = 0,
                          gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))
p1 <- ggplot(stt, aes(x=adherence, y=age_first_purch)) + 
  geom_point(alpha = .3) + 
  ggtitle("Adherence vs Age at 1st purchase") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "Adherence") + 
  scale_y_continuous(name = "Age at 1st purchase") + 
  annotation_custom(grob) + 
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black"))

# Correlation adherence - tot days
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(stt$adherence, stt$tot_days), 4) ), 
                          x = 0.1, y = 0.97, hjust = 0,
                          gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))
p2 <- ggplot(stt, aes(x=adherence, y=tot_days/max(stt$tot_days))) + 
  geom_point(alpha = .3) + 
  ggtitle("Adherence vs Tot days") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "Adherence") + 
  scale_y_continuous(name = "Tot days") + 
  annotation_custom(grob) + 
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black"))

# Correlation adherence - SD days
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(stt$adherence, stt$sd_days_norm), 4) ), 
                          x = 0.1, y = 0.97, hjust = 0,
                          gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))
p3 <- ggplot(stt, aes(x=adherence, y=sd_days_norm)) + 
  geom_point(alpha = .3) + 
  ggtitle("Adherence vs SD(days between purchases)") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "Adherence") + 
  scale_y_continuous(name = "SD days") + 
  annotation_custom(grob) + 
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black"))
grid.arrange(p1, p2, p3, ncol = 3)

# # # # Infer daily dose, examples
# par(mfrow=c(1,3))
# case_oneb <- st[st$FINNGENID == 'FGJN6YVWVN', ]
# d <- density(case_oneb$pills_norm, na.rm = T)
# plot(d)
# case_half <- st[st$FINNGENID == 'FGCRS4JHU2', ]
# d <- density(case_half$pills_norm, na.rm = T)
# plot(d)
# p <- c(runif(30, min=0.35, max=0.6), runif(26, min=0.8, max=1.01))
# d <- density(p, na.rm = T)
# plot(d)


# # # 
# Merge phenotypes with covariates for GWAS
ad <- stt %>%
  mutate(statins = adherence*10,
         age_first_statins = age_first_purch) %>%
  select(FINNGENID,statins,age_first_statins)

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt', sep = '\t', quote = F)


# # # # # # # # # # # # # # # # # #  
#   BLOOD PRESSURE MEDICATIONS    #
# # # # # # # # # # # # # # # # # #

# Define 'chronic' users: events that lead to a strong need of bp:

ep_chronic <- c('I9_ASO', 'I9_CHD', 'I9_ATHSCLE', 'I9_CEREBVASC', 'I9_INTRACRA', 'I9_SAH',
                'I9_ICH', 'I9_OTHINTRACRA', 'I9_STR', 'I9_STR_SAH', 'I9_STR_EXH', 'I9_STENOSIS')

# Individual trajectories
st <- getTrajectories(purch,'^C10AA')

# Summarised trajectories
stt <- summarizeTrajectories(st)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

stt <- stt %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

# Plots
library(cowplot)
library(dplyr)
library(readr)
library(grid)
library(gridExtra)
source("R/RainCloudPlots/tutorial_R/R_rainclouds.R")

p1 <- ggplot(data = stt, aes(y = adherence, x = -.2)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_continuous(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ggtitle('Overall') +
  theme_minimal()

p2 <- ggplot(data = stt, aes(y = adherence, x = factor(chronic), fill = factor(chronic))) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  scale_x_discrete(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ggtitle('Stratified - starting after CVD event') +
  theme_minimal()
grid.arrange(p1, p2, nrow = 2)

palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00")
p3 <- ggplot(stt, aes(x=age_bin, y=adherence, fill=age_bin)) + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_boxplot(width = .1, guides = FALSE, alpha = 0.5) +
  coord_flip() +
  ggtitle('Adherence') +
  theme_minimal() +
  scale_fill_manual(values=palette) +
  theme(legend.position = "none")

p4 <- ggplot(stt, aes(x=age_bin, fill=age_bin)) +
  geom_bar(stat = 'count') +
  ggtitle('N individuals') +
  labs(fill = "Age at first purch.") +
  scale_x_discrete(name = "") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_fill_manual(values=palette) +
  theme_minimal()
plot_grid(p4, p3, align = "h", ncol = 2, rel_widths = c(2/3, 1/3))

# Correlation adherence - age 1st purchase
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(stt$adherence, stt$age_first_purch), 4) ), 
                         x = 0.1, y = 0.97, hjust = 0,
                         gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))
p1 <- ggplot(stt, aes(x=adherence, y=age_first_purch)) + 
  geom_point(alpha = .3) + 
  ggtitle("Adherence vs Age at 1st purchase") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "Adherence") + 
  scale_y_continuous(name = "Age at 1st purchase") + 
  annotation_custom(grob) + 
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black"))

# Correlation adherence - tot days
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(stt$adherence, stt$tot_days), 4) ), 
                         x = 0.1, y = 0.97, hjust = 0,
                         gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))
p2 <- ggplot(stt, aes(x=adherence, y=tot_days/max(stt$tot_days))) + 
  geom_point(alpha = .3) + 
  ggtitle("Adherence vs Tot days") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "Adherence") + 
  scale_y_continuous(name = "Tot days") + 
  annotation_custom(grob) + 
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black"))

# Correlation adherence - SD days
grob = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(stt$adherence, stt$sd_days_norm), 4) ), 
                         x = 0.1, y = 0.97, hjust = 0,
                         gp = gpar(col = "blue", fontsize = 11, fontface = "bold")))
p3 <- ggplot(stt, aes(x=adherence, y=sd_days_norm)) + 
  geom_point(alpha = .3) + 
  ggtitle("Adherence vs SD(days between purchases)") + 
  geom_smooth(method=lm, se=FALSE) + 
  scale_x_continuous(name = "Adherence") + 
  scale_y_continuous(name = "SD days") + 
  annotation_custom(grob) + 
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), axis.line = element_line(color="black"), axis.line.x = element_line(color="black"))
grid.arrange(p1, p2, p3, ncol = 3)

# # # # Infer daily dose, examples
# par(mfrow=c(1,3))
# case_oneb <- st[st$FINNGENID == 'FGJN6YVWVN', ]
# d <- density(case_oneb$pills_norm, na.rm = T)
# plot(d)
# case_half <- st[st$FINNGENID == 'FGCRS4JHU2', ]
# d <- density(case_half$pills_norm, na.rm = T)
# plot(d)
# p <- c(runif(30, min=0.35, max=0.6), runif(26, min=0.8, max=1.01))
# d <- density(p, na.rm = T)
# plot(d)


# # # 
# Merge phenotypes with covariates for GWAS
ad <- stt %>%
  mutate(statins = adherence*10,
         age_first_statins = age_first_purch) %>%
  select(FINNGENID,statins,age_first_statins)

covs <- covs %>%
  left_join(ad)

fwrite(covs, '/home/cordioli/drugs/data/R5_cov_pheno_adherence.txt', sep = '\t', quote = F)



phenolist <- c('adherence')
fwrite(list(phenolist), '/home/cordioli/drugs/data/adherence_phenolist.txt', col.names = F)