rm(list = ls())

library(data.table)
library(ggplot2)
library(dplyr)

st <- fread('~/drugs/data/R7_statins_summarized_150.txt')

prs_edu <- fread('~/drugs/prs/r7/finngen_R7_Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted.txt.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_risk_tol <- fread('~/drugs/prs/r7/finngen_R7_RISK_GWAS_MA_UKB+23andMe+replication.txt.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_ldl <- fread('~/drugs/prs/r7/finngen_R7_continuous-LDLC-both_sexes-medadj_irnt.tsv.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_sbp <- fread('~/drugs/prs/r7/finngen_R7_UKB-ICBPmeta750k_SBPsummaryResults.txt.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_t2d <- fread('~/drugs/prs/r7/finngen_R7_METAANALYSIS_DIAGRAM_SE1.txt.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_bmi <- fread('~/drugs/prs/r7/finngen_R7_SNP_gwas_mc_merge_nogc.tbl.uniq.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_depression <- fread('~/drugs/prs/r7/finngen_R7_MDD2018_ex23andMe.19fields.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_smoke <- fread('~/drugs/prs/r7/finngen_R7_GSCAN_SmokingInitiation.txt.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_scz <- fread('~/drugs/prs/r7/finngen_R7_daner_PGC_SCZ52_0513a.resultfiles_PGC_SCZ52_0513.sh2_nofin.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))

prs_bip <- fread('~/drugs/prs/r7/finngen_R7_Bipolar_vs_control_PGC_Cell_2018_formatted.txt.sscore') %>% 
  mutate(prs = SCORE1_AVG,
         prs_std = scale(prs))



for(df in dfs){
  st_prs <- st %>%
    left_join(get(df), by = c("FINNGENID" = "IID"))
  
  s <- summary(lm(formula = as.formula(paste0("adherence_std ~ AGE_AT_DEATH_OR_END_OF_FOLLOWUP + SEX_IMPUTED +", paste(paste0("PC",seq(1,10)), collapse = "+"), "+ prs_std")), data = st_prs))

  p1 <- ggplot(st_prs, aes(x=prs, y=adherence)) + 
    geom_point(alpha=0.4, colour = "dodgerblue4") +
    geom_smooth(method='lm', formula= y~x) +
    xlab("PRS")+
    ylab("Adherence to statins")+
    theme_minimal()
  
  library(gridExtra)
  png(paste0("~/drugs/results/adherence_vs_prs/finngen_R7_statins_", df, ".png"), width = 5, height = 5, res = 200, units = "in")
  print(p1)
  dev.off()
}