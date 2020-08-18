# conda activate ldsc
# 
# sumstats_list=$(ls -d /home/cordioli/drugs/sumstats_for_rg/*gz | tr '\n' ',')
# sumstats_list=${sumstats_list::-1}
# 
# /home/cordioli/ldsc/ldsc.py \
# --rg /home/cordioli/drugs/results/statins.sumstats.gz,${sumstats_list} \
# --ref-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --w-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/statins_RGS
# 
# 
# /home/cordioli/ldsc/ldsc.py \
# --rg /home/cordioli/drugs/results/blood_pressure.sumstats.gz,${sumstats_list} \
# --ref-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --w-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/blood_pressure_RGS
# 
# 
# /home/cordioli/ldsc/ldsc.py \
# --rg /home/cordioli/drugs/results/clopi_dipy.sumstats.gz,${sumstats_list} \
# --ref-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --w-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/clopi_dipy_RGS
# 
# 
# /home/cordioli/ldsc/ldsc.py \
# --rg /home/cordioli/drugs/results/breast_canc.sumstats.gz,${sumstats_list} \
# --ref-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --w-ld-chr /home/cordioli/ldsc/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/breast_canc_RGS


# # # # # # # # #

### READ RGS ###
rm(list=ls())

readRG <- function(log) {
  rg_log <- readLines(log)
  n_pheno <- as.integer(strsplit(rg_log[grep("Computing rg for phenotype 2/",rg_log)],
                                 "/")[[1]][2])
  
  h2 <- rg_log[grep("Heritability of phenotype 1$",rg_log)+2]
  h2 <- strsplit(h2," ")
  h2_beta <- as.numeric(h2[[1]][5])
  h2_se <- as.numeric(gsub("\\(|\\)","",h2[[1]][6]))
  h2_p <- 2*pnorm(-abs(h2_beta/h2_se))
  
  rg_idx <- as.integer(grep("Summary of",rg_log)+1)
  rg <- rg_log[rg_idx:(rg_idx+n_pheno)]
  
  rg <- strsplit(rg, " ")
  
  getNotNull <- function(x) {
    x[x!=""]
  }
  
  rrg <- sapply(rg, getNotNull)
  df <- as.data.frame(t(matrix(unlist(rrg),ncol=n_pheno)),stringsAsFactors = F)
  colnames(df) <- df[1,]
  df <- df[-1,]
  rownames(df) <- NULL
  df <- df[-(n_pheno-1),]
  row_h2 <- c(df$p1[1],"h2 | SE | p", h2_beta, h2_se, h2_p,"-","-","-","-","-","-","-")
  df <- rbind(row_h2,df)
  return(df)
}

phenolist <- readLines('/home/cordioli/drugs/data/adherence_phenolist.txt')
d <- NULL
for (p in phenolist) {
  d <- rbind(d,readRG(paste0('/home/cordioli/drugs/results/',p,'_RGS.log')))
}

d$p1 <- sub('/home/cordioli/drugs/results/','',d$p1)
d$p2 <- sub('/home/cordioli/drugs/sumstats_for_rg/','',d$p2)
cols.num <- seq(3,12)
d[cols.num] <- sapply(d[cols.num],as.numeric)

dd <- d[-which(d$p2=="h2 | SE | p"),]

dd$p_label <- paste0('p: ', round(dd$p, 4))

dd <- dd[which(dd$p2 %in% c("Educational_Attainment_excl23andme_Lee_2018_NatGen_formatted.sumstats.gz",
                           "systolic_BP_93_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz",
                           "LDL_UKBB_lifted.sumstats.gz",
                           "depressive_symptoms.sumstats.gz",
                           "risk_taking.sumstats.gz",
                           "schizophrenia.sumstats.gz",
                           "subjective_well_being.sumstats.gz")), ]
 
library(ggplot2)
ggplot(dd, aes(x=rg, y=factor(p2))) +
  geom_point() +
  geom_errorbarh(aes(xmin=rg-1.96*se, xmax=rg+1.96*se), height=0.1) +
  # geom_text(aes(label=p_label), vjust=-1) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(.~factor(p1))