# # # # # # # # #
### READ RGS ###
rm(list=ls())

get_rg <- function(x)
{
  # Get RGs from the log
  tt <- system( paste0('grep "^/home/cordioli/drugs/results/20210422_R7/munged/" ', x), intern = T)
  
  df <- NULL
  
  if (length(tt)>0) {
  tt2 <- lapply(strsplit(tt," "), function(x) x[x != ""])
  
  df <- data.frame(matrix(unlist(tt2), nrow=length(tt2), byrow=T),stringsAsFactors=FALSE)
  
  # Get header
  h <- system( paste0('grep "^p1" ', x), intern = T)
  h <- unlist(strsplit(h," "))
  h <- h[h!=""]
  
  colnames(df) <- h
  }
  
  return(df)
}

get_h2 <- function(x)
{
  # Get RGs from the log
  tt <- system( paste0('grep "^Total Observed scale h2:" ', x), intern = T)
  
  r <- c(NA,NA,NA)
  
  if (length(tt)>0) {
    tt2 <- lapply(strsplit(tt," "), function(x) x[x != ""])
    
    h2 <- as.numeric(tt2[[1]][5])
    se <- as.numeric(gsub("[()]","",tt2[[1]][6]))
    p <- 2*pnorm(-abs(h2/se))
    
    r <- c(h2,se,p)
  }
  
  return(r)
}


phenolist <- readLines('/home/cordioli/drugs/data/adherence_phenolist.txt')
phenolist <- c(phenolist, 'combo_adherence_std')
# phenolist <- c('A10B', 'BP_MED', 'C10')

h2 <- NULL
d <- NULL
for (p in phenolist) {
  temp <- get_rg(paste0('/home/cordioli/drugs/results/20210422_R7/rg/finngen_R7_',p,'.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz_RGS.log'))
  d <- rbind(d,temp)
  
  temp <- get_h2(paste0('/home/cordioli/drugs/results/20210422_R7/h2/finngen_R7_',p,'.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz_H2.log'))
  h2 <- rbind(h2, c(p,temp))
}

h2 <- as.data.frame(h2, stringsAsFactors = F)
colnames(h2) <- c("p","h2","se","pval")

d$p1 <- sub('/home/cordioli/drugs/results/20210422_R7/munged/','',d$p1)
d$p2 <- sub('/home/cordioli/drugs/sumstats_for_rg/processed/','',d$p2)
cols.num <- seq(3,12)
d[cols.num] <- sapply(d[cols.num],as.numeric)

dd <- d
#dd <- d[-which(d$p2=="h2 | SE | p"),]

# traits to keep
to.keep <- c("ADHD_2017.LDSC.sumstats.gz", "BMI.LDSC.sumstats.gz", "CAD.LDSC.sumstats.gz", "EA.LDSC.sumstats.gz", "IBD.LDSC.sumstats.gz",
             "LDL.LDSC.sumstats.gz", "MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz", "Neuroticism_Full.LDSC.sumstats.gz", "SCZ2.LDSC.sumstats.gz",
             "SWB_Full.LDSC.sumstats.gz", "T2D_2018.LDSC.sumstats.gz", "TG.LDSC.sumstats.gz", "UKB.self_rated_health.LDSC.sumstats.gz",
             "alcohol_clarke.LDSC.sumstats.gz","anorexia_2019.LDSC.sumstats.gz","anxiety.LDSC.sumstats.gz","autism_2017.ipsych.pgc.LDSC.sumstats.gz",
             "bipolar.LDSC.sumstats.gz","breast_cancer.LDSC.sumstats.gz","cannabis_ever_2018.no23andMe.LDSC.sumstats.gz",
             "fasting_insuline.LDSC.sumstats.gz","loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz","openness.LDSC.sumstats.gz",
             "risk_behavior.LDSC.sumstats.gz","smoking_cigs_per_day.LDSC.sumstats.gz",
             "smoking_ever_vs_never.LDSC.sumstats.gz", "systolic_BP_93_irnt.ldsc.imputed_v3.both_sexes.tsv.gz", "UKB-ICBPmeta750k_SBP.sumstats.gz",
             "Xue_et_al_T2D_META_Nat_Commun_2018.sumstats.gz", "BOLT_sex_combined_s340951_unrelated_WhiteBritish.LDL.sumstats.gz")

dd <- dd[dd$p2 %in% to.keep,]

dd$p2[dd$p2=="ADHD_2017.LDSC.sumstats.gz"] <- "ADHD"
dd$p2[dd$p2=="BMI.LDSC.sumstats.gz"] <- "BMI"
dd$p2[dd$p2=="CAD.LDSC.sumstats.gz"] <- "CAD"
dd$p2[dd$p2=="EA.LDSC.sumstats.gz"] <- "Educational Attainment"
dd$p2[dd$p2=="IBD.LDSC.sumstats.gz"] <- "IBD"
dd$p2[dd$p2=="LDL.LDSC.sumstats.gz"] <- "LDL"
dd$p2[dd$p2=="SWB_Full.LDSC.sumstats.gz"] <- "Subjective Well-being"
dd$p2[dd$p2=="T2D_2018.LDSC.sumstats.gz"] <- "T2D"
dd$p2[dd$p2=="TG.LDSC.sumstats.gz"] <- "Triglycerides"
dd$p2[dd$p2=="UKB.self_rated_health.LDSC.sumstats.gz"] <- "Self-rated Health"
dd$p2[dd$p2=="MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz"] <- "Depression"
dd$p2[dd$p2=="Neuroticism_Full.LDSC.sumstats.gz"] <- "Neuroticism"
dd$p2[dd$p2=="SCZ2.LDSC.sumstats.gz"] <- "Schizophrenia"
dd$p2[dd$p2=="alcohol_clarke.LDSC.sumstats.gz"] <- "Alcohol consumption"
dd$p2[dd$p2=="anorexia_2019.LDSC.sumstats.gz"] <- "Anorexia"
dd$p2[dd$p2=="anxiety.LDSC.sumstats.gz"] <- "Anxiety"
dd$p2[dd$p2=="autism_2017.ipsych.pgc.LDSC.sumstats.gz"] <- "Autism"
dd$p2[dd$p2=="bipolar.LDSC.sumstats.gz"] <- "Bipolar Disorder"
dd$p2[dd$p2=="breast_cancer.LDSC.sumstats.gz"] <- "Breast Cancer"
dd$p2[dd$p2=="cannabis_ever_2018.no23andMe.LDSC.sumstats.gz"] <- "Cannabis use (ever)"
dd$p2[dd$p2=="fasting_insuline.LDSC.sumstats.gz"] <- "Fasting Insuline"
dd$p2[dd$p2=="loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz"] <- "Loneliness"
dd$p2[dd$p2=="openness.LDSC.sumstats.gz"] <- "Openness"
dd$p2[dd$p2=="risk_behavior.LDSC.sumstats.gz"] <- "Risk Behaviour"
dd$p2[dd$p2=="smoking_ever_vs_never.LDSC.sumstats.gz"] <- "Smoking (ever)"
dd$p2[dd$p2=="smoking_cigs_per_day.LDSC.sumstats.gz"] <- "Smoking (CPD)"
dd$p2[dd$p2=="systolic_BP_93_irnt.ldsc.imputed_v3.both_sexes.tsv.gz"] <- "Systolic Blood Pressure"



dd$p_fdr <- p.adjust(dd$p, "fdr")
dd$p_label <- ""
dd$p_label[dd$p_fdr < 0.05] <- "*"
dd$sig <- 0
dd$sig[dd$p_fdr < 0.05] <- 1


dd$p1[dd$p1=="finngen_R7_statins.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "Statins"
dd$p1[dd$p1=="finngen_R7_blood_pressure.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "BP Medications"
dd$p1[dd$p1=="finngen_R7_clopi_dipy.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "Clopidogrel+Dipyridamol"
dd$p1[dd$p1=="finngen_R7_breast_canc.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "Beast Cancer Medications"
dd$p1[dd$p1=="finngen_R7_doac.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "DOAC"
dd$p1[dd$p1=="finngen_R7_glauc.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "Glaucoma Medications"
dd$p1[dd$p1=="finngen_R7_combo_adherence_std.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "Combined"

dd$p2 <- factor(dd$p2, levels=sort(unique(dd$p2), decreasing = T))

ddd <- dd[dd$p1 %in% c("Statins","BP Medications", "Combined"),]
#ddd <- dd

ddd$p1 <- factor(ddd$p1, levels=c("BP Medications", "Statins", "Combined"))

library(ggplot2)

pdf("/home/cordioli/drugs/results/20210422_R7/RGS.pdf", width = 13, height = 8)
ggplot(ddd, aes(x=rg, y=p2, colour = factor(sig))) +
  geom_point() +
  geom_errorbarh(aes(xmin=rg-1.96*se, xmax=rg+1.96*se), height=0.1) +
  # geom_text(aes(label=p_label), color = "red", size = 5) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(.~factor(p1)) +
  ylab("") +
  xlab("Genetic correlation") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("", "Significant after FDR correction")) +
  theme(legend.position = "none",
        text = element_text(size=18))
dev.off()


d_fdr <- dd[dd$p_label=="*", c("p1","p2","rg","p","p_fdr")]
