# # # # # # # # #
### READ RGS ###
rm(list=ls())

get_rg <- function(x)
{
  # Get RGs from the log
  tt <- system( paste0('grep "^ /home/ivm/drugs/results/R8_20220112/munged_stop/" ', x), intern = T)
  
  df <- NULL
  
  if (length(tt)>0) {
  tt2 <- lapply(strsplit(tt," "), function(x) x[x != ""])
  
  df <- data.frame(matrix(unlist(tt2), nrow=length(tt2), byrow=T),stringsAsFactors=FALSE)
  
  # Get header
  h <- system( paste0('grep "p1" ', x), intern = T)
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

phenolist <- c(paste0(c("statins", "blood_pressure", "breast_cancer", "clopi_dipy", "doac", "glauc"),"_stop_1"),
               paste0(c("statins", "blood_pressure", "breast_cancer", "clopi_dipy", "doac", "glauc"),"_stop_2"),
               paste0(c("statins", "blood_pressure", "breast_cancer", "clopi_dipy", "doac", "glauc"),"_stop_3"))


h2 <- NULL
d <- NULL

for (p in phenolist) {
  temp <- get_rg(paste0('/home/ivm/drugs/results/R8_20220112/rg_stop/',p,'.af_0.01.info_0.6.gz.rsid.gz.sumstats.gz_RGS.log'))
  d <- rbind(d,temp)
  
  temp <- get_h2(paste0('/home/ivm/drugs/results/R8_20220112/h2_stop/',p,'.af_0.01.info_0.6.gz.rsid.gz.sumstats.gz_H2.log'))
  h2 <- rbind(h2, c(p,temp))
}

h2 <- as.data.frame(h2, stringsAsFactors = F)
colnames(h2) <- c("p","h2","se","pval")

d$p1 <- sub('/home/ivm/drugs/results/R8_20220112/munged_stop/','',d$p1)
d$p2 <- sub('/home/ivm/drugs/sumstats_for_rg/processed/','',d$p2)
cols.num <- seq(3,12)
d[cols.num] <- sapply(d[cols.num],as.numeric)

dd <- d
#dd <- d[-which(d$p2=="h2 | SE | p"),]

# traits to keep
# to.keep <- c("ADHD_2017.LDSC.sumstats.gz", "BMI.LDSC.sumstats.gz", "CAD.LDSC.sumstats.gz", "EA.LDSC.sumstats.gz", "IBD.LDSC.sumstats.gz",
#              "LDL.LDSC.sumstats.gz", "MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz", "Neuroticism_Full.LDSC.sumstats.gz", "SCZ2.LDSC.sumstats.gz",
#              "SWB_Full.LDSC.sumstats.gz", "T2D_2018.LDSC.sumstats.gz", "TG.LDSC.sumstats.gz", "UKB.self_rated_health.LDSC.sumstats.gz",
#              "alcohol_clarke.LDSC.sumstats.gz","anorexia_2019.LDSC.sumstats.gz","anxiety.LDSC.sumstats.gz","autism_2017.ipsych.pgc.LDSC.sumstats.gz",
#              "bipolar.LDSC.sumstats.gz","breast_cancer.LDSC.sumstats.gz","cannabis_ever_2018.no23andMe.LDSC.sumstats.gz",
#              "fasting_insuline.LDSC.sumstats.gz","loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz","openness.LDSC.sumstats.gz",
#              "risk_behavior.LDSC.sumstats.gz","smoking_cigs_per_day.LDSC.sumstats.gz",
#              "smoking_ever_vs_never.LDSC.sumstats.gz", "systolic_BP_93_irnt.ldsc.imputed_v3.both_sexes.tsv.gz", "UKB-ICBPmeta750k_SBP.sumstats.gz",
#              "Xue_et_al_T2D_META_Nat_Commun_2018.sumstats.gz", "BOLT_sex_combined_s340951_unrelated_WhiteBritish.LDL.sumstats.gz")

to.keep <- c("ADHD_2017.LDSC.sumstats.gz", "BMI.LDSC.sumstats.gz", "CAD.LDSC.sumstats.gz", "EA.LDSC.sumstats.gz", "IBD.LDSC.sumstats.gz",
             "LDL.LDSC.sumstats.gz", "MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz", "Neuroticism_Full.LDSC.sumstats.gz", "SCZ2.LDSC.sumstats.gz",
             "SWB_Full.LDSC.sumstats.gz", "T2D_2018.LDSC.sumstats.gz", "TG.LDSC.sumstats.gz", "UKB.self_rated_health.LDSC.sumstats.gz",
             "alcohol_clarke.LDSC.sumstats.gz","anorexia_2019.LDSC.sumstats.gz","anxiety.LDSC.sumstats.gz","autism_2017.ipsych.pgc.LDSC.sumstats.gz",
             "bipolar.LDSC.sumstats.gz","breast_cancer.LDSC.sumstats.gz","cannabis_ever_2018.no23andMe.LDSC.sumstats.gz",
             "fasting_insuline.LDSC.sumstats.gz","loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz","openness.LDSC.sumstats.gz",
             "risk_behavior.LDSC.sumstats.gz","smoking_cigs_per_day.LDSC.sumstats.gz",
             "smoking_ever_vs_never.LDSC.sumstats.gz", "UKB-ICBPmeta750k_SBP.sumstats.gz",
             "GCST90012790_buildGRCh38.tsv.gz.sumstats.gz", "GCST90012793_buildGRCh38.tsv.gz.sumstats.gz")

# ,"GCST90012792_buildGRCh38.tsv.gz.sumstats.gz", "GCST90012793_buildGRCh38.tsv.gz.sumstats.gz")

dd <- dd[dd$p2 %in% to.keep,]

# Assign class to each trait
# Psychiatric
dd$class[dd$p2=="ADHD_2017.LDSC.sumstats.gz"] <- "Psychiatric"
dd$class[dd$p2=="anorexia_2019.LDSC.sumstats.gz"] <- "Psychiatric"
dd$class[dd$p2=="anxiety.LDSC.sumstats.gz"] <- "Psychiatric"
dd$class[dd$p2=="autism_2017.ipsych.pgc.LDSC.sumstats.gz"] <- "Psychiatric"
dd$class[dd$p2=="bipolar.LDSC.sumstats.gz"] <- "Psychiatric"
dd$class[dd$p2=="MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz"] <- "Psychiatric"
dd$class[dd$p2=="SCZ2.LDSC.sumstats.gz"] <-"Psychiatric"

# Biomarkers
dd$class[dd$p2=="BMI.LDSC.sumstats.gz"] <- "Biomarkers"
dd$class[dd$p2=="fasting_insuline.LDSC.sumstats.gz"] <- "Biomarkers"
dd$class[dd$p2=="UKB-ICBPmeta750k_SBP.sumstats.gz"] <- "Biomarkers"
dd$class[dd$p2=="LDL.LDSC.sumstats.gz"] <- "Biomarkers"
dd$class[dd$p2=="TG.LDSC.sumstats.gz"] <- "Biomarkers"

# Disease Liability
dd$class[dd$p2=="CAD.LDSC.sumstats.gz"] <- "Disease Liability"
dd$class[dd$p2=="IBD.LDSC.sumstats.gz"] <- "Disease Liability"
dd$class[dd$p2=="breast_cancer.LDSC.sumstats.gz"] <- "Disease Liability"
dd$class[dd$p2=="T2D_2018.LDSC.sumstats.gz"] <- "Disease Liability"

# Psychological/Beahvioural
dd$class[dd$p2=="EA.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[dd$p2=="risk_behavior.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[dd$p2=="loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[dd$p2=="openness.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[dd$p2=="Neuroticism_Full.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[dd$p2=="SWB_Full.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[dd$p2=="UKB.self_rated_health.LDSC.sumstats.gz"] <- "Psychological/Beahvioural"
dd$class[which(startsWith(dd$p2, "GCS"))] <- "Psychological/Beahvioural"

# Substance Use
dd$class[dd$p2=="alcohol_clarke.LDSC.sumstats.gz"] <- "Substance Use"
dd$class[dd$p2=="cannabis_ever_2018.no23andMe.LDSC.sumstats.gz"] <- "Substance Use"
dd$class[dd$p2=="smoking_ever_vs_never.LDSC.sumstats.gz"] <- "Substance Use"
dd$class[dd$p2=="smoking_cigs_per_day.LDSC.sumstats.gz"] <- "Substance Use"


# Rename traits
dd$p2[dd$p2=="ADHD_2017.LDSC.sumstats.gz"] <- "ADHD"
dd$p2[dd$p2=="BMI.LDSC.sumstats.gz"] <- "BMI"
dd$p2[dd$p2=="CAD.LDSC.sumstats.gz"] <- "Coronary Artery Disease"
dd$p2[dd$p2=="EA.LDSC.sumstats.gz"] <- "Educational Attainment"
dd$p2[dd$p2=="IBD.LDSC.sumstats.gz"] <- "Inflammatory Bowel Disease"
dd$p2[dd$p2=="LDL.LDSC.sumstats.gz"] <- "LDL Cholesterol"
dd$p2[dd$p2=="SWB_Full.LDSC.sumstats.gz"] <- "Subjective Well-being"
dd$p2[dd$p2=="T2D_2018.LDSC.sumstats.gz"] <- "Type 2 Diabetes"
dd$p2[dd$p2=="TG.LDSC.sumstats.gz"] <- "Triglycerides"
dd$p2[dd$p2=="UKB.self_rated_health.LDSC.sumstats.gz"] <- "Self-rated Health"
dd$p2[dd$p2=="MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz"] <- "Depression"
dd$p2[dd$p2=="Neuroticism_Full.LDSC.sumstats.gz"] <- "Neuroticism"
dd$p2[dd$p2=="SCZ2.LDSC.sumstats.gz"] <- "Schizophrenia"
dd$p2[dd$p2=="alcohol_clarke.LDSC.sumstats.gz"] <- "Alcohol Consumption"
dd$p2[dd$p2=="anorexia_2019.LDSC.sumstats.gz"] <- "Anorexia"
dd$p2[dd$p2=="anxiety.LDSC.sumstats.gz"] <- "Anxiety"
dd$p2[dd$p2=="autism_2017.ipsych.pgc.LDSC.sumstats.gz"] <- "Autism"
dd$p2[dd$p2=="bipolar.LDSC.sumstats.gz"] <- "Bipolar Disorder"
dd$p2[dd$p2=="breast_cancer.LDSC.sumstats.gz"] <- "Breast Cancer"
dd$p2[dd$p2=="cannabis_ever_2018.no23andMe.LDSC.sumstats.gz"] <- "Cannabis Use (ever)"
dd$p2[dd$p2=="fasting_insuline.LDSC.sumstats.gz"] <- "Fasting Insulin"
dd$p2[dd$p2=="loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz"] <- "Loneliness"
dd$p2[dd$p2=="openness.LDSC.sumstats.gz"] <- "Openness"
dd$p2[dd$p2=="risk_behavior.LDSC.sumstats.gz"] <- "Risk Tolerance"
dd$p2[dd$p2=="smoking_ever_vs_never.LDSC.sumstats.gz"] <- "Smoking (ever)"
dd$p2[dd$p2=="smoking_cigs_per_day.LDSC.sumstats.gz"] <- "Cigarettes Per Day"
dd$p2[dd$p2=="UKB-ICBPmeta750k_SBP.sumstats.gz"] <- "Systolic Blood Pressure"
dd$p2[dd$p2=="GCST90012790_buildGRCh38.tsv.gz.sumstats.gz"] <- "Participation to UKB Food Q."
dd$p2[dd$p2=="GCST90012791_buildGRCh38.tsv.gz.sumstats.gz"] <- "Participation to UKB Phisical Activity Monitoring"
dd$p2[dd$p2=="GCST90012792_buildGRCh38.tsv.gz.sumstats.gz"] <- "Participation to UKB Mental Health Q."
dd$p2[dd$p2=="GCST90012793_buildGRCh38.tsv.gz.sumstats.gz"] <- "Participation to UKB Aide-memoire"



dd$p_fdr <- p.adjust(dd$p, "fdr")
dd$p_label <- ""
dd$p_label[dd$p_fdr < 0.05] <- "*"
dd$sig <- 0
dd$sig[dd$p_fdr < 0.05] <- 1

dd$p1 <- sub(".af_0.01.info_0.6.gz.rsid.gz.sumstats.gz","",dd$p1)
dd$stop <- substr(dd$p1, nchar(dd$p1), nchar(dd$p1))
dd$p1 <- sub("_.*", "", dd$p1)

dd$p1[dd$p1=="statins"] <- "Statins"
dd$p1[dd$p1=="blood"] <- "BP Medications"
dd$p1[dd$p1=="clopi"] <- "Clopidogrel+Dipyridamol"
dd$p1[dd$p1=="breast"] <- "Beast Cancer Medications"
dd$p1[dd$p1=="doac"] <- "DOAC"
dd$p1[dd$p1=="glauc"] <- "Glaucoma Medications"
#dd$p1[dd$p1=="_combo_adherence_std.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz"] <- "Combined"

dd$p2 <- factor(dd$p2, levels=sort(unique(dd$p2), decreasing = T))

ddd <- dd[dd$p1 %in% c("Statins","BP Medications"),]
#ddd <- dd

ddd$p1 <- factor(ddd$p1, levels=c("BP Medications", "Statins"))

library(ggplot2)

png("/home/ivm/drugs/results/R8_20220112/RGS_stop_all_meds_R8.png", width = 18, height = 11, units = "in", res = 300)
ggplot(dd, aes(x=rg, y=p2, colour = factor(sig), shape = factor(stop), group = factor(stop, levels = c(3,2,1)))) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=rg-1.96*se, xmax=rg+1.96*se), height=0, position=position_dodge(width=0.5)) +
  # geom_text(aes(label=p_label), color = "red", size = 5) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(class)~factor(p1), scales = "free_y", space = "free_y", switch = "y") +
  ylab("") +
  xlab("Genetic correlation with stopping after n purchases\n[-1: less likely to stop, 1: more likely to stop]") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("", "Significant after FDR correction")) +
  theme(#legend.position = "none",
        text = element_text(size=17),
        strip.placement = "outside",
        #strip.background =element_rect(fill=alpha('dodgerblue4', 0.5)),
        strip.text.y.left = element_text(angle = 0, hjust = 1))
dev.off()

png("/home/ivm/drugs/results/R8_20220112/RGS_stop_R8.png", width = 13, height = 10, units = "in", res = 300)
ggplot(ddd, aes(x=rg, y=p2, colour = factor(sig), shape = factor(stop), group = factor(stop, levels = c(3,2,1)))) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbarh(aes(xmin=rg-1.96*se, xmax=rg+1.96*se), height=0, position=position_dodge(width=0.5)) +
  # geom_text(aes(label=p_label), color = "red", size = 5) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(factor(class)~factor(p1), scales = "free_y", space = "free_y", switch = "y") +
  ylab("") +
  xlab("Genetic correlation with stopping after n purchases\n[-1: less likely to stop, 1: more likely to stop]") +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("", "Significant after FDR correction")) +
  theme(#legend.position = "none",
    text = element_text(size=17),
    strip.placement = "outside",
    #strip.background =element_rect(fill=alpha('dodgerblue4', 0.5)),
    strip.text.y.left = element_text(angle = 0, hjust = 1))
dev.off()

library(data.table)
fwrite(ddd, '/home/ivm/drugs/results/R8_20220112/RGS_finngen_R8.tsv', sep = "\t", quote = F)

d_fdr <- dd[dd$p_label=="*", c("p1","p2","rg","p","p_fdr")]
