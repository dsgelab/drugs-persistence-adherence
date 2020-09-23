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
  d <- rbind(d,readRG(paste0('/home/cordioli/drugs/results/rg_04092020/',p,'_RGS.log')))
}

d$p1 <- sub('/home/cordioli/drugs/results/','',d$p1)
d$p2 <- sub('/home/cordioli/drugs/sumstats_for_rg/processed/','',d$p2)
cols.num <- seq(3,12)
d[cols.num] <- sapply(d[cols.num],as.numeric)

dd <- d[-which(d$p2=="h2 | SE | p"),]

# traits to keep
to.keep <- c("ADHD_2017.LDSC.sumstats.gz", "BMI.LDSC.sumstats.gz", "CAD.LDSC.sumstats.gz", "EA.LDSC.sumstats.gz", "IBD.LDSC.sumstats.gz",
             "LDL.LDSC.sumstats.gz", "MDD2018_ex23andMe.sumstats.LDSC.sumstats.gz", "Neuroticism_Full.LDSC.sumstats.gz", "SCZ2.LDSC.sumstats.gz",
             "SWB_Full.LDSC.sumstats.gz", "T2D_2018.LDSC.sumstats.gz", "TG.LDSC.sumstats.gz", "UKB.self_rated_health.LDSC.sumstats.gz",
             "alcohol_clarke.LDSC.sumstats.gz","anorexia_2019.LDSC.sumstats.gz","anxiety.LDSC.sumstats.gz","autism_2017.ipsych.pgc.LDSC.sumstats.gz",
             "bipolar.LDSC.sumstats.gz","breast_cancer.LDSC.sumstats.gz","cannabis_ever_2018.no23andMe.LDSC.sumstats.gz",
             "fasting_insuline.LDSC.sumstats.gz","loneliness.rapid_UKB.sumstats.LDSC.sumstats.gz","openness.LDSC.sumstats.gz",
             "risk_behavior.LDSC.sumstats.gz","smoking_cigs_per_day.LDSC.sumstats.gz",
             "smoking_ever_vs_never.LDSC.sumstats.gz")

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




dd$p_fdr <- p.adjust(dd$p, "fdr")
dd$p_label <- ""
dd$p_label[dd$p_fdr < 0.05] <- "*"
dd$sig <- 0
dd$sig[dd$p_fdr < 0.05] <- 1

dd$p1[dd$p1=="statins.sumstats.gz"] <- "Statins"
dd$p1[dd$p1=="blood_pressure.sumstats.gz"] <- "Blood Pressure Medications"

dd$p2 <- factor(dd$p2, levels=sort(unique(dd$p2), decreasing = T))


library(ggplot2)

pdf("/home/cordioli/drugs/results/rg_04092020/rgs2.pdf", 5, 6)
ggplot(dd, aes(x=rg, y=p2, colour = factor(sig))) +
  geom_point() +
  geom_errorbarh(aes(xmin=rg-1.96*se, xmax=rg+1.96*se), height=0.1) +
  # geom_text(aes(label=p_label), color = "red", size = 5) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(.~factor(p1)) +
  ylab("Trait") +
  xlab("Genetic correlation with adherence") +
  scale_colour_manual(values=c("#999999", "dodgerblue3")) +
  theme(legend.position = "none")
dev.off()


d_fdr <- dd[dd$p_label=="*", c("p1","p2","rg","p","p_fdr")]
