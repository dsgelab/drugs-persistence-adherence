# # # # # # # # #
### READ RGS ###
rm(list=ls())

get_rg <- function(x)
{
  # Get RGs from the log
  tt <- system( paste0('grep "^ /home/ivm/drugs/results/R8_20220112/munged/" ', x), intern = T)
  
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

phenolist <- c("statins", "blood_pressure", "breast_canc", "clopi_dipy", "doac", "glauc")


h2 <- NULL
d <- NULL

for (p in phenolist) {
  temp <- get_rg(paste0('/home/ivm/drugs/results/R8_20220112/rg_adh_stop/',p,'.af_0.01.info_0.6.gz.rsid.gz.sumstats.gz_RGS.log'))
  d <- rbind(d,temp)
}

d$p1 <- sub('/home/ivm/drugs/results/R8_20220112/munged/','',d$p1)
d$p2 <- sub('/home/ivm/drugs/results/R8_20220112/munged_stop/','',d$p2)
cols.num <- seq(3,12)
d[cols.num] <- sapply(d[cols.num],as.numeric)

dd <- d

dd$sig <- 0
dd$sig[dd$p < 0.05] <- 1

dd$p1 <- sub(".af_0.01.info_0.6.gz.rsid.gz.sumstats.gz","",dd$p1)
dd$p2 <- sub(".af_0.01.info_0.6.gz.rsid.gz.sumstats.gz","",dd$p2)
dd$stop <- substr(dd$p2, nchar(dd$p2), nchar(dd$p2))
dd$p1 <- sub("_.*", "", dd$p1)
dd$p2 <- sub("_.*", "", dd$p2)

dd$p1[dd$p1=="statins"] <- "Statins"
dd$p1[dd$p1=="blood"] <- "BP Medications"
dd$p1[dd$p1=="clopi"] <- "Clopidogrel+Dipyridamol"
dd$p1[dd$p1=="breast"] <- "Beast Cancer Medications"
dd$p1[dd$p1=="doac"] <- "DOAC"
dd$p1[dd$p1=="glauc"] <- "Glaucoma Medications"

dd$p2[dd$p2=="statins"] <- "Statins"
dd$p2[dd$p2=="blood"] <- "BP Medications"
dd$p2[dd$p2=="clopi"] <- "Clopidogrel+Dipyridamol"
dd$p2[dd$p2=="breast"] <- "Beast Cancer Medications"
dd$p2[dd$p2=="doac"] <- "DOAC"
dd$p2[dd$p2=="glauc"] <- "Glaucoma Medications"

dd$p1 <- factor(dd$p1, levels=sort(unique(dd$p1), decreasing = F))

dd$p2 <- paste0(dd$p2,"_stop_",dd$stop)
dd$p2 <- factor(dd$p2, levels=sort(unique(dd$p2), decreasing = T))

library(dplyr)
dd <- dd %>% 
  mutate(lower = rg-1.96*se,
         upper = rg+1.96*se,
         lower_cut = ifelse(lower < -1.5, -1.5, lower),
         upper_cut = ifelse(upper > 1.5, 1.5, upper))

library(ggplot2)

png("/home/ivm/drugs/results/R8_20220112/RGS_adh_stop_R8.png", width = 14, height = 11, units = "in", res = 300)
ggplot(dd, aes(x=rg, y=p2, colour = factor(sig) )) +
  geom_point() +
  geom_errorbarh(aes(xmin=lower_cut, xmax=upper_cut), height=0) +
  # geom_text(aes(label=p_label), color = "red", size = 5) +
  geom_vline(aes(xintercept = 0), alpha=.5) +
  theme_minimal() +
  facet_grid(.~factor(p1), scales = "free_y", space = "free_y", switch = "y") +
  ylab("") +
  xlab("Genetic correlation with adherence") +
  xlim(c(-1.6,1.6)) +
  scale_colour_manual(values=c("#999999", "dodgerblue3"), name = "", labels = c("", "Significant")) +
  theme(#legend.position = "none",
        text = element_text(size=17),
        #strip.placement = "outside",
        #strip.background =element_rect(fill=alpha('dodgerblue4', 0.5)),
        strip.text.y.left = element_text(angle = 0, hjust = 1))
dev.off()

