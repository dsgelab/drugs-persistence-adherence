rm(list = ls())

library(data.table)
library(ggplot2)
library(dplyr)

setwd('/home/ivm/drugs/results/R10_GWAS/persistence/sumstats/')

stats <- system('ls *rsid.AF0.01.gz', intern = T)

df <- NULL

for (s in stats){
  file <- s
  print(paste("reading file:", file))
  
  d <- gsub(".gz.rsid.AF0.01.gz","",s)
  
  drug_names <- data.frame(d = c("blood_pressure", "breast_canc", "clopi_dipy", "doac", "statins"),
                           name = c("Blood pressure medications", "Breast cancer medications", "Antiplatelets", "Anticoagulants", "Statins"),
                           Nca = c("126,126", "9,760", "12,618", "2,491", "111,309"),
                           Nco = c("8,489", "144", "4,145", "312", "5,130"))
  
  tmp1 <- fread(file, header=T)
  
  tmp <- tmp1 %>%
    filter(is.finite(pval), pval > 0, pval <= 0.01) %>%
    mutate(chrom=gsub('X', '23', `#chrom`),
           pos=as.integer(pos),
           pval_t=-log10(pval),
           odd=as.numeric(chrom) %% 2,
           chromnum=as.numeric(chrom),
           drug_n=paste0(drug_names$name[drug_names$d == d], " - ",drug_names$Nca[drug_names$d == d], " cases, ",drug_names$Nco[drug_names$d == d], " controls")) %>% 
    select(chromnum, pos, pval_t, drug_n, odd)
  
 tmp$variant <- as.character(1:nrow(tmp))
 
 df <- bind_rows(df, tmp)
}

posmin <- tapply(df$pos,df$chromnum, min)
posmax <- tapply(df$pos,df$chromnum, max)
posshift <- head(c(0,cumsum(as.numeric(posmax))),-1)
names(posshift) <- names(posmin)

for (k in unique(df$chrom))
{
  df$pos_new[df$chrom==k] <-  df$pos[df$chrom==k] + posshift[names(posshift) == k]
}

# dfmsplit <- split(df, df$chrom)
# xbreaks <- sapply(dfmsplit,function(x) x$pos_new[length(x$pos_new)/2])

xbreaks <- df %>% 
  group_by(chromnum) %>% 
  summarise(mean_pos = mean(pos_new)) %>% 
  pull(mean_pos)
names(xbreaks) <- seq(1:23)

xlabels <- names(xbreaks)
xlabels[seq(14,22,2)] <- "" 

ymax <- 8.5

df <- df %>% filter(!startsWith(drug_n, "Antiplatelets"))


png("/home/ivm/drugs/results/manuscript/R10_manhattan_persistence.png", width=1600, height=1800, res = 300)
ggplot(df, aes(x = pos_new,y = pval_t)) +
  geom_point(aes(colour=as.factor(odd)), size = .05) +
  scale_x_continuous(breaks = xbreaks, labels = xlabels, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(-log10(0.01),ymax), breaks=c(seq(2,round(ymax),length.out=7)),labels=c(round(seq(2,round(ymax),length.out=7),0)))+
  expand_limits(x = 23.3) +
  guides(colour = FALSE, alpha=FALSE, size=FALSE, fill=FALSE) +
  labs(x = "chromosome", y = expression(-log[10](italic(P)))) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +
  geom_hline(aes(yintercept= -log10(1e-08)),colour = "red", lwd=0.3, linetype = "dotted")  + 
  scale_colour_manual(values = c("deepskyblue4","deepskyblue2")) + 
  facet_wrap(~drug_n, ncol = 1) +
  ggtitle("a.") +
  theme(plot.title.position = "plot")
dev.off()
