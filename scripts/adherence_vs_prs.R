#PRS vs DRUGtraj
rm(list = ls())
library(data.table)
library(ggplot2)
library(dplyr)
#LOAD DATA and PRS and MERGE
st <- fread('/home/cordioli/drugs/data/statins_summarized.txt')
bp <- fread('/home/cordioli/drugs/data/bp_summarized.txt')
cl <- fread('/home/cordioli/drugs/data/clopi_dipy_summarized.txt')


prs_ldl <- fread("/home/cordioli/drugs/prs/ldlstatinadj_FinnGen_LDpred_GRS.txt", hea=T, data.table=F)
head(prs_ldl)
prs_ldl <- prs_ldl[,c(1, 4)]
colnames(prs_ldl) <- c("FINNGENID", "prs_ldl")

prs_sbp <- fread("/home/cordioli/drugs/prs/sbp_FinnGen_LDpred_GRS.txt", data.table = F)
head(prs_sbp)
prs_sbp <- prs_sbp[, c(1, 4)]
colnames(prs_sbp) <- c("FINNGENID", "prs_sbp")

st <- merge(st, prs_ldl, by ="FINNGENID", all.x = T)
bp <- merge(bp, prs_sbp, by = "FINNGENID", all.x = T)

summary(lm(formula = adherence_std ~ prs_ldl, data = st))

#SCALE PRS
#data$PRS_dm2 <- scale(data$prs_dm2)
#data$PRS_ldl <- scale(data$prs_ldl)
#data$PRS_sbp <- scale(data$prs_sbp)
df <- bp
df$prs <- df$prs_sbp
tenths<- quantile(df$prs, seq(0.1,0.9, 0.1))
df$tenths <- ifelse(df$prs < tenths[1], "0-10%",
                    ifelse(df$prs < tenths[2], "10-20%", 
                           ifelse(df$prs< tenths[3], "20-30%",
                                  ifelse(df$prs < tenths[4], "30-40%",
                                         ifelse(df$prs < tenths[5], "40-50%",
                                                ifelse(df$prs < tenths[6], "50-60%",
                                                       ifelse(df$prs < tenths[7], "60-70%",
                                                              ifelse(df$prs < tenths[8], "70-80%",
                                                                     ifelse(df$prs < tenths[9], "80-90%", "90-100%")))))))))

sumdf <- df %>%
  group_by(tenths) %>%
  summarise(N = n(),
            adherence_mean = mean(adherence, na.rm = T),
            adherence_sd = sd(adherence, na.rm = T)) %>%
  # Calculate standard error of the mean
  mutate(mean_se = adherence_sd/sqrt(N))

# Confidence interval multiplier for standard error
# Calculate t-statistic for confidence interval: 
# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
ciMult <- qt(.975/2 + .5, sumdf$N-1)    
sumdf <- sumdf %>%
  mutate(ci = mean_se * ciMult)

ymin <- min(sumdf[1,3]-sumdf[1,6]*2,sumdf[10,3]-sumdf[10,6]*2)
ymax <- max(sumdf[1,3]+sumdf[1,6]*2,sumdf[10,3]+sumdf[10,6]*2)


pd <- position_dodge(0.1)
ggplot(sumdf, aes(x=tenths, y=adherence_mean)) + 
  geom_errorbar(aes(ymin=adherence_mean-ci, ymax=adherence_mean+ci), colour="black", width=.2, position=pd) +
  geom_line(position=pd, size = 2) +
  geom_point(position=pd, size=3)+
  ylim(ymin,ymax)+
  xlab("PRS systolic blood pressure")+
  ylab("Mean adherence to blood pressure medications")+
  #scale_colour_hue(name="Sex",    # Legend label, use darker colors
  #                 breaks=c("Female", "Male"),
  #                 labels=c("Females", "Males"),
  #                 l=40) +   
  #ggtitle("The Effect of alcohol consumption PRS on\nweekly alcohol consumption estimate") +
  expand_limits(y=0) +                        # Expand y range
  #scale_y_continuous(breaks=0:30*5) +         # Set tick every 4
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.justification=c(1,0),
        legend.position=c(0.85,0.15))
