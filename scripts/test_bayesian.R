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

# Complete trajectories
st <- getTrajectories(purch,'^C10AA')

# Summarised trajectories
stt <- summarizeTrajectories(st)

# Age first related event
age_first <- getAgeFirstEndpoint(ep, stt$FINNGENID, ep_chronic)

stt <- stt %>%
  left_join(age_first) %>%
  mutate(chronic = case_when(age_first_ev <= age_first_purch ~ 1,
                             TRUE ~ 0))

# # Sample random trajectories (first 75 from yougest age_bin where it seems more likely to have 0.5 pills doses)
# fids <- unique(c(unlist(stt[stt$age_bin=='(0.744,20.1]','FINNGENID']), sample(stt$FINNGENID, 1000)))
# 
# ts <- st %>%
#   filter(FINNGENID %in% fids)

# Apply bayesian classification to each individual trajectory
st <- st %>%
  group_by(FINNGENID) %>%
  mutate(dose = classifyDailyDose(pills_norm))

table(st$dose)

# Bayesian classification of daily dose
# classifyDailyDose <- function(trajectory) {
trajectory <- st$pills_norm[st$FINNGENID == 'FG222CNUBA']
x <- trajectory[which(!is.na(trajectory))]
n <- length(x)
# Estimation starts from here
#Set parameter values
mus = c(1, 0.5) #two groups have fixed mean: the two doses we assume for that drug
sigmas = c(0.15, 0.15) #and fixed SD
n.iter = 1000 #Gibbs sampler iterations
n.burnin = 50 #Burnin period will be discarded from final results (iterations 1,...,n.burnin)

#Initialize variables for the sampler
gr = sample(c(1,2), prob = c(0.8, 0.2), size = n, replace = T) #initialize group memberships
prop = mean(gr == 1) #initialize proportion to their empirical estimate
sum.gr = matrix(0, ncol = 2, nrow = n) #col1 counts membership in group 1, col2 in group2
res.prop = rep(NA, n.iter - n.burnin) #posterior distribution of prop

#Gibbs sampler
for(ii in 1:n.iter){
  
  prop = rbeta(1, 0.5 + sum(gr == 1), 0.5 + sum(gr == 2)) #prior for prop is Beta(0.5, 0.5)
  if(ii > n.burnin) res.prop[ii-n.burnin] = prop
  
  loglk = cbind( log(prop) + dnorm(x, mus[1], sigmas[1], log = T), log(1-prop) + dnorm(x, mus[2], sigmas[2], log = T) )
  pr.1 = 1/(1 + exp(loglk[,2]-loglk[,1])) #probability of group 1 for each observation
  gr = 2 - rbinom(n, prob = pr.1, size =  1) #sample group membership, where pr.1 is probability of 1, and other outcome value is 2
  
  if(ii > n.burnin){
    sum.gr[gr == 1, 1] = sum.gr[gr == 1, 1] + 1
    sum.gr[gr == 2, 2] = sum.gr[gr == 2, 2] + 1
  }
}


p.1 <- mean(res.prop) # probability of dose 1. p.2 = 1-p.1

if (p.1 >= .6) {dose <- 1
}else if (p.1 <=.4) {dose <- 2
}else {dose <- 0
}

# return(dose)

# Classify the individual as taking dose 1 if p.1 >60%, dose 2 if p.2 < 40%, NA otherwise

layout(matrix(c(1,1,2,2), nrow = 1))
barplot(t(sum.gr)/(n.iter - n.burnin), col = c("blue","white"), main = "Membership in group 1 (in blue)" )
plot(density(x))
hist(res.prop, breaks = 40, col = "limegreen", main ="posterior of proportion of 1")
  
# }