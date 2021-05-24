rm(list=ls())

library(dplyr)
library(data.table)


dat <- fread('/home/cordioli/R5_pheno/finngen_R5_v2_detailed_longitudinal.gz')
head(dat)

dat <- subset(dat, SOURCE == "PURCH")


bp <- subset(dat, grepl('^C0[2|3|7|8|9]', CODE1))
dm2 <- subset(dat, grepl('^A10B', CODE1))
insulins <- subset(dat, grepl('^A10A', CODE1))
statins <- subset(dat, grepl('^C10AA', CODE1))

#cardiometabolic <- rbind(bp,dm2,insulins,statins)

#   Field |Register        |Register long name         |
#   ------|----------------|---------------------------|
#   SOURCE|PURCH           |Kela drug purchase register|
#   
#   Field |Code in register|Code description.         |
#   ------|---------------|--------------------------|
#   CODE1 | ATC_CODE      | ATC code                  |
#   CODE2 | SAIR 	        | Kela reimbursement code   |
#   CODE3 | VNRO	        | Product number		        |
#   CODE4 | PLKM          | Number of packages        |

# > length(which(is.na(statins$CODE2)))
# [1] 1743916
# Number of NA VRNOs
# > length(which(is.na(statins$CODE3)))
# [1] 0
# Number of NA numb. packages
# > length(which(is.na(statins$CODE4)))
# [1] 0
# Number packages = 0
# > length(which(statins$CODE4==0))
# [1] 35
# Number of different VRNOs
# > length(unique(statins$CODE3))
# [1] 447

nrow(statins)
# [1] 2202346 total purchases

statins$APPROX_EVENT_DAY <- as.Date(statins$APPROX_EVENT_DAY)

min(statins$APPROX_EVENT_DAY)
max(statins$APPROX_EVENT_DAY)
# [1] "1995-12-20"
# [1] "2019-01-15"

unique(statins$CODE1)
st_codes <- c("C10AA01", "C10AA02", "C10AA03", "C10AA04", "C10AA05", "C10AA06", "C10AA07")
st_names <- c("Simvastatin", "Lovastatin", "Pravastatin", "Fluvastatin", "Atorvastatin", "Cerivastatin", "Rosuvastatin")

statins$NAME <- ""
for (i in 1:length(st_codes)) {
  statins$NAME[statins$CODE1 == st_codes[i]] <- st_names[i]
}

statins$NAME <- factor(statins$NAME, levels = st_names)

length(unique(statins$FINNGENID))
# [1] 78291

purch_ind <- data.frame(table(statins$FINNGENID))

# Basic histograms
library(ggplot2)
ggplot(purch_ind, aes(x=Freq)) + 
  geom_histogram(bins=50) +
  ggtitle('Purchases of all statins (C10AA[01-07]) per individual') +
  theme_minimal()

ggplot(statins, aes(x=NAME)) + 
  geom_histogram(stat = "count") +
  ggtitle('Purchases of different type statins') +
  theme_minimal()


# Time-window: 2005-2010
# Include people already taking statins in 2005
# Exclude who die in the time window
# Exclude specific complications

d <- fread('/home/cordioli/R5_pheno/finngen_R5_V2_endpoint.gz')

died <- d %>%
  filter(between(DEATH_YEAR, 2005, 2010)) %>%
  pull(FINNGENID)

users_04 <- statins %>%
  filter(between(APPROX_EVENT_DAY, "2004-07-01", "2004-12-31")) %>%
  count(FINNGENID) %>%
  filter(n > 1) %>%
  pull(FINNGENID)

s_05_10 <- statins %>%
  filter(between(APPROX_EVENT_DAY, "2005-01-01", "2010-12-31") &
           !FINNGENID %in% died &
           FINNGENID %in% users_04) %>%
  group_by(FINNGENID) %>% 
  summarise(m = min(APPROX_EVENT_DAY)) %>%
  filter(m < "2005-03-01") %>%
  pull(FINNGENID)

dd <- statins %>%
  filter(FINNGENID %in% s_05_10)


# length(unique(s_05_10$FINNGENID))
# [1] 23033

# plot some trajectories

df <- dd %>%
  filter(FINNGENID %in% s_05_10[1:150])

library(ggplot2)

ggplot(df,aes(x=EVENT_AGE,y=factor(FINNGENID),colour=FINNGENID,group=FINNGENID)) + 
  geom_line() +
  theme_minimal()

ggplot(df,aes(x=APPROX_EVENT_DAY,y=factor(FINNGENID),colour=FINNGENID,group=FINNGENID)) + 
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none")
