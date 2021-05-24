rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

dir_name <- 'drugs/data/'
file_name <- 'R6_bp_w_gaps.txt' 
df <- fread(paste0(dir_name, file_name))

# Get rows where days_next_purch==0 and the following one
dup_idx <- which(df$days_next_purch == 0 & !is.na(df$days_next_purch))
length(dup_idx)/nrow(df)*100
# 0.24%
dup_idx <- sort(c(dup_idx, dup_idx+1))
dup <- df[dup_idx,]
dup$pid <- seq(1,nrow(dup))

# Case 1: same age_event, same approx_day, different VNRs and strengths.
# gap to the next purchase suggest patient is taking 2 pills
# 6116 purchases (0.23%)
dup1 <- dup %>%
  filter(EVENT_AGE == lead(EVENT_AGE) & APPROX_EVENT_DAY == lead(APPROX_EVENT_DAY) & vnr != lead(vnr))
plus.one <- dup1$pid+1
dup1 <- rbind(dup1, dup[dup$pid %in% plus.one,])
dup1 <- dup1[order(dup1$pid),]

# Case 2: same event_age, different approx_event_day (we cannot resolve gap<3.6 days having age with 2 decimals), same VNR
# 292 (0.01%)
dup2 <- dup %>%
  filter(EVENT_AGE == lead(EVENT_AGE) & APPROX_EVENT_DAY != lead(APPROX_EVENT_DAY) & vnr == lead(vnr))
plus.one <- dup2$pid+1
dup2 <- rbind(dup2, dup[dup$pid %in% plus.one,])

# Case 3: same event_age, different approx_event_day (we cannot resolve gap<3.6 days having age with 2 decimals), differnt VNR
# 367 (0.01%)
dup3 <- dup %>%
  filter(EVENT_AGE == lead(EVENT_AGE) & APPROX_EVENT_DAY != lead(APPROX_EVENT_DAY) & vnr != lead(vnr))
plus.one <- dup3$pid+1
dup3 <- rbind(dup3, dup[dup$pid %in% plus.one,])

# N pills per purchase
summary(df$n_pills)
n_pills <- df$n_pills[!is.na(df$n_pills)]
summary(n_pills)

pdf(paste0(dir_name,file_name,'.hist.pills.pdf'), 5, 3)
hist(n_pills, main=file_name)
dev.off()

# Days without pills
gaps <- df$days_next_purch - df$n_pills
gaps <- gaps[gaps>0 & gaps<601 & !is.na(gaps)]

br <- seq(0,600,by=50)
pdf(paste0(dir_name,file_name,'.hist.gaps.pdf'), 5, 3)
hist(gaps, breaks=br, include.lowest=TRUE, main=file_name)
dev.off()
