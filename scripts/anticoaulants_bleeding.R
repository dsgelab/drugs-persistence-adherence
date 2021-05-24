rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

setwd('/home/cordioli/')

to.keep <-  c("FINNGENID", "sd_days", "mean_days_norm", "sd_days_norm", "tot_purch", "age_first_purch", "adherence", "age_bin", 
              "chronic", "SEX", "HEIGHT", "WEIGHT", "SMOKE2", "SMOKE3", "SMOKE5", "SMOKE_AGE","BMI")

do <- fread('drugs/data/R6_doac_summarized_150.txt') %>%
  select(all_of(to.keep)) %>%
  mutate(drug="DOAC")

do_traj <- fread('drugs/data/R6_doac_gap_150.txt')

# Mount finngen data bucket
system("gcsfuse --only-dir finngen_R6/phenotype_2.0/data/ --file-mode 444 finngen-production-library-red /home/cordioli/mount_fg/")
ep <- fread('/home/cordioli/mount_fg/finngen_R6_v2_endpoint.gz') %>%
  filter(FINNGENID %in% do$FINNGENID)
# Unmount finngen data
system("fusermount -u /home/cordioli/mount_fg")

d <- do_traj %>%
  group_by(FINNGENID) %>%
  summarize(st = min(EVENT_AGE), end = max(EVENT_AGE), end_plus_1m = round(max(EVENT_AGE)+(30/3.65)*0.01, 2)) %>%
  filter(FINNGENID %in% do$FINNGENID)

bleeding <- c("D3_HAEMORRHAGCIRGUANTICO","H7_CONJUHAEMOR","H7_RETINAHAEMORR","H7_VITRHAEMORR","H7_CHORHAEMORRHAGE","ST19_EPIDU_HAEMORRHAGE",
              "I9_INTRACRA","I9_OTHINTRACRA","ST19_EPIDU_HAEMORRHAGE","ST19_TRAUMAT_SUBDU_HAEMORRHAGE","ST19_TRAUMAT_SUBAR_HAEMORRHAGE","
              R18_HAEMORRHAGE_RESPI_PASSA","K11_GIBLEEDING")

getAgeEndpoints <- function(dat,ids,endpoints) {
  df <- dat
  # Discard endpoints not present in the file
  endpoints <- intersect(endpoints, colnames(df))
  # Attach _age to get the age of onset columns
  endpoints_age <- paste0(endpoints, '_AGE')
  # Select people with onset of one of the endpoints and respective age of onset
  d <- df %>%
    select(FINNGENID, FU_END_AGE, all_of(endpoints), all_of(endpoints_age)) %>%
    filter_at(vars(-FINNGENID), any_vars(.==1)) %>%
    select(FINNGENID, FU_END_AGE, all_of(endpoints_age)) %>%
    filter(FINNGENID %in% ids)
  
  d <- d %>%
    mutate_at(.vars = vars(endpoints_age),
              .funs = funs(case_when(. == d$FU_END_AGE ~ NA_real_,
                                     TRUE ~ .))) %>%
    select(-FU_END_AGE)
  
  return(d)
}

d_age <- getAgeEndpoints(ep, d$FINNGENID, bleeding)

d_age_t <- melt(d_age)

d_age_t <- d_age_t %>%
  filter(!is.na(value)) %>%
  left_join(d)

d_age_t$bleeding <- 0
d_age_t$bleeding[between(d_age_t$value, d_age_t$st, d_age_t$end_plus_1m)] <- 1

sum(d_age_t$bleeding)
