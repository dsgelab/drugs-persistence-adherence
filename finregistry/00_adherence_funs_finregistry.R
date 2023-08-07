require(data.table)
require(dplyr)
require(lubridate)

# Exctract individual level purchase trajectories using event date
getTrajectoriesDates <- function(purch,dose=1) {
  # ATCs as a regex string (e.g. '^C10AA' for statins) 
  df <- purch %>%
    group_by(HETU) %>%
    # filter trajectories with at least 4 purchases
    filter(n()>=4) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(PKOKO))) %>%
    # remove purchases dispensed through the ANJA system 
    filter(!ANJA %in% c('U','K','B')) %>% 
    # for each individual, order by event_age
    arrange(OSTOPV, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = round( (lead(OSTOPV) - OSTOPV)),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific dose
      n_pills = case_when(PLKM == 0 ~ (days_next_purch - (days_next_purch %% PKOKO))/dose,
                          TRUE ~ PKOKO*PLKM/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = ATC != lag(ATC),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      n_pills = case_when(lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      # identify gap >=150 days without pills
      gap = ifelse(days_next_purch > n_pills+150, 1, 0),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      gap_days = ifelse(gap == 1, days_next_purch - n_pills, 0),
      # set days to NA for gaps
      days_next_purch = case_when(gap == 1 ~ NA_real_,
                                  TRUE ~ days_next_purch),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                          TRUE ~ n_pills),
      pills_norm = n_pills/days_next_purch,
      days_norm = days_next_purch/n_pills)
  return(df)
}


getTrajectoriesDates2 <- function(purch,ATCs) {
  # ATCs as a df with daily dose of each ATC
  df <- purch %>%
    inner_join(ATCs, by = c("ATC"="atc")) %>%
    group_by(HETU) %>%
    # filter trajectories with at least four purchases
    filter(n()>=4) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(PKOKO))) %>%
    # remove purchases dispensed through the ANJA system 
    filter(!ANJA %in% c('U','K','B')) %>% 
    # for each individual, order by event_age
    arrange(OSTOPV, .by_group = T)  %>%
    # calculate n pills for each purchase, adjusting for the dose
    mutate(# compute difference between purchases
      days_next_purch = round( (lead(OSTOPV) - OSTOPV)),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific doses
      n_pills = case_when(PLKM == 0 ~ (days_next_purch - (days_next_purch %% PKOKO))/dose,
                          TRUE ~ PKOKO*PLKM/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = ATC != lag(ATC),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      n_pills = case_when(lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      # identify gap >=150 days without pills
      gap = ifelse(days_next_purch > n_pills+150, 1, 0),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      gap_days = ifelse(gap == 1, days_next_purch - n_pills, 0),
      # set days to NA for gaps
      days_next_purch = case_when(gap == 1 ~ NA_real_,
                                  TRUE ~ days_next_purch),
      # set to NA last purchase and 'gap' purchases (where days_nest is.na)
      n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                          TRUE ~ n_pills),
      pills_norm = n_pills/days_next_purch,
      days_norm = days_next_purch/n_pills)
  return(df)
}


summarizeTrajectories <- function(traj, min_age = 0, thr = 1.1) {
  t <- traj %>%
    group_by(HETU) %>%
    summarise(tot_pills = sum(n_pills, na.rm = T),
              tot_days = sum(days_next_purch, na.rm = T),
              mean_days = mean(days_next_purch, na.rm = T),
              sd_days = sd(days_next_purch, na.rm = T),
              tot_purch = length(which(!is.na(n_pills))),
              age_first_purch = first(EVENT_AGE),
              n_break = sum(gap, na.rm = T)) %>%
    mutate(adherence = tot_pills/tot_days) %>%
    # filter out trajectories which result to have NA for the previous metrics (e.g. if there is one purchase left)
    filter_all(all_vars(!is.na(.))) %>%
    # filter out outliers
    filter(adherence <= thr,
           age_first_purch >= min_age,
           tot_days >= 365) %>%
    mutate(age_bin = cut(age_first_purch, breaks = c(18, 48, 60, 80, Inf), labels = c("18-39", "40-59", "60-79", "80+"), right= FALSE)) %>%
    mutate(adherence_std = as.numeric(scale(adherence)))
}


# Extract early stopping purchases
getTrajectoriesPersistence <- function(df) {
  # Define "persistent" users:
  # - at least one year treatment
  # - without breaks
  # - no overbuyers (adherence < 1.1))
   df <- df %>% 
    filter(tot_days >= 365) %>% 
    mutate(persistent = 1) %>% 
    select(HETU, tot_pills, age_first_purch, age_bin, tot_purch, persistent)
  return(df)
}


getTrajectoriesDiscontinuation <- function(purch,cov,dose=1,min_age) {
  # Define early stopping users
  df <- purch %>%
    left_join(cov, by = c("HETU" = "FINREGISTRYID")) %>% 
    group_by(HETU) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(PKOKO))) %>%
    # for each individual, order by event_date
    arrange(OSTOPV, .by_group = T)  %>%
    # keep trajectories with only one purchase at least 2 years before end of followup (EOF = death OR emigration OR 2020-01-01)
    filter(n()==1) %>%
    mutate(END_OF_FOLLOWUP = min(death_date, emigration_date, as.Date('2020-01-01'), na.rm = T)) %>%
    filter(as.numeric(difftime(as.Date(END_OF_FOLLOWUP), as.Date(OSTOPV), units = "days"))/365.25 > 2) %>%
    # calculate n pills for each purchase
    mutate(tot_pills = PKOKO*PLKM/dose,
           persistent = 0,
           tot_purch = 1,
           persistent = 0) %>%
    mutate(age_first_purch = as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25,
           age_bin = cut(age_first_purch, breaks = c(18, 48, 60, 80, Inf), labels = c("18-39", "40-59", "60-79", "80+"), right= FALSE)) %>%
    filter(age_first_purch >= min_age) %>%
    select(HETU, tot_pills, age_first_purch, age_bin, tot_purch, persistent)
  return(df)
}


getTrajectoriesDiscontinuation2 <- function(purch,cov,ATCs,min_age) {
  # ATCs as a df with daily dose of each ATC
    df <- purch %>%
      inner_join(ATCs, by = c("ATC"="atc")) %>%
      left_join(cov, by = c("HETU" = "FINREGISTRYID")) %>% 
      group_by(HETU) %>%
      # filter trajectories with VNR info (package size) for all the events
      filter(all(!is.na(PKOKO))) %>%
      # for each individual, order by event_date
      arrange(OSTOPV, .by_group = T)  %>%
      # keep trajectories with only one purchase at least 2 years before end of followup (EOF = death OR emigration OR 2020-01-01)
      filter(n()==1) %>%
      mutate(END_OF_FOLLOWUP = min(death_date, emigration_date, as.Date('2020-01-01'), na.rm = T)) %>%
      filter(as.numeric(difftime(as.Date(END_OF_FOLLOWUP), as.Date(OSTOPV), units = "days"))/365.25 > 2) %>%
      # calculate n pills for each purchase
      mutate(tot_pills = PKOKO*PLKM/dose,
             persistent = 0,
             tot_purch = 1,
             persistent = 0) %>%
      mutate(age_first_purch = as.numeric(difftime(OSTOPV, date_of_birth, units = "days"))/365.25,
             age_bin = cut(age_first_purch, breaks = c(18, 48, 60, 80, Inf), labels = c("18-39", "40-59", "60-79", "80+"), right= FALSE)) %>%
      filter(age_first_purch >= min_age) %>%
      select(HETU, tot_pills, age_first_purch, age_bin, tot_purch, persistent)
    return(df)
}