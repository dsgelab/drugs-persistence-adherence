require(data.table)
require(dplyr)


# Exctract individual level purchase trajectories

getTrajectories <- function(dat,ATCs) {
  purch <- dat
  # TODO: impute n_packages where 0 (comes like that from the registries). Just set it to 1 now
  purch$CODE4[purch$CODE4==0] <- 1
  df <- purch %>%
    mutate(APPROX_EVENT_DAY = as.Date(APPROX_EVENT_DAY)) %>%
    filter(grepl(ATCs, CODE1)) %>%
    group_by(FINNGENID) %>%
    # remove rows duplicated because they have different INDEX values in the registries
    select(-INDEX) %>%
    distinct() %>%
    # filter trajectories with at least for purchases
    filter(n()>=4) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(pkoko_num))) %>%
    # for each individual, order by event_age
    arrange(FINNGENID, EVENT_AGE)  %>%
    # calculate n pills for each purchse
    mutate(n_pills = pkoko_num*CODE4,
           # compute difference between purchases
           days_next_purch = round( (lead(EVENT_AGE) - EVENT_AGE)*365.25)) %>%
    # TODO: deal with duplicated where days_next is 0. just remove them for now
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
      # exclude gap >=1.5 years
      days_next_purch = case_when(days_next_purch >= 365*1.5 ~ NA_real_,
                                  TRUE ~ days_next_purch),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      # set to NA last purchase and 'gap' purchases (where days_nest is.na)
      n_pills = case_when(is.na(days_next_purch) ~ NA_real_,
                          lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      pills_norm = n_pills/days_next_purch,
      days_norm = days_next_purch/n_pills)
  return(df)
}


summarizeTrajectories <- function(traj) {
  t <- traj %>%
    group_by(FINNGENID) %>%
    summarise(tot_pills = sum(n_pills, na.rm = T),
              tot_days = sum(days_next_purch, na.rm = T),
              mean_days = mean(days_next_purch, na.rm = T),
              sd_days = sd(days_next_purch, na.rm = T),
              mean_days_norm = mean(pills_norm, na.rm = T),
              sd_days_norm = sd(days_norm, na.rm = T),
              mean_pills_norm = mean(pills_norm, na.rm = T),
              sd_pills_norm = sd(pills_norm, na.rm = T),
              tot_purch = length(which(!is.na(n_pills))),
              age_first_purch = first(EVENT_AGE)) %>%
    mutate(adherence = tot_pills/tot_days) %>%
    # filter out trajectories which might have NA measures (e.g. if there is one purchase left)
    filter_all(all_vars(!is.na(.))) %>%
    # filter out outliers
    filter(adherence <= 1.1) %>%
    mutate(age_bin = cut(age_first_purch, 5))
}


# Get age of onset of the first of related endpoints
getAgeFirstEndpoint <- function(dat,ids,endpoints) {
  df <- dat
  # Discard endpoints not present in the file
  endpoints <- intersect(endpoints, colnames(df))
  # Attach _age to get the age of onset columns
  endpoints_age <- paste0(endpoints, '_AGE')
  # Select people with onset of one of the endpoints and respective age of onset
  d <- df %>%
    select(FINNGENID, all_of(endpoints), all_of(endpoints_age)) %>%
    filter_at(vars(-FINNGENID), any_vars(.==1)) %>%
    select(FINNGENID, all_of(endpoints_age))
  # For each row, the min accross columns is the age of onset of the first of the events
  min_age <- apply(d[,2:length(endpoints)], 1, min)
  d <- d %>%
    mutate(age_first_ev = min_age) %>%
    select(FINNGENID, age_first_ev)
  return(d)
}