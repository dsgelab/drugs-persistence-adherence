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
    arrange(EVENT_AGE, .by_group = T)  %>%
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


summarizeTrajectories <- function(traj, thr = 1.1) {
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
    filter(adherence <= thr) %>%
    mutate(age_bin = cut(age_first_purch, 5)) %>%
    mutate(adherence_norm = as.numeric(scale(adherence)))
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


# Bayesian classification of daily dose
classifyDailyDose <- function(trajectory) {
  n <- length(trajectory)
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
  
  return(dose)
  
  # Classify the individual as taking dose 1 if p.1 >60%, dose 2 if p.2 < 40%, NA otherwise
  
  
  #layout(matrix(c(1,1,2,2), nrow = 1))
  # barplot(t(sum.gr)/(n.iter - n.burnin), col = c("blue","white"), main = "Membership in group 1 (in blue)" )
  #plot(density(x))
  #hist(res.prop, breaks = 40, col = "limegreen", main ="posterior of proportion of 1")
  
}