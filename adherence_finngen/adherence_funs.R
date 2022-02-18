require(data.table)
require(dplyr)
require(lubridate)


# Exctract individual level purchase trajectories
getTrajectories <- function(purch,ATCs,dose=1) {
  # ATCs as a regex string (e.g. '^C10AA' for statins) 
  df <- purch %>%
    filter(grepl(ATCs, CODE1)) %>%
    group_by(FINNGENID) %>%
    # remove rows duplicated because they have different INDEX values in the registries
    select(-INDEX, -SOURCE, -ATC) %>%
    distinct() %>%
    # filter trajectories with at least 4 purchases
    filter(n()>=4) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(pkoko_num))) %>%
    # for each individual, order by event_age
    arrange(EVENT_AGE, .by_group = T)  %>%
    mutate(# compute difference between purchases
           days_next_purch = round( (lead(EVENT_AGE) - EVENT_AGE)*365.25),
           # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
           # adjusting specific dose
           n_pills = case_when(CODE4 == 0 ~ (days_next_purch - (days_next_purch %% pkoko_num))/dose,
                               TRUE ~ pkoko_num*CODE4/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
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


# Exctract individual level purchase trajectories using event date
getTrajectoriesDates <- function(purch,ATCs,dose=1) {
  # ATCs as a regex string (e.g. '^C10AA' for statins) 
  df <- purch %>%
    filter(grepl(ATCs, CODE1)) %>%
    group_by(FINNGENID) %>%
    # remove rows duplicated because they have different INDEX values in the registries
    select(-INDEX, -SOURCE, -ATC) %>%
    distinct() %>%
    # filter trajectories with at least 4 purchases
    filter(n()>=4) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(pkoko_num))) %>%
    # for each individual, order by event_age
    arrange(APPROX_EVENT_DAY, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = round( (lead(APPROX_EVENT_DAY) - APPROX_EVENT_DAY)),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific dose
      n_pills = case_when(CODE4 == 0 ~ (days_next_purch - (days_next_purch %% pkoko_num))/dose,
                          TRUE ~ pkoko_num*CODE4/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
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


# Extract individual level purchase trajectories to look at early stopping
getTrajectoriesStopping <- function(purch,ATCs,dose=1) {
  # ATCs as a regex string (e.g. '^C10AA' for statins) 
  df <- purch %>%
    filter(grepl(ATCs, CODE1)) %>%
    group_by(FINNGENID) %>%
    # remove rows duplicated because they have different INDEX values in the registries
    select(-INDEX, -SOURCE, -ATC) %>%
    distinct() %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(pkoko_num))) %>%
    # for each individual, order by event_age
    arrange(APPROX_EVENT_DAY, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = round( (lead(APPROX_EVENT_DAY) - APPROX_EVENT_DAY)),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific dose
      n_pills = case_when(CODE4 == 0 ~ (days_next_purch - (days_next_purch %% pkoko_num))/dose,
                          TRUE ~ pkoko_num*CODE4/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      n_pills = case_when(lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      # identify gap >=150 days without pills
      gap = ifelse(days_next_purch > n_pills+150, 1, 0),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      gap_days = ifelse(gap == 1, days_next_purch - n_pills, 0))
  return(df)
}


getTrajectories2 <- function(purch,ATCs) {
  # ATCs as a df with daily dose of each ATC
  df <- purch %>%
    inner_join(ATCs, by = c("CODE1"="atc")) %>%
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
    # calculate n pills for each purchase, adjusting for the dose
    mutate(# compute difference between purchases
           days_next_purch = round( (lead(EVENT_AGE) - EVENT_AGE)*365.25),
           # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
           # adjusting specific doses
           n_pills = case_when(CODE4 == 0 ~ (days_next_purch - (days_next_purch %% pkoko_num))/dose,
                               TRUE ~ pkoko_num*CODE4/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
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


getTrajectoriesDates2 <- function(purch,ATCs) {
  # ATCs as a df with daily dose of each ATC
  df <- purch %>%
    inner_join(ATCs, by = c("CODE1"="atc")) %>%
    group_by(FINNGENID) %>%
    # remove rows duplicated because they have different INDEX values in the registries
    select(-INDEX) %>%
    distinct() %>%
    # filter trajectories with at least four purchases
    filter(n()>=4) %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(pkoko_num))) %>%
    # for each individual, order by event_age
    arrange(APPROX_EVENT_DAY, .by_group = T)  %>%
    # calculate n pills for each purchase, adjusting for the dose
    mutate(# compute difference between purchases
      days_next_purch = round( (lead(APPROX_EVENT_DAY) - APPROX_EVENT_DAY)),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific doses
      n_pills = case_when(CODE4 == 0 ~ (days_next_purch - (days_next_purch %% pkoko_num))/dose,
                          TRUE ~ pkoko_num*CODE4/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
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


getTrajectoriesStopping2 <- function(purch,ATCs) {
  # ATCs as a df with daily dose of each ATC
  df <- purch %>%
    inner_join(ATCs, by = c("CODE1"="atc")) %>%
    group_by(FINNGENID) %>%
    # remove rows duplicated because they have different INDEX values in the registries
    select(-INDEX) %>%
    distinct() %>%
    # filter trajectories with VNR info (package size) for all the events
    filter(all(!is.na(pkoko_num))) %>%
    # for each individual, order by event_age
    arrange(APPROX_EVENT_DAY, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = round( (lead(APPROX_EVENT_DAY) - APPROX_EVENT_DAY)),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific doses
      n_pills = case_when(CODE4 == 0 ~ (days_next_purch - (days_next_purch %% pkoko_num))/dose,
                          TRUE ~ pkoko_num*CODE4/dose)) %>%
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = CODE1 != lag(CODE1),
      # adjust n pills: if I change formulation after e.g. 20 days and I've bought 100 pills, I count only 20 pills for that purch
      n_pills = case_when(lead(change_type) == TRUE & days_next_purch <= n_pills ~ days_next_purch,
                          TRUE ~ n_pills),
      # identify gap >=150 days without pills
      gap = ifelse(days_next_purch > n_pills+150, 1, 0),
      # set n_pills to NA for last purchase and 'gap' purchases (where days_next is.na)
      gap_days = ifelse(gap == 1, days_next_purch - n_pills, 0))
  return(df)
}

# Exctract individual level purchase trajectories from UKB data
getTrajectoriesUKB <- function(scripts,drug_names_list,dose=1) {
  # BNFs as a regex string (e.g. '^C10AA' for statins) 
  df <- scripts %>%
    filter(drug_name %in% drug_names_list) %>%
    group_by(eid) %>%
    # remove duplicated
    distinct() %>%
    # filter trajectories with at least 4 purchases, and with the same provider for all
    filter(n()>=4,
           n_distinct(data_provider) == 1,
           issue_date >= date_recruitment) %>%
    # keep only trajectories with quantities for all the events
    filter(all(!is.na(quantity_num))) %>%
    # for each individual, order by event_age
    arrange(issue_date, .by_group = T)  %>%
    mutate(# compute difference between purchases
      days_next_purch = as.numeric(difftime(lead(issue_date), issue_date, units = "d")),
      # calculate n pills for each purchase, imputing when n pkg is 0 (e.g. 238 days, 100 pills: 238 %% 100 = 38, n_pills = 200)
      # adjusting specific dose
      n_pills = as.numeric(as.numeric(quantity_num)/dose)) %>%
    # TODO: deal with duplicated where days_next is 0. just remove them for now
    filter(days_next_purch != 0 | is.na(days_next_purch)) %>%
    mutate(
      change_type = drug_name != lag(drug_name),
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


summarizeTrajectories <- function(traj, min_age = 0, thr = 1.1) {
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
              age_first_purch = first(EVENT_AGE),
              n_break = sum(gap, na.rm = T)) %>%
    mutate(adherence = tot_pills/tot_days) %>%
    # filter out trajectories which result to have NA for the previous metrics (e.g. if there is one purchase left)
    filter_all(all_vars(!is.na(.))) %>%
    # filter out outliers
    filter(adherence <= thr,
           age_first_purch >= min_age) %>%
    mutate(age_bin = cut(age_first_purch, 5)) %>%
    mutate(adherence_std = as.numeric(scale(adherence)))
}


summarizeTrajectoriesUKB <- function(traj, thr = 1.1) {
  t <- traj %>%
    group_by(eid) %>%
    summarise(tot_pills = sum(n_pills, na.rm = T),
              tot_days = sum(days_next_purch, na.rm = T),
              mean_days = mean(days_next_purch, na.rm = T),
              sd_days = sd(days_next_purch, na.rm = T),
              mean_days_norm = mean(pills_norm, na.rm = T),
              sd_days_norm = sd(days_norm, na.rm = T),
              mean_pills_norm = mean(pills_norm, na.rm = T),
              sd_pills_norm = sd(pills_norm, na.rm = T),
              tot_purch = length(which(!is.na(n_pills))),
              age_first_purch = first(lubridate::time_length(difftime(first(issue_date), date_recruitment), "years") + age_recruitment) ) %>%
    mutate(adherence = tot_pills/tot_days) %>%
    # filter out trajectories which result to have NA for the previous metrics (e.g. if there is one purchase left)
    filter_all(all_vars(!is.na(.))) %>%
    # filter out outliers
    filter(adherence <= thr) %>%
    mutate(age_bin = cut(age_first_purch, 5)) %>%
    mutate(adherence_std = as.numeric(scale(adherence)))
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
    mutate(age_first_ev = as.numeric(min_age)) %>%
    select(FINNGENID, age_first_ev)
  return(d)
}




# Bayesian classification of daily dose
classifyDailyDose <- function(trajectory) {
  x <- trajectory
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
  
}


# # # Plots
plotAdherenceDensity <- function(df, drug) {
  
  require(grid)
  require(gridExtra)
  require(cowplot)

  p1 <- ggplot(data = df, aes(y = adherence)) +
    geom_density(alpha = .8) +
    geom_boxplot(width = .5, alpha = 0.8, outlier.shape = NA) +
    coord_flip() +
    xlab("Density") +
    ylab("Adherence") +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.1)) +
    ggtitle(paste0('Adherence to ', drug)) +
    theme_minimal()
  
  # Create labels for type of prevention
  df$Prevention <- ifelse(df$chronic == 0, "Primary", "Secodary")
  
  p2 <- ggplot(data = df, aes(y = adherence, fill = factor(Prevention))) +
    geom_density(alpha = .8) +
    geom_boxplot(width = .5, alpha = 0.8, outlier.shape = NA) +
    coord_flip() +
    xlab("Density") +
    ylab("Adherence") +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.1)) +
    ggtitle(paste0('Adherence to ', drug, ' by type of prevention')) +
    facet_grid(Prevention~.) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr", rel_heights = c(0.35, 0.65)))
}


plotAdherenceByAge <- function(df, drug) {
  
  require(cowplot)
  require(grid)
  require(gridExtra)

  palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00")
  p1 <- ggplot(df, aes(y=adherence, fill = age_bin)) + 
    geom_density(alpha = .8) +
    geom_boxplot(width = .5, alpha = 0.8, outlier.shape = NA) +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.1)) +
    coord_flip() +
    ggtitle(paste0('Adherence to ', drug,' by age at initiation')) +
    xlab("Density") +
    ylab("Adherence") +
    facet_grid(age_bin~.) +
    theme_minimal() +
    scale_fill_manual(values=palette) +
    theme(legend.position = "none")
  
  p2 <- ggplot(df, aes(x=age_bin, fill=age_bin)) +
    geom_bar(stat = 'count') +
    ggtitle(paste0('Number of ', drug,' users per age bin')) +
    labs(fill = "Age at initiation") +
    ylab("N users") +
    xlab("Age at initiation") +
    scale_fill_manual(values=palette) +
    theme_minimal()
  
  return(plot_grid(p2, p1, align = "h", ncol = 2, axis = "b", rel_widths = c(2/3, 1/3)))
}




plotAdherenceByAgeViolin <- function(df,drug) {
  require(ggplot2)
  require(cowplot)
  require(grid)
  require(gridExtra)
  
  labs <- table(factor(df$age_bin)) %>%
    as.numeric()
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m - sd(x)
    ymax <- m + sd(x)
    return(c(y=m, ymin=ymin, ymax=ymax))}
  
  p <-ggplot(traj, aes(x = actor(age_bin), 
                       y = adherence,
                       fill = factor(age_bin))) +
    geom_violin(trim=F) +
    scale_y_continuous(limits = c(0, 1.1)) +
    geom_hline(aes(yintercept = .8), 
               color = "black", 
               alpha = .75) +
    stat_summary(fun.data = data_summary, 
                 shape = 18, 
                 size =.5,
                 geom = "pointrange", 
                 color = "red") +
    scale_fill_brewer(palette = "Blues", 
                      name = "Group size",
                      labels = labs) + 
    ggtitle(paste0("Adherence to ", drug, "by age at initiation")) +
    xlab("Age bin") +
    ylab("Adherence") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 15),
          legend.text =  element_text(size = 10))
}
