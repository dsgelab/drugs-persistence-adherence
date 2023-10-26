rm(list=ls())

library(data.table)
library(dplyr)
library(broom)
library(gridExtra)
library(ggplot2)
library(forcats)
library(gtools)

# # # Load adherence and persistence phenotypes

# Adherence
adh <- fread("/home/ivm/drugs/data/R10_cov_pheno_adherence.txt") %>% mutate(FINNGENID = IID)

dd <- fread("/home/ivm/drugs/data/R10_blood_pressure_summarized.txt") %>%
  mutate(blood_pressure_adh = adherence) %>%
  select(FINNGENID, blood_pressure_adh)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_breast_cancer_summarized.txt") %>%
  mutate(breast_canc_adh = adherence) %>%
  select(FINNGENID, breast_canc_adh, tamoxifen_prop=prop_atc, tamoxifen_atc=atc)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_clopi_dipy_summarized.txt") %>%
  mutate(clopi_dipy_adh = adherence) %>%
  select(FINNGENID, clopi_dipy_adh, clopidogrel_prop=prop_atc, clopidogrel_atc=atc)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_doac_summarized.txt") %>% 
  mutate(doac_adh = adherence) %>%
  select(FINNGENID,doac_adh)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_statins_summarized.txt") %>% 
  mutate(statins_adh = adherence) %>%
  select(FINNGENID, statins_adh)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_statins_rosu_summarized.txt") %>% 
  mutate(statins_adh = adherence) %>%
  select(FINNGENID, rosuvastatin_prop = prop_atc, rosuvastatin_atc = atc)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_statins_flu_summarized.txt") %>% 
  mutate(statins_adh = adherence) %>%
  select(FINNGENID, fluvastatin_prop = prop_atc, fluvastatin_atc = atc)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

adh <- adh %>% mutate(age_first_statins = scale(age_first_statins),                                
                      age_first_blood_pressure = scale(age_first_blood_pressure),                  
                      age_first_clopi_dipy = scale(age_first_clopi_dipy),                         
                      age_first_breast_canc = scale(age_first_breast_canc),                           
                      age_first_doac = scale(age_first_doac),
                      PC1 = scale(PC1),
                      PC2 = scale(PC2),
                      PC3 = scale(PC3),
                      PC4 = scale(PC4),
                      PC5 = scale(PC5),
                      PC6 = scale(PC6),
                      PC7 = scale(PC7),
                      PC8 = scale(PC8),
                      PC9 = scale(PC9),
                      PC10 = scale(PC10))

# Persistence
per <- fread("/home/ivm/drugs/data/R10_cov_pheno_persistence.txt") %>% mutate(FINNGENID = IID, breast_canc = breast_cancer)

dd <- fread("/home/ivm/drugs/data/R10_breast_cancer_persistence.txt") %>%
  select(FINNGENID, tamoxifen_prop=prop_atc, tamoxifen_atc=atc)
per <- merge(x = per, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_clopi_dipy_persistence.txt") %>% 
  select(FINNGENID, clopidogrel_prop=prop_atc, clopidogrel_atc=atc)
per <- merge(x = per, y = dd , by = "FINNGENID", all.x=TRUE)

# Add persistence for fluvastatin and rosuvastatin specifically
pp <- fread("/home/ivm/drugs/data/R10_statins_flu_summarized.txt") %>% 
  filter(prop_atc == 1) %>% 
  mutate(fluvastatins = 1) %>% 
  select(FINNGENID, fluvastatins, age_first_fluvastatins = age_first_purch)
ss <- fread("/home/ivm/drugs/data/R10_statins_persistence.txt") %>% 
  filter(persistent == 0,
         atc == "C10AA04")%>% 
  select(FINNGENID, fluvastatins = persistent, age_first_fluvastatins = age_first_purch)
dd <- bind_rows(ss, pp)
per <- merge(x = per, y = dd , by = "FINNGENID", all.x=TRUE)

pp <- fread("/home/ivm/drugs/data/R10_statins_rosu_summarized.txt") %>% 
  filter(prop_atc == 1) %>% 
  mutate(rosuvastatins = 1) %>% 
  select(FINNGENID, rosuvastatins, age_first_rosuvastatins = age_first_purch)
ss <- fread("/home/ivm/drugs/data/R10_statins_persistence.txt") %>% 
  filter(persistent == 0,
         atc == "C10AA07")%>% 
  select(FINNGENID, rosuvastatins = persistent, age_first_rosuvastatins = age_first_purch)
dd <- bind_rows(ss, pp)
per <- merge(x = per, y = dd , by = "FINNGENID", all.x=TRUE)

# Scale age and PCs
per <- per %>% mutate(age_first_statins = scale(age_first_statins),
                      age_first_fluvastatins = scale(age_first_fluvastatins),
                      age_first_rosuvastatins = scale(age_first_rosuvastatins),
                      age_first_blood_pressure = scale(age_first_blood_pressure),                  
                      age_first_clopi_dipy = scale(age_first_clopi_dipy),                         
                      age_first_breast_canc = scale(age_first_breast_cancer),                           
                      age_first_doac = scale(age_first_doac),
                      PC1 = scale(PC1),
                      PC2 = scale(PC2),
                      PC3 = scale(PC3),
                      PC4 = scale(PC4),
                      PC5 = scale(PC5),
                      PC6 = scale(PC6),
                      PC7 = scale(PC7),
                      PC8 = scale(PC8),
                      PC9 = scale(PC9),
                      PC10 = scale(PC10),
                      statins = factor(statins),
                      fluvastatins = factor(fluvastatins),
                      rosuvastatins = factor(rosuvastatins),
                      blood_pressure = factor(blood_pressure),
                      breast_canc = factor(breast_canc),
                      doac = factor(doac),
                      clopi_dipy = factor(clopi_dipy))

# # # Drug-genes specific analysis
drugs <- c("statins", "fluvastatins", "rosuvastatins", "breast_canc", "clopi_dipy")

for (d in drugs){
  
  if (d == "statins"){
    g <- "slco1b1"
  } else if (d == "fluvastatins") {
    g <- "cyp2c9"
  } else if (d == "rosuvastatins") {
    g <- "abcg2"
  } else if (d == "breast_canc") {
    g <- "cyp2d6"
  } else if (d == "clopi_dipy") {
    g <- "cyp2c19"
  }
  
  if (grepl("statin", d)) {
    drug <- "statins"
  } else {
    drug <- d
  }
  
  dip <- fread(paste0('/home/ivm/drugs/data/pharmacogenes/stargazer_call/',g,'_report.tsv'))
  
  # # # ADHERENCE
  print(paste0("Adherence to ", d))
  
  # Merge adherence and star allele df and filter for only clopidogrel, tamoxifen or specific type of statins
  adh_dip <- left_join(adh, dip, by = c("FINNGENID" = "Sample"))
  
  subdata <- adh_dip %>%
    select(paste0(drug, "_adh"), paste0("age_first_",drug), SEX_IMPUTED, PC1:PC10, Score, Phenotype, ends_with("prop")) %>% 
    filter(Score >= 0) %>%
    mutate(SEX_IMPUTED = factor(SEX_IMPUTED))
  
  if (g == "slco1b1") {
    subdata <- subdata %>%
      select(-tamoxifen_prop, -clopidogrel_prop, -rosuvastatin_prop, -fluvastatin_prop) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_function", "poor_function", "decreased_function", "increased_function")))
  } else if (g == "cyp2c9") {
    subdata <- subdata %>%
      select(-tamoxifen_prop, -clopidogrel_prop, -rosuvastatin_prop) %>% 
      filter(fluvastatin_prop == 1) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_metabolizer", "intermediate_metabolizer", "poor_metabolizer")))
  } else if (g == "abcg2") {
    subdata <- subdata %>%
      select(-tamoxifen_prop, -clopidogrel_prop, -fluvastatin_prop) %>% 
      filter(rosuvastatin_prop == 1) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = case_when(Score == 1 ~ "poor_function",
                                   Score == 1.5 ~ "decreased_function",
                                   Score == 2 ~ "normal_function"),
             Phenotype = factor(Phenotype, levels = c("normal_function", "poor_function", "decreased_function")))
  } else if (g == "cyp2d6") {
    subdata <- subdata %>%
      select(-clopidogrel_prop, -rosuvastatin_prop, -fluvastatin_prop) %>% 
      filter(tamoxifen_prop == 1) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_metabolizer", "intermediate_metabolizer", "poor_metabolizer")))
  } else if (g == "cyp2c19"){
    subdata <- subdata %>%
      select(-tamoxifen_prop, -rosuvastatin_prop, -fluvastatin_prop) %>% 
      filter(clopidogrel_prop == 1) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_metabolizer", "intermediate_metabolizer", "poor_metabolizer", "ultrarapid_metabolizer")))
  }
  
  if (d == "breast_canc"){
    model <- paste0(drug,"_adh ~  + age_first_",drug," +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Phenotype")
  } else {
    model <- paste0(drug,"_adh ~ SEX_IMPUTED + age_first_",drug," +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Phenotype")
  }
  
  m <- lm(formula = model, data = subdata)
  
  pheno_counts <- as.data.frame(table(subdata$Phenotype))
  
  tt <- m %>%
    tidy() %>% 
    filter(startsWith(term, "Phenotype")) %>% 
    select(term, estimate, std.error, p.value) %>% 
    mutate(term = gsub("Phenotype", "", term)) %>% 
    left_join(pheno_counts, by = c("term" = "Var1")) %>% 
    mutate(gene = g,
           drug = d,
           N_tot = nrow(subdata)) %>% 
    select(gene, drug, N_tot, term, N = Freq, everything())
  
  fwrite(tt, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence_",d,"_",g,".tsv"), sep = "\t")
  
  # # # PERSISTENCE
  print(paste0("Persistence to ", d))
  
  # Merge adherence and star allele df and filter for only clopidogrel, tamoxifen or specific type of statins
  per_dip <- left_join(per, dip, by = c("FINNGENID" = "Sample"))
  
  subdata <- per_dip %>%
    select(paste0(d), paste0("age_first_",d), SEX_IMPUTED, PC1:PC10, Score, Phenotype, ends_with("prop")) %>% 
    filter(Score >= 0) %>%
    mutate(SEX_IMPUTED = factor(SEX_IMPUTED))
  
  if (g == "slco1b1") {
    subdata <- subdata %>%
      select(-tamoxifen_prop, -clopidogrel_prop) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_function", "poor_function", "decreased_function", "increased_function")))
  } else if (g == "cyp2c9") {
    subdata <- subdata %>%
      select(-tamoxifen_prop, -clopidogrel_prop) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_metabolizer", "intermediate_metabolizer", "poor_metabolizer")))
  } else if (g == "abcg2") {
    subdata <- subdata %>%
      select(-tamoxifen_prop, -clopidogrel_prop) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = case_when(Score == 1 ~ "poor_function",
                                   Score == 1.5 ~ "decreased_function",
                                   Score == 2 ~ "normal_function"),
             Phenotype = factor(Phenotype, levels = c("normal_function", "poor_function", "decreased_function")))
  } else if (g == "cyp2d6") {
    subdata <- subdata %>%
      select(-clopidogrel_prop) %>% 
      filter(tamoxifen_prop == 1) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_metabolizer", "intermediate_metabolizer", "poor_metabolizer")))
  } else if (g == "cyp2c19"){
    subdata <- subdata %>%
      select(-tamoxifen_prop) %>% 
      filter(clopidogrel_prop == 1) %>% 
      filter(complete.cases(.)) %>% 
      mutate(Phenotype = factor(Phenotype, levels = c("normal_metabolizer", "intermediate_metabolizer", "poor_metabolizer", "ultrarapid_metabolizer")))
  }
  
  if (d == "breast_canc"){
    model <- paste0(d," ~  + age_first_", d," +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Phenotype")
  } else {
    model <- paste0(d," ~ SEX_IMPUTED + age_first_",d," +
      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Phenotype")
  }
  
  m <- glm(formula = model, data = subdata, family = "binomial")
  
  pheno_counts <- as.data.frame(table(subdata$Phenotype, subdata[[d]]))
  
  tt <- m %>%
    tidy() %>% 
    filter(startsWith(term, "Phenotype")) %>% 
    select(term, estimate, std.error, p.value) %>% 
    mutate(term = gsub("Phenotype", "", term)) %>% 
    left_join(pheno_counts %>% filter(Var2 == 0) %>% select(Var1, N_dis = Freq), by = c("term" = "Var1")) %>% 
    left_join(pheno_counts %>% filter(Var2 == 1) %>% select(Var1, N_per = Freq), by = c("term" = "Var1")) %>% 
    mutate(gene = g,
           drug = d,
           N_tot = nrow(subdata)) %>% 
    select(gene, drug, N_tot, term, N_dis, N_per, everything())
  
  fwrite(tt, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence_",d,"_",g,".tsv"), sep = "\t")
  
}

RES_ADH <- NULL
RES_PER <- NULL

for (d in drugs) {
  if (d == "statins"){
    g <- "slco1b1"
  } else if (d == "fluvastatins") {
    g <- "cyp2c9"
  } else if (d == "rosuvastatins") {
    g <- "abcg2"
  } else if (d == "breast_canc") {
    g <- "cyp2d6"
  } else if (d == "clopi_dipy") {
    g <- "cyp2c19"
  }
  
  a <- fread(paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence_",d,"_",g,".tsv"))
  RES_ADH <- bind_rows(RES_ADH, a)
  RES_ADH[which(RES_ADH$N < 5), c("estimate", "std.error", "p.value")] <- NA
  
  p <- fread(paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence_",d,"_",g,".tsv"))
  RES_PER <- bind_rows(RES_PER, p)
  RES_PER[which(RES_PER$N_dis < 5 | RES_PER$N_per < 5), c("estimate", "std.error", "p.value")] <- NA
}

fwrite(RES_ADH, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence_all.tsv"), sep = "\t", na = "NA")
fwrite(RES_PER, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence_all.tsv"), sep = "\t", na = "NA")
