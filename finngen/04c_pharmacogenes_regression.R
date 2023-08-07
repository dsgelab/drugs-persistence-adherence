# # # Association between each diplotype and persistence and adherence

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

dd <- fread("/home/ivm/drugs/data/R10_blood_pressure_summarized.txt") %>% mutate(blood_pressure_adh = adherence) %>%
  select(FINNGENID,blood_pressure_adh)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_breast_cancer_summarized.txt") %>% mutate(breast_canc_adh = adherence) %>%
  select(FINNGENID,breast_canc_adh,tamoxifen_prop=prop_atc,tamoxifen_atc=atc)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_clopi_dipy_summarized.txt") %>% mutate(clopi_dipy_adh = adherence) %>%
  select(FINNGENID,clopi_dipy_adh, clopidogrel_prop=prop_atc,clopidogrel_atc=atc)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_doac_summarized.txt") %>% mutate(doac_adh = adherence) %>%
  select(FINNGENID,doac_adh)
adh <- merge(x = adh, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_statins_summarized.txt") %>% mutate(statins_adh = adherence) %>%
  select(FINNGENID,statins_adh)
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
per <- per %>% mutate(age_first_statins = scale(age_first_statins),                                
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
                      blood_pressure = factor(blood_pressure),
                      breast_canc = factor(breast_canc),
                      doac = factor(doac),
                      clopi_dipy = factor(clopi_dipy))

dd <- fread("/home/ivm/drugs/data/R10_breast_cancer_persistence.txt") %>%
  select(FINNGENID,tamoxifen_prop=prop_atc,tamoxifen_atc=atc)
per <- merge(x = per, y = dd , by = "FINNGENID", all.x=TRUE)

dd <- fread("/home/ivm/drugs/data/R10_clopi_dipy_persistence.txt") %>% 
  select(FINNGENID,clopidogrel_prop=prop_atc,clopidogrel_atc=atc)
per <- merge(x = per, y = dd , by = "FINNGENID", all.x=TRUE)


# # # Diplotype frequencies
frq <- fread('/home/ivm/drugs/results/manuscript/pharmacogenes/diplotype_freq.tsv')

# # # Drug-genes specific analysis

drugs <- c("statins", "breast_canc", "clopi_dipy")

for (d in drugs){
  
  if (d == "statins"){
    genes <- c("slco1b1", "cyp2c9", "abcg2")  
  } else if (d == "breast_canc"){
    genes <- c("cyp2d6")
  } else if (d == "clopi_dipy"){
    genes <- c("cyp2c19")
  }
  
  for (g in genes){
    dip <- fread(paste0('/home/ivm/drugs/data/pharmacogenes/stargazer_call/',g,'_report.tsv'))
    frq_g <- frq %>% filter(gene == g)
    # Reference (*1/*1 or normal function VS all others)
    dip_reference <- dip %>% filter(Score == 2) %>% pull(Diplotype) %>% unique()
    # diplotypes to test (score >= 0 to exclude those marked as "unknown function" or "unknonw metabolizer")
    # test only diplotypes w frq > 1%
    dip_to_test <- dip %>%
      filter(Score != 2, Score >= 0) %>%
      #left_join(frq_g, by = c("Diplotype" = "diplotype")) %>% 
      #filter(freq > 0.01) %>% 
      pull(Diplotype) %>%
      unique()
    
    
    # # # ADHERENCE
    print(paste0("Adherence to ", d))
    # merge adherence and star allele df and filter for only clopidogrel or tamoxifen
    adh_dip <- left_join(adh, dip, by = c("FINNGENID" = "Sample"))

    if (d == "breast_canc"){
      adh_dip <- adh_dip %>% filter(tamoxifen_prop == 1)
    } else if (d == "clopi_dipy"){
      adh_dip <- adh_dip %>% filter(clopidogrel_prop == 1)
    }
    
    betas <- NULL
    for (di in dip_to_test){
      adh_dip <- adh_dip %>% 
        mutate(diplo = case_when(Diplotype %in% dip_reference ~ 0,
                                 Diplotype == di ~ 1,
                                 TRUE ~ NA_integer_))
      pheno <- adh_dip %>% filter(diplo == 1) %>% pull(Phenotype) %>% unique()
      
      subdata <- adh_dip %>%
        mutate(SEX_IMPUTED = factor(SEX_IMPUTED),
               diplo = factor(diplo)) %>% 
        select(paste0(d,"_adh"), paste0("age_first_",d), SEX_IMPUTED, PC1:PC10, diplo) %>% 
        filter(complete.cases(.))
      
      dip_groups <- length(table(subdata$diplo))
      min_dip_count <- min(table(subdata$diplo))
      
      if (dip_groups == 2 & min_dip_count > 5){
        print(di)
        if (d == "breast_canc") {
          model <- paste0(d,"_adh ~  + age_first_",d," +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
        } else {
          model <- paste0(d,"_adh ~ SEX_IMPUTED + age_first_",d," +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
        }
        
        m <- lm(formula = model, subdata)
        
        tt <- m %>%
          tidy() %>%
          filter(term == "diplo1") %>%
          mutate(drug = d,
                 gene = g,
                 star_allele = di,
                 pheno = pheno,
                 N = nrow(subdata),
                 N_ref = table(subdata$diplo)["0"],
                 N_dip = table(subdata$diplo)["1"]) %>% 
          select(drug, gene, star_allele, pheno, N, N_ref, N_dip, estimate, std.error, p.value)
      } else {
        tt <- data.frame(drug = d,
                         gene = g,
                         star_allele = di,
                         pheno = ifelse(length(pheno)==0, NA, pheno),
                         N = nrow(subdata),
                         N_ref = table(subdata$diplo)["0"],
                         N_dip = table(subdata$diplo)["1"],
                         estimate = NA,
                         std.error = NA,
                         p.value = NA)
      }
      betas <- rbind(betas,tt)
    }
    betas <- betas %>% mutate(fdr = p.adjust(p.value, "fdr"))
    fwrite(betas, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence_",d,"_",g,".tsv"), sep = "\t")
    
    # # # PERSISTENCE
    print(paste0("Persistence to ", d))
    per_dip <- left_join(per, dip, by = c("FINNGENID" = "Sample"))
   
    if (d == "breast_canc"){
      per_dip <- per_dip %>% filter((breast_canc == 1 & tamoxifen_prop == 1) | (breast_canc == 0 & tamoxifen_atc == "L02BA01"))
    } else if (d == "clopi_dipy"){
      per_dip <- per_dip %>% filter((clopi_dipy == 1 & clopidogrel_prop == 1) | (clopi_dipy == 0 & clopidogrel_atc == "B01AC04"))
    }
    
    betas <- NULL
    for (di in dip_to_test){
      per_dip <- per_dip %>% 
        mutate(diplo = case_when(Diplotype %in% dip_reference ~ 0,
                                 Diplotype == di ~ 1,
                                 TRUE ~ NA_integer_))
      pheno <- per_dip %>% filter(diplo == 1) %>% pull(Phenotype) %>% unique()
      
      subdata <- per_dip %>%
        mutate(SEX_IMPUTED = factor(SEX_IMPUTED),
               diplo = factor(diplo)) %>% 
        select(paste0(d), paste0("age_first_",d), SEX_IMPUTED, PC1:PC10, diplo) %>% 
        filter(complete.cases(.))
      
      dip_groups <- length(table(subdata$diplo))
      min_dip_count <- min(table(subdata[[d]], subdata$diplo))
      
      if (dip_groups == 2 & min_dip_count > 5){
        print(di)
        if (d == "breast_canc") {
          model <- paste0(d," ~  + age_first_",d," +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
        } else {
          model <- paste0(d," ~ SEX_IMPUTED + age_first_",d," +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
        }
        
        m <- glm(formula = model, data = subdata, family = "binomial")
        
        tt <- m %>%
          tidy() %>%
          filter(term == "diplo1") %>%
          mutate(drug = d,
                 gene = g,
                 star_allele = di,
                 pheno = pheno,
                 N = nrow(subdata),
                 N_pers = table(subdata[[d]])["1"],
                 N_ref = table(subdata$diplo)["0"],
                 N_dip = table(subdata$diplo)["1"],
                 N_dip_pers = table(subdata[[d]], subdata$diplo)["1","1"],
                 N_dip_non_pers = table(subdata[[d]], subdata$diplo)["0","1"]
                 ) %>% 
          select(drug, gene, star_allele, pheno, N, N_pers, N_ref, N_dip, N_dip_pers, N_dip_non_pers, estimate, std.error, p.value)
      } else {
        tt <- data.frame(drug = d,
                         gene = g,
                         star_allele = di,
                         pheno = ifelse(length(pheno)!=0, pheno,NA),
                         N = nrow(subdata),
                         N_pers = table(subdata[[d]])["1"],
                         N_ref = table(subdata$diplo)["0"],
                         N_dip = table(subdata$diplo)["1"],
                         N_dip_pers = NA,
                         N_dip_non_pers = NA,
                         estimate = NA,
                         std.error = NA,
                         p.value = NA)
      }
      betas <- rbind(betas,tt)
    }
    betas <- betas %>% mutate(fdr = p.adjust(p.value, "fdr"))
    fwrite(betas, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence_",d,"_",g,".tsv"), sep = "\t")
  }
}

# # # # CYP2D6 across medications analysis
# drugs <- c("statins", "breast_canc", "clopi_dipy", "blood_pressure", "doac")
# 
# for (d in drugs){
#   
#   g <- "cyp2d6"
#   
#   dip <- fread(paste0('/home/ivm/drugs/data/pharmacogenes/stargazer_call/',g,'_report.tsv'))
#   # Reference (*1/*1 or normal function VS all others)
#   dip_reference <- dip %>% filter(Score == 2) %>% pull(Diplotype) %>% unique()
#   # diplotypes to test (score >= 0 to exclude those marked as "unknown function" or "unknonw metabolizer")
#   dip_to_test <- dip %>% filter(Score != 2, Score >= 0) %>% pull(Diplotype) %>% unique()
#   
#   # Adherence
#   adh_dip <- left_join(adh, dip, by = c("FINNGENID" = "Sample"))
#   
#   betas <- NULL
#   for (di in dip_to_test){
#     adh_dip <- adh_dip %>% 
#       mutate(diplo = case_when(Diplotype %in% dip_reference ~ 0,
#                                Diplotype == di ~ 1,
#                                TRUE ~ NA_integer_))
#     pheno <- adh_dip %>% filter(diplo == 1) %>% pull(Phenotype) %>% unique()
#     
#     subdata <- adh_dip %>%
#       mutate(SEX_IMPUTED = factor(SEX_IMPUTED),
#              diplo = factor(diplo)) %>% 
#       select(paste0(d,"_adh"), paste0("age_first_",d), SEX_IMPUTED, PC1:PC10, diplo) %>% 
#       filter(complete.cases(.))
#     
#     dip_count <- table(subdata$diplo)[[2]]
#     
#     if (dip_count != 0){
#       if (d == "breast_canc") {
#         model <- paste0(d,"_adh ~  + age_first_",d," +
#       PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
#       } else {
#         model <- paste0(d,"_adh ~ SEX_IMPUTED + age_first_",d," +
#       PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
#       }
#       
#       m <- lm(formula = model, subdata)
#       
#       tt <- m %>%
#         tidy() %>%
#         filter(term == "diplo1") %>%
#         mutate(drug = d,
#                gene = g,
#                star_allele = di,
#                pheno = pheno,
#                N = nrow(subdata)) %>% 
#         select(drug, gene, star_allele, pheno, N, estimate, std.error, p.value)
#     } else {
#       tt <- data.frame(drug = d,
#                        gene = g,
#                        star_allele = di,
#                        pheno = pheno,
#                        N = 0,
#                        estimate = NA,
#                        std.error = NA,
#                        p.value = NA)
#     }
#     betas <- rbind(betas,tt)
#   }
#   betas <- betas %>% mutate(fdr = p.adjust(p.value, "fdr"))
#   fwrite(betas, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_adherence_",d,"_",g,".tsv"), sep = "\t")
#   
#   # Persistence
#   per_dip <- left_join(per, dip, by = c("FINNGENID" = "Sample"))
#   
#   betas <- NULL
#   
#   for (di in dip_to_test){
#     per_dip <- per_dip %>% 
#       mutate(diplo = case_when(Diplotype %in% dip_reference ~ 0,
#                                Diplotype == di ~ 1,
#                                TRUE ~ NA_integer_))
#     pheno <- per_dip %>% filter(diplo == 1) %>% pull(Phenotype) %>% unique()
#     
#     subdata <- per_dip %>%
#       mutate(SEX_IMPUTED = factor(SEX_IMPUTED),
#              diplo = factor(diplo)) %>% 
#       select(paste0(d), paste0("age_first_",d), SEX_IMPUTED, PC1:PC10, diplo) %>% 
#       filter(complete.cases(.))
#     
#     min_cont_table <- min(table(subdata[[d]], subdata$diplo))
#     
#     if (min_cont_table > 5){
#       if (d == "breast_canc") {
#         model <- paste0(d," ~  + age_first_",d," +
#       PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
#       } else {
#         model <- paste0(d," ~ SEX_IMPUTED + age_first_",d," +
#       PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + diplo")
#       }
#       
#       m <- glm(formula = model, data = subdata, family = "binomial")
#       
#       tt <- m %>%
#         tidy() %>%
#         filter(term == "diplo1") %>%
#         mutate(drug = d,
#                gene = g,
#                star_allele = di,
#                pheno = pheno,
#                N = nrow(subdata)) %>% 
#         select(drug, gene, star_allele, pheno, N, estimate, std.error, p.value)
#     } else {
#       tt <- data.frame(drug = d,
#                        gene = g,
#                        star_allele = di,
#                        pheno = pheno,
#                        N = 0,
#                        estimate = NA,
#                        std.error = NA,
#                        p.value = NA)
#     }
#     betas <- rbind(betas,tt)
#   }
#   betas <- betas %>% mutate(fdr = p.adjust(p.value, "fdr"))
#   fwrite(betas, paste0("/home/ivm/drugs/results/manuscript/pharmacogenes/betas/betas_persistence_",d,"_",g,".tsv"), sep = "\t")
# }