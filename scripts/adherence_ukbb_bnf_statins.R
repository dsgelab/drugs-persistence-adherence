rm(list=ls())

library(ggplot2)
library(data.table)
library(tidyverse)
library(lubridate)
library(readxl)

# # # Pre-process
# read_excel('drugs/data/ukbb/BNF Snomed Mapping data 20210218.xlsx') %>%
#   as.data.frame() %>%
#   select(bnf_code = `BNF Code`, bnf_name = `BNF Name`, SNOMED_code = `SNOMED Code`, strength = Strength,
#          unit = `Unit Of Measure`, pkg_size = `Pack`) %>%
#   fwrite('drugs/data/ukbb/bnf_snomed_mapping_processed.tsv', sep = '\t')

# read in files
scripts <- fread('~/drugs/data/ukbb/ukb31063.gp_scripts.20200706w.processed.bnf_ok.date_ok.txt')

# Complete trajectories
# st <- getTrajectoriesUKB(scripts,'^0212000(AC|Y0|X0|M0|B0|C0|AA)') --> <200k purchases when grep with bnf 


# extract purchases with BNF codes
st_all_bnf <- scripts %>% 
  filter(grepl('^0212000(AC|Y0|X0|M0|B0|C0|AA)', bnf_code))

# extract purchases with statins names
st_all_str <- scripts %>% 
  filter(grepl('SIMVASTATIN', toupper(drug_name)) |
           grepl('LOVASTATIN', toupper(drug_name)) |
           grepl('PRAVASTATIN', toupper(drug_name)) |
           grepl('FLUVASTATIN', toupper(drug_name)) |
           grepl('ATORVASTATIN', toupper(drug_name)) |
           grepl('CERIVASTATIN', toupper(drug_name)) |
           grepl('ROSUVASTATIN', toupper(drug_name)))

# check if drug_name refers to the substanceor the product name --> it is the product commercial name (Igeny = simsvastatin)
st_inegy_str <- scripts %>% 
  filter(grepl('INEGY', toupper(drug_name)))

# N purchases extracted using BNF
st_all_bnf %>% 
  distinct() %>% 
  nrow()
# 198,407

# N purchases extracted using name
st_all_str %>% 
  distinct() %>% 
  nrow()
# 3,107,032

# N purchases we can add using both methods
st_all_bnf %>% 
  bind_rows(st_all_str) %>% 
  distinct() %>% 
  nrow()
# ~ +400

# compare BNF and drug names not picked by BNF
st_all_bnf_codes <- st_all_bnf %>% 
  select(bnf_code, drug_name) %>% 
  distinct()

st_all_str_codes <- st_all_str %>% 
  select(bnf_code, drug_name) %>% 
  filter(!bnf_code %in% st_all_bnf$bnf_code) %>% 
  group_by(bnf_code, drug_name) %>% 
  summarise(n=n())

unique(st_all_str_codes$bnf_code)
table(st_all_str$bnf_code)

statins_021204 <- st_all_str %>% 
  filter(startsWith(bnf_code, "02120400")) %>% 
  mutate(name = tolower(sub("([A-Za-z]+).*", "\\1", drug_name))) %>% 
  group_by(name) %>% 
  summarise(n=n()) %>% 
  arrange(-n)

# SOLUTION to pick the most possible number of records:
# - get all possible commercial names from the lookup table, grep-ing by substance name
# - get purchases using product names

bnf_lkp <- fread('~/drugs/data/ukbb/all_lkps_maps_v2__bnf_lkp.csv')

product_names <- bnf_lkp %>% 
  filter(grepl('SIMVASTATIN', toupper(BNF_Chemical_Substance)) |
           grepl('LOVASTATIN', toupper(BNF_Chemical_Substance)) |
           grepl('PRAVASTATIN', toupper(BNF_Chemical_Substance)) |
           grepl('FLUVASTATIN', toupper(BNF_Chemical_Substance)) |
           grepl('ATORVASTATIN', toupper(BNF_Chemical_Substance)) |
           grepl('CERIVASTATIN', toupper(BNF_Chemical_Substance)) |
           grepl('ROSUVASTATIN', toupper(BNF_Chemical_Substance))) %>% 
  select(BNF_Chemical_Substance, BNF_Product) %>% 
  distinct() %>% 
  arrange(BNF_Chemical_Substance, BNF_Product) %>% 
  mutate(name = tolower(sub("([A-Za-z]+).*", "\\1", BNF_Product)))


# Combinations?
names <- paste(unique(product_names$name),collapse = "|")

st_all <- scripts %>% 
  filter(# grepl based on all possible commercial names
    grepl(names, tolower(drug_name)),
    # ensure we only pick up statins, and not other drugs containing part of the names 
    grepl("^0212",bnf_code))

st_all_names <- st_all %>% 
  select(bnf_code, drug_name) %>% 
  group_by(bnf_code, drug_name) %>% 
  summarise(n=n()) %>% 
  arrange(drug_name)

# Remove fenofibrate
st_all_names <- st_all_names %>% 
  filter(!grepl('fenofibrate', tolower(drug_name)))

fwrite(st_all_names, '~/drugs/data/ukbb/ukbb_statins.tsv', sep = "\t")

# extract all statins scripts and reformat quantities
st_list <- fread('~/drugs/data/ukbb/ukbb_statins.tsv') %>% 
  pull(drug_name)

# pull quantities
q <- scripts %>% 
  filter(drug_name %in% st_list) %>% 
  select(quantity) %>% 
  group_by(quantity) %>% 
  summarise(n=n())

# functions to transform quantities
tablet <- function(x){
  return(as.integer(gsub(".*?([0-9]+).*", "\\1", x)))
}
pack_per_tablet <- Vectorize(function(x) {
  qts <- as.integer(regmatches(x, gregexpr('[0-9]+',x))[[1]])
  return(ifelse(length(qts) == 1, qts[1], qts[1]*qts[2]))
}, vectorize.args = "x")

# extract all statins prescriptions and transform quantities to integer
st_all <- scripts %>%
  filter(drug_name %in% st_list) %>%
  mutate(quantity_num = case_when(
    # e.g. "28 tablet", "56 tabs"
    grepl('tablet|tab|capsule', tolower(quantity)) &
      !grepl('pack', tolower(quantity)) &
      !grepl('\\*', tolower(quantity)) &
      !grepl('x', tolower(quantity)) ~ tablet(quantity),
    # e.g. "2 packs of 28 tablet(s)", "2*28 tablet - 40 mg", "2x28tablet(s)"
    grepl('tablet|tab', tolower(quantity)) &
      (grepl('pack', tolower(quantity)) |
         grepl('\\*', tolower(quantity)) |
         grepl('x', tolower(quantity)) ) ~ pack_per_tablet(quantity),
    TRUE ~ as.integer(quantity)
  ))

fwrite(st_all, '~/drugs/data/ukbb/UKB_all_statins.tsv', quote = F, sep = "\t")

st_scripts <- fread('~/drugs/data/ukbb/UKB_all_statins.tsv')

nrow(st_scripts)
sum(is.na(st_scripts$quantity_num))
sum(!is.na(st_scripts$quantity_num))

check_data_provider <- purch %>% 
  group_by(eid) %>% 
  summarize(data_p = sum(data_provider)/n()) %>% 
  filter(data_p != 2,
         data_p != 3)
purch$eid <- as.character(purch$eid)

