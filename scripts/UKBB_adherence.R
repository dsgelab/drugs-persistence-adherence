#### SET ####
rm(list=ls())

library(ggplot2)
library(data.table)
library(lubridate)
library(readxl)
library(dplyr)
library(Hmisc)

# functions to transform quantities
tablet <- function(x){
  return(as.integer(gsub(".*?([0-9]+).*", "\\1", x)))
}
pack_per_tablet <- Vectorize(function(x) {
  qts <- as.integer(regmatches(x, gregexpr('[0-9]+',x))[[1]])
  return(ifelse(length(qts) == 1, qts[1], qts[1]*qts[2]))
}, vectorize.args = "x")


#### READ ####
# read in files
scripts <- fread('~/data/gp_scripts.20200706w.processed.bnf_ok.date_ok.txt')
#read map atc-nbf
map <- fread("~/data/bnf_atc_map.tsv",sep= "\t")
#read lookup table
bnf_lkp <- fread('~/data/all_lkps_maps_v2__bnf_lkp.csv')
#read atc to substance table
atc_map <- fread("~/data/atc_codes_wikipedia.csv")

colnames(map) <- c("bnf","code")
mappa <- merge(map,atc_map,by="code",all=T)


#### STATINS ####

#get atc and substance name
st_atc <- c("C10AA01",
"C10AA01",
"C10AA02",
"C10AA03",
"C10AA04",
"C10AA05",
"C10AA06",
"C10AA07") 

st_sub <- atc_map[atc_map$code %in% st_atc,]$desc %>% tolower() %>% as.character() %>% sort()
st_sub_str <- paste(st_sub,collapse = "|")

#get commercial name from lookup table
st_lkp <- bnf_lkp %>% filter(grepl(st_sub_str,tolower(BNF_Chemical_Substance)))
st_com <- st_lkp$BNF_Product %>% unique %>% tolower() %>% as.character() %>% sort()


#merge commercial and substance names
st_names <- unique(c(st_com,st_sub)) %>% sort()
st_names <- st_names[!grepl("/",st_names)]
st_names <- gsub( " .*$", "", st_names)

st_names_str <- paste(st_names,collapse = "|")

start_time <- Sys.time()
st<- scripts %>% filter(grepl("^0212",bnf_code))  %>% filter(grepl(st_names_str,tolower(drug_name)))
end_time <- Sys.time()
end_time- start_time


length(st$eid)
#3236105

# N s.u. extracted using BNF
st$eid %>% 
  unique() %>%
  length()
# 63441 s.u.




# convert quantity to numeric
st <- st %>%
  mutate(quantity_num = case_when(
    # e.g. "28 tablet", "56 tabs"
    grepl('tablet|tab|cap', tolower(quantity)) &
      !grepl('pack', tolower(quantity)) &
      !grepl('\\*', tolower(quantity)) &
      !grepl('x', tolower(quantity)) &
      !grepl('*', tolower(quantity)) ~ tablet(quantity),
    grepl('tablet|tab|cap', tolower(quantity)) &
      (grepl('pack', tolower(quantity)) |
         grepl('\\*', tolower(quantity)) |
         grepl('x', tolower(quantity)) |
         grepl('*', tolower(quantity)) ) ~ pack_per_tablet(quantity),
    TRUE ~ as.integer(quantity)
  ))



table(is.na(st$quantity_num))
tail(st[is.na(st$quantity_num),])

#write table
fwrite(st,"~/data/st_scripts.tsv",quote = F, row.names = F,sep="\t")


#### BLOOD PRESSURE ####

#get atc and substance name
bp_atc <- c("C02AB01",
            "C02AC01",
            "C02AC02",
            "C02AC05",
            "C02CA01",
            "C02DC01",
            "C02KX01",
            "C02KX02",
            "C02KX04",
            "C02KX05",
            "C02LA01",
            "C03AA03",
            "C03BA11",
            "C03CA01",
            "C03CA02",
            "C03DA01",
            "C03DB01",
            "C03DB02",
            "C03EA01",
            "C03EA02",
            "C03EB01",
            "C08CA01",
            "C08CA02",
            "C08CA03",
            "C08CA05",
            "C08CA06",
            "C08CA07",
            "C08CA10",
            "C08CA13",
            "C08CX01",
            "C08DA01",
            "C08DB01",
            "C09AA01",
            "C09AA02",
            "C09AA03",
            "C09AA04",
            "C09AA05",
            "C09AA06",
            "C09AA08",
            "C09AA10",
            "C09AA16",
            "C09BA02",
            "C09BA03",
            "C09BA04",
            "C09BA05",
            "C09BA06",
            "C09BB02",
            "C09BB04",
            "C09BB05",
            "C09BB10",
            "C09BX01",
            "C09BX02",
            "C09CA",
            "C09CA01",
            "C09CA02",
            "C09CA03",
            "C09CA04",
            "C09CA06",
            "C09CA07",
            "C09CA08",
            "C09DA01",
            "C09DA02",
            "C09DA03",
            "C09DA06",
            "C09DA07",
            "C09DA08",
            "C09DB01",
            "C09DB02",
            "C09DX01",
            "C09DX04",
            "C09XA02") 
bp_sub <- atc_map[atc_map$code %in% bp_atc,]$desc %>% tolower() %>% as.character() %>% sort()
bp_sub <- bp_sub[!grepl("diuretics| and ",bp_sub)]
bp_sub <- bp_sub[!bp_sub %in%c("methyldopa (levorotatory)",
                               "angiotensin ii receptor blockers (arbs), plain",
                               "olmesartan medoxomil")]
bp_sub <- append(bp_sub,c("sacubitril","methyldopa", "olmesartan"))
bp_sub_str <- paste(bp_sub,collapse = "|")

#get commercial name from lookup table
bp_lkp <- bnf_lkp %>% filter(grepl(bp_sub_str,tolower(BNF_Chemical_Substance)))
bp_com <- bp_lkp$BNF_Product %>% unique %>% tolower() %>% as.character() %>% sort()

comm <- strsplit(bp_com,"[ ]")
com_name <- c()
for (i in 1:364) {
  com_name <- append(com_name,comm[[i]][1])
}
bp_com <-  unique(com_name)

#merge commercial and substance names
bp_names <- unique(c(bp_com,bp_sub)) %>% sort()
bp_names <- bp_names[!grepl("/",bp_names)]
bp_names <- bp_names[! bp_names %in% c("burinex-a",
                                       "burinex-k",
                                       "gx",
                                       "hydrosaluric-k")]

bp_names <- append(bp_names,"sacubitril")
bp_names_str <- paste(bp_names,collapse = "|")

start_time <- Sys.time()
bp <- scripts %>% filter(grepl("^02",bnf_code)) %>% filter(grepl(bp_names_str,tolower(drug_name)))
end_time <- Sys.time()
end_time- start_time

length(bp$eid)
#5011220

# N s.u. extracted using BNF
bp$eid %>% 
  unique() %>%
  length()
# 67430 s.u.




# convert quantity to numeric
bp <- bp %>%
  mutate(quantity_num = case_when(
    # e.g. "28 tablet", "56 tabs"
    grepl('tablet|tab|cap', tolower(quantity)) &
      !grepl('pack', tolower(quantity)) &
      !grepl('\\*', tolower(quantity)) &
      !grepl('x', tolower(quantity)) &
      !grepl('*', tolower(quantity)) ~ tablet(quantity),
    grepl('tablet|tab|cap', tolower(quantity)) &
      (grepl('pack', tolower(quantity)) |
         grepl('\\*', tolower(quantity)) |
         grepl('x', tolower(quantity)) |
         grepl('*', tolower(quantity)) ) ~ pack_per_tablet(quantity),
    TRUE ~ as.integer(quantity)
  ))


table(is.na(bp$quantity_num))
head(bp[is.na(bp$quantity_num),])
tail(bp[is.na(bp$quantity_num),])


#write table
fwrite(bp,"~/data/bp_scripts.tsv",quote = F, row.names = F,sep="\t")


#### ANTIPLATELET ####

#get atc and substance name
ap_atc <- c("B01AC04", "B01AC07") 
ap_sub <- atc_map[atc_map$code %in% ap_atc,]$desc %>% tolower() %>% as.character()
ap_sub_str <- paste(ap_sub,collapse = "|")

#get commercial name from lookup table
ap_lkp <- bnf_lkp %>% filter(grepl(ap_sub_str,tolower(BNF_Chemical_Substance)))
ap_com <- ap_lkp$BNF_Product %>% unique %>% tolower() %>% as.character() 

#merge commercial and substance names
ap_names <- unique(c(ap_com,ap_sub))
ap_names <- ap_names[! ap_names %in% c("dipyridamole/aspirin","ofcram pr")]
ap_names <- append(ap_names,"ofcram")
ap_names_str <- paste(ap_names,collapse = "|")

#get doac prescriptions
ap<- scripts %>%  filter(grepl("^0209",bnf_code)) %>%  filter(grepl(ap_names_str,tolower(drug_name)))
length(ap$eid)
#247157

# N s.u. extracted using BNF
ap$eid %>% 
  unique() %>%
  length()
# 10183 s.u.




#convert quantity to numeric
ap <- ap %>%
  mutate(quantity_num = case_when(
    # e.g. "28 tablet", "56 tabs"
    grepl('tablet|tab|cap', tolower(quantity)) &
      !grepl('pack', tolower(quantity)) &
      !grepl('\\*', tolower(quantity)) &
      !grepl('x', tolower(quantity)) &
      !grepl('*', tolower(quantity)) ~ tablet(quantity),
    grepl('tablet|tab|cap', tolower(quantity)) &
      (grepl('pack', tolower(quantity)) |
         grepl('\\*', tolower(quantity)) |
         grepl('x', tolower(quantity)) |
         grepl('*', tolower(quantity)) ) ~ pack_per_tablet(quantity),
    TRUE ~ as.integer(quantity)
  ))

table(is.na(ap$quantity_num))
head(ap[is.na(ap$quantity_num),])
tail(ap[is.na(ap$quantity_num),])


#write table
fwrite(ap,"~/data/ap_scripts.tsv",quote = F, row.names = F,sep="\t")


#### BREAST CANCER ####

#get atc and substance name
bc_atc <- c("L02BA01", "L02BG03", "L02BG04", "L02BG06") 
bc_sub <- atc_map[atc_map$code %in% bc_atc,]$desc %>% tolower() %>% as.character()
bc_sub_str <- paste(bc_sub,collapse = "|")

#get commercial name from lookup table
bc_lkp <- bnf_lkp %>% filter(grepl(bc_sub_str,tolower(BNF_Chemical_Substance)))
bc_com <- bc_lkp$BNF_Product %>% unique %>% tolower() %>% as.character() 

#merge commercial and substance names
bc_names <- unique(c(bc_com,bc_sub))
bc_names <- bc_names[! bc_names %in% c("nu-cross tamoxifen cit","tamoxifen cit")]
bc_names_str <- paste(bc_names,collapse = "|")

#get doac prescriptions
bc<- scripts %>%  filter(grepl("^0803|^0605",bnf_code)) %>% filter(grepl(bc_names_str,tolower(drug_name)))
length(bc$eid)
#129 150

# N s.u. extracted 
bc$eid %>% 
  unique() %>%
  length()
# 4 269 s.u.


bc <- bc %>%
  mutate(quantity_num = case_when(
    # e.g. "28 tablet", "56 tabs"
    grepl('tablet|tab|cap', tolower(quantity)) &
      !grepl('pack', tolower(quantity)) &
      !grepl('\\*', tolower(quantity)) &
      !grepl('x', tolower(quantity)) &
      !grepl('*', tolower(quantity)) ~ tablet(quantity),
    grepl('tablet|tab|cap', tolower(quantity)) &
      (grepl('pack', tolower(quantity)) |
         grepl('\\*', tolower(quantity)) |
         grepl('x', tolower(quantity)) |
         grepl('*', tolower(quantity)) ) ~ pack_per_tablet(quantity),
    TRUE ~ as.integer(quantity)
  ))


table(is.na(bc$quantity_num))
head(bc[is.na(bc$quantity_num),])


#write table
fwrite(bc,"~/data/bc_scripts.tsv",quote = F,row.names = F,sep="\t")




#### DOAC ####

#get atc and substance name
doac_atc <- c("B01AE07", "B01AF01", "B01AF02", "B01AF03") 
doac_sub <- tolower(as.character(atc_map[atc_map$code %in% doac_atc,]$desc))
doac_sub <- gsub( " .*$", "", doac_sub)
doac_sub_str <- paste(doac_sub,collapse = "|")


#get commercial name from lookup table
doac_lkp <- bnf_lkp %>% filter(grepl(doac_sub_str,tolower(BNF_Chemical_Substance)))
doac_com <- unique(doac_lkp$BNF_Product) %>% tolower() %>% as.character() 

#merge commercial and substance names
doac_names <- unique(c(doac_com,doac_sub))
doac_names_str <- paste(doac_names,collapse = "|")

#get doac prescriptions
doac <- scripts %>%
  filter(grepl("^0208",bnf_code)) %>%
  filter(grepl(doac_names_str,tolower(drug_name)))

# N prescription using both commercial and chemical substance names
doac$eid %>%
  length()
# 27 271

# N s.u. extracted using BNF
doac$eid %>% 
  unique() %>%
  length()
# 2 648 s.u.



# extract all statins prescriptions and transform quantities to integer
doac <- doac %>%
  mutate(quantity_num = case_when(
    # e.g. "28 tablet", "56 tabs"
    grepl('tablet|tab|cap', tolower(quantity)) &
      !grepl('pack', tolower(quantity)) &
      !grepl('\\*', tolower(quantity)) &
      !grepl('x', tolower(quantity)) &
      !grepl('*', tolower(quantity)) ~ tablet(quantity),
    grepl('tablet|tab|cap', tolower(quantity)) &
      (grepl('pack', tolower(quantity)) |
         grepl('\\*', tolower(quantity)) |
         grepl('x', tolower(quantity)) |
         grepl('*', tolower(quantity)) ) ~ pack_per_tablet(quantity),
    TRUE ~ as.integer(quantity)
  ))

table(is.na(doac$quantity_num))
tail(doac[is.na(doac$quantity_num),])


#write table
fwrite(doac,"~/data/doac_scripts.tsv",quote = F,row.names = F,sep="\t")


#### GLAUCOMA ####

#bnf must start with 11
mappa %>% filter(grepl("^S01E(E|D)",mappa$code))

#get atc and substance name
gla_atc <- c("S01ED01", "S01ED02", "S01ED51",
             "S01EE01", "S01EE03", "S01EE04","S01EE05") 
gla_sub <- tolower(as.character(atc_map[atc_map$code %in% gla_atc,]$desc))
gla_sub <- gla_sub[gla_sub!="timolol, combinations"]
gla_sub_str <- paste(gla_sub,collapse = "|")

#get commercial name from lookup table
gla_lkp <- bnf_lkp %>% filter(grepl(gla_sub_str,tolower(BNF_Chemical_Substance)))
gla_com <- unique(gla_lkp$BNF_Product) %>% tolower() %>% as.character()  %>% sort()
gla_com <- gla_com[!grepl("/| ",gla_com)]
gla_com <- append(gla_com,c("brinzolamide","dorzolamide","brimonidine"))



#merge commercial and substance names
gla_names <- unique(c(gla_com,gla_sub))
gla_names_str <- paste(gla_names,collapse = "|")

#get glaucoma prescriptions
gla <- scripts %>% filter(grepl("^1106",bnf_code)) %>%  filter(grepl(gla_names_str,tolower(drug_name)))


# N prescription using both commercial and chemical substance names
gla$eid %>%
  length()
# 311635

# N s.u. extracted using BNF
gla$eid %>% 
  unique() %>%
  length()
# 5374 s.u.


#write table
fwrite(gla,"~/data/gla_scripts.tsv",quote = F, row.names = F, sep= "\t")
