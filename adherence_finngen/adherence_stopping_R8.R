rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(dplyr)

source('drugs/scripts/adherence_funs.R')

purch <- fread('/home/ivm/drugs/data/finngen_R8_purch_vnr_98.gz')

# # # For each medication, get complete trajectories (i.e. without filtering for min N of purchases) to then calculate stopping

# # # # # # # #
#   STATINS   #
# # # # # # # #
st <- getTrajectoriesStopping(purch,'^C10AA')
fwrite(st, 'drugs/data/R8_statins_stopping.txt', sep = '\t', quote = F)


# # # # # # # # # # # # # # # # # #  
#   BLOOD PRESSURE MEDICATIONS    #
# # # # # # # # # # # # # # # # # #
bp <- getTrajectoriesStopping(purch,'^C0(2|3|8|9)')
fwrite(bp, 'drugs/data/R8_blood_pressure_stopping.txt', sep = '\t', quote = F)


# # # # # # # # # # # # # # # # # #  
#   CLOPIDOGREL & ASA+DIP   #
# # # # # # # # # # # # # # # # # #
# clopi '^B01AC04'
ATCs <- data.frame(atc = c("B01AC04","B01AC30"),
                   dose = c(1,2),
                   stringsAsFactors = F)

cd <- getTrajectoriesStopping2(purch,ATCs)
fwrite(cd, 'drugs/data/R8_clopi_dipy_stopping.txt', sep = '\t', quote = F)


# # # # # # # # # # #
#   BREAST CANCER   #
# # # # # # # # # # #
bc <- getTrajectoriesStopping(purch,'^L02B(A01|G04|G06|G03)')
fwrite(bc, 'drugs/data/R8_breast_cancer_stopping.txt', sep = '\t', quote = F)


# # # # # # # #
#   DOAC      #
# # # # # # # #
ATCs <- data.frame(atc = c("B01AF01","B01AF02","B01AF03","B01AE07"),
                   dose = c(1,2,1,2),
                   stringsAsFactors = F)

do <- getTrajectoriesStopping2(purch,ATCs)
fwrite(do, 'drugs/data/R8_doac_stopping.txt', sep = '\t', quote = F)


# # # # # # # #
#   GLAUCOMA  #
# # # # # # # #

# # # Some pkokoNum need manual check becasue they are wrong
gl <- purch %>%
  filter(grepl('^S01E(E|D)', CODE1))

pkoko_list <- c("1X2,5 ML", "3 x 3 ml", "3X5ML", "30 x 0.4 g", "90 x 0.3 ml", "30 x 0.2 ml", "90 x 0.2 ml", "90 (18 x 5) x 0.2 ml", "3 x 2.5 ml")

getPkokoNum <- function(pkoko) {
  b <- strsplit(gsub(" |ml|g","",tolower(gsub(",", "\\.", pkoko))), "x")
  return( unlist(lapply(b, function(i) as.numeric(i[1])*as.numeric(i[2]))) )
}

gl$pkoko_num[gl$pkoko %in% pkoko_list] <- getPkokoNum(gl$pkoko[gl$pkoko %in% pkoko_list])

pkoko <- gl %>%
  ungroup() %>%
  select(pkoko, pkoko_num) %>%
  distinct

fwrite(gl, 'drugs/data/R8_glaucoma_all_pkoko_ok.txt', sep = '\t', quote = F)

g <- fread('drugs/data/R8_glaucoma_all_pkoko_ok.txt')
gl <- getTrajectoriesStopping(g,'^S01E(E|D)',0.1)
fwrite(gl, 'drugs/data/R8_glauc_stopping.txt', sep = '\t', quote = F)

getmode <- function(v) {
  uniqv <- unique(v[!is.na(v)])
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(gl$pills_norm)


