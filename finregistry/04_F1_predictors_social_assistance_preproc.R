rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyr)


# Individual level social assistance data
soc <- fread('/data/processed_data/thl_social_assistance/3214_FinRegistry_toitu_MattssonHannele07122020.csv.finreg_IDsp')

# We want to model social aid as "ever received any social support in the year prior start of treatment"
soc <- soc %>%
    select(FINREGISTRYID = TNRO, YEAR = TILASTOVUOSI, TAMMI:JOULU) 

soc_long <- soc %>%
    pivot_longer(cols = TAMMI:JOULU,
                 names_to = "month",
                 values_to = "received") %>% 
    mutate(MONTH = case_when(month == "TAMMI" ~ "01",
                             month == "HELMI" ~ "02",
                             month == "MAALIS" ~ "03",
                             month == "HUHTI" ~ "04",
                             month == "TOUKO" ~ "05",
                             month == "KESA" ~ "06",
                             month == "HEINA" ~ "07",
                             month == "ELO" ~ "08",
                             month == "SYYS" ~ "09",
                             month == "LOKA" ~ "10",
                             month == "MARRAS" ~ "11",
                             month == "JOULU" ~ "12"),
           date = as.Date(paste0(YEAR, "-", MONTH, "-01"), "%Y-%m-%d"))

saveRDS(soc_long, '/data/projects/project_adherence/data/soc_assistance_date.rds')
