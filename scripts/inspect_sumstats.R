rm(list=ls())

library(ggplot2)
library(data.table)
library(dplyr)

gwas <- fread('/home/cordioli/drugs/results/breast_canc.pheweb.gz')
gws_snps <- gwas %>%
  filter(pval <= 5e-08)
