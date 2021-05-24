rm(list = ls())

library(data.table)
library(dplyr)
library(tidyr)
library(qqman)

setwd('~/drugs/results/20210520_UKB/')

d <- fread('~/drugs/results/20210520_UKB/ukbb_statins.gz')
d <- fread('~/drugs/results/20210520_UKB/ukbb_statins.AF0.001.reheader.gz')


# split AF
split <- strsplit(d$AF, ",")
d$AF_A0 <- unlist( lapply( split, function(x) as.numeric(gsub("\\[","",x[1])) ))
d$AF_A1 <- unlist( lapply( split, function(x) as.numeric(gsub("\\]","",x[2])) ))

data <- d %>% 
  select(rsid, chr, pos, A0, A1, AF_A0, AF_A1, info, n, beta, standard_error, p_value) %>% 
  filter(!is.nan(p_value),
         AF_A0 > 0.001,
         AF_A1 > 0.001)

pcol <- "p_value"
bp_col <- "pos"
chr_col <- "chr"

fwrite(data, "ukbb_statins.AF0.001.gz", quote = F, na = "NA", sep = "\t", compress = "gzip")

output_prefix <- paste0("UKB_statins.AF0.001")

quants <- c(0.7,0.5,0.1,0.01, 0.001)
subdata <- data[ !is.na(data[[pcol]]) & is.numeric( data[[pcol]]  ) ]
lambda  <- round(  quantile(  (qchisq(1-subdata[[pcol]], 1) ), probs=quants ) / qchisq(quants,1), 3)

png( paste(output_prefix,"_", pcol ,"_qqplot.png", sep="" ))
qq(subdata[[pcol]])
title(c(file, paste("\n", "\nlambda ", quants, ": ", lambda, sep="" )), line=-0.2)
dev.off()

print("subsetting p-vals < 0.01 for manhattan...")
subdata <- subdata[ subdata[[pcol]]<0.01 & subdata[[pcol]]>0 ]

print( paste0("Plotting manhattan with ", nrow(subdata), " variants") )
print( summary(subdata[[pcol]] ))

png( paste(output_prefix,"_",pcol,"_manhattan.png", sep=""), width=1000, height=400)
logs <- -log10(subdata[[pcol]])
manhattan( data.table(subdata[,c(bp_col,pcol,chr_col), with=F]) , chr=chr_col, bp=bp_col, p=pcol, ylim=c( 2,max(logs)+1))
dev.off()