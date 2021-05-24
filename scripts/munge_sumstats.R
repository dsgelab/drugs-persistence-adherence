rm(list=ls())

library(data.table)
library(dplyr)

setwd('~/drugs/results/20210422_R7/')

sumstats <- system("ls *.gz", intern = TRUE)

rsids <- fread('~/finngen_rsids.tsv') %>% 
  select(rsid = V1, chr_pos = V2)

for (s in sumstats) {
  d <- fread(s)
  d %>% 
    mutate(chr_pos = paste(ifelse(`#chrom` == 23, "X", `#chrom`), pos, sep = "_")) %>% 
    left_join(rsids) %>% 
    select(rsid, `#chrom`, pos, ref, alt, pval, beta, sebeta, af_alt, info, n) %>% 
    fwrite(paste0(s,'.rsid.txt'), sep = "\t", quote = F)
  
  system(paste0('bgzip ',s,'.rsid.txt'))
}

# # # # 
# conda activate ldsc
# 
# for s in *.rsid.txt.gz
# do
# /home/cordioli/ldsc/munge_sumstats.py \
# --sumstats /home/cordioli/drugs/results/20210422_R7/${s} \
# --out /home/cordioli/drugs/results/20210422_R7/munged/${s} \
# --snp rsid \
# --a2 ref \
# --a1 alt \
# --p pval \
# --signed-sumstats beta,0 \
# --merge-alleles /home/cordioli/ldsc/w_hm3.snplist
# done

# /home/cordioli/ldsc/munge_sumstats.py \
# --sumstats ukbb_statins.AF0.001.gz \
# --out ukbb_statins.AF0.001.gz \
# --snp rsid \
# --a2 A0 \
# --a1 A1 \
# --p p_value \
# --signed-sumstats beta,0 \
# --merge-alleles /home/cordioli/ldsc/w_hm3.snplist
