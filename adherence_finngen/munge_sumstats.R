rm(list=ls())

library(data.table)
library(dplyr)

setwd('~/drugs/results/R8_20220112/sumstats_stop/')

sumstats <- system("ls *.af_0.01.info_0.6.gz", intern = TRUE)

rsids <- fread('~/finngen_rsids.tsv') %>% 
  select(rsid = V1, chr_pos = V2)

for (s in sumstats) {
  d <- fread(s)
  d %>% 
    mutate(chr_pos = paste(ifelse(`#chrom` == 23, "X", `#chrom`), pos, sep = "_")) %>% 
    left_join(rsids) %>% 
    select(rsid, `#chrom`, pos, ref, alt, pval, beta, sebeta, af_alt, info, n) %>% 
    fwrite(paste0(s,'.rsid'), sep = "\t", quote = F)
  
  system(paste0('bgzip ',s,'.rsid'))
}

sumstats <- system("ls *.rsid.gz", intern = TRUE)

system("mkdir ../munged")

# cd /home/ivm/drugs/results/R8_20220112/sumstats/
# for s in *rsid.gz
# do
# python2 /usr/local/ldsc/munge_sumstats.py \
# --sumstats ${s} \
# --out ../munged/${s} \
# --snp rsid \
# --a2 ref \
# --a1 alt \
# --p pval \
# --signed-sumstats beta,0 \
# --merge-alleles /home/ivm/w_hm3.snplist
# done

cd /home/ivm/drugs/results/R8_20220112/sumstats_stop/
for s in *rsid.gz
do
python2 /usr/local/ldsc/munge_sumstats.py \
--sumstats ${s} \
--out ../munged_stop/${s} \
--snp rsid \
--a2 ref \
--a1 alt \
--p pval \
--signed-sumstats beta,0 \
--merge-alleles /home/ivm/w_hm3.snplist
done


  
# /home/cordioli/ldsc/munge_sumstats.py \
# --sumstats ukbb_statins.AF0.001.gz \
# --out ukbb_statins.AF0.001.gz \
# --snp rsid \
# --a2 A0 \
# --a1 A1 \
# --p p_value \
# --signed-sumstats beta,0 \
# --merge-alleles /home/cordioli/ldsc/w_hm3.snplist


# munge sumstats optional components

cd /home/ivm/drugs/sumstats_for_rg/

python2 /usr/local/ldsc/munge_sumstats.py \
--sumstats GCST90012790_buildGRCh38.tsv.gz \
--out processed/GCST90012790_buildGRCh38.tsv.gz \
--snp variant_id \
--a2 other_allele \
--a1 effect_allele \
--p p_value \
--signed-sumstats beta,0 \
--N 300639 \
--merge-alleles /home/ivm/w_hm3.snplist

python2 /usr/local/ldsc/munge_sumstats.py \
--sumstats GCST90012791_buildGRCh38.tsv.gz \
--out processed/GCST90012791_buildGRCh38.tsv.gz \
--snp variant_id \
--a2 other_allele \
--a1 effect_allele \
--p p_value \
--signed-sumstats beta,0 \
--N-cas 96035 \
--N-con 119092 \
--merge-alleles /home/ivm/w_hm3.snplist

python2 /usr/local/ldsc/munge_sumstats.py \
--sumstats GCST90012792_buildGRCh38.tsv.gz \
--out processed/GCST90012792_buildGRCh38.tsv.gz \
--snp variant_id \
--a2 other_allele \
--a1 effect_allele \
--p p_value \
--signed-sumstats beta,0 \
--N-cas 146074 \
--N-con 148713 \
--merge-alleles /home/ivm/w_hm3.snplist

python2 /usr/local/ldsc/munge_sumstats.py \
--sumstats GCST90012793_buildGRCh38.tsv.gz \
--out processed/GCST90012793_buildGRCh38.tsv.gz \
--snp variant_id \
--a2 other_allele \
--a1 effect_allele \
--p p_value \
--signed-sumstats beta,0 \
--N-cas 361501 \
--N-con 89535 \
--merge-alleles /home/ivm/w_hm3.snplist

python2 /usr/local/ldsc/munge_sumstats.py \
--sumstats GCST90012794_buildGRCh38.tsv.gz \
--out processed/GCST90012794_buildGRCh38.tsv.gz \
--snp variant_id \
--a2 other_allele \
--a1 effect_allele \
--p p_value \
--signed-sumstats beta,0 \
--N-cas 336633 \
--N-con 114464 \
--merge-alleles /home/ivm/w_hm3.snplist