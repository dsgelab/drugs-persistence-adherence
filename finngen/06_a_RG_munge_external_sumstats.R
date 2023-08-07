rm(list = ls())
library(data.table)
library(dplyr)

to_use <- fread('/home/ivm/drugs/data/genetic_correlations/sumstats_list.csv', skip = 1, header = F)
all <- fread('/finngen/library-green/prs/PRS_data.txt')
all <- all %>%
  filter(filename %in% to_use$V3) %>% 
  select(filename, pheno, n_total, n_cases, n_ctrls, effect_type, variant, effect_allele, other_allele, effect, pval, build)

all$pheno[all$filename == "Meta-analysis_Locke_et_al+UKBiobank_2018_bmi.txt.gz"] <- "bmi"

fwrite(all, '/home/ivm/drugs/data/genetic_correlations/sumstats_for_rg.tsv', sep ="\t", na = "NA", quote = F, col.names = F)

# # # bash
readarray stats < /home/ivm/drugs/data/genetic_correlations/sumstats_for_rg.tsv
for row in "${stats[@]}";do      
row_array=(${row})
f=${row_array[0]}
s=${row_array[1]}
n=${row_array[2]}
nca=${row_array[3]}
nco=${row_array[4]}
eff_type=${row_array[5]}
snp=${row_array[6]}
a1=${row_array[7]}
a2=${row_array[8]}
eff=${row_array[9]}
if [ "${eff_type}" =  "BETA" ]
then
  n_val=0
else
  n_val=1  
fi
p=${row_array[10]}

if [ "${nca}" =  "NA" ]
then
  python2 /usr/local/ldsc/munge_sumstats.py \
  --sumstats /finngen/library-green/prs/sumstats/${f} \
  --out ./${s} \
  --snp ${snp} \
  --a2 ${a2} \
  --a1 ${a1} \
  --p ${p} \
  --signed-sumstats ${eff},${n_val} \
  --merge-alleles /home/ivm/w_hm3.snplist \
  --N ${n}
else
  python2 /usr/local/ldsc/munge_sumstats.py \
  --sumstats /finngen/library-green/prs/sumstats/${f} \
  --out ./${s} \
  --snp ${snp} \
  --a2 ${a2} \
  --a1 ${a1} \
  --p ${p} \
  --signed-sumstats ${eff},${n_val} \
  --merge-alleles /home/ivm/w_hm3.snplist \
  --N-cas ${nca} \
  --N-con ${nco}
fi
done