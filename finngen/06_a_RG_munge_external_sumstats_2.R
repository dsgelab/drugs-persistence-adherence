rm(list = ls())
library(data.table)
library(dplyr)

left.1 <- c("biomarkers-30750-both_sexes-irnt.tsv.gz",
            "biomarkers-30780-both_sexes-irnt.tsv.gz",
            "biomarkers-30870-both_sexes-irnt.tsv.gz",
            "sumstats_neuroticism_ctg_format_formatted.txt.gz",
            "continuous-LDLC-both_sexes-medadj_irnt.tsv.gz",
            "UKB-ICBPmeta750k_SBPsummaryResults.txt.gz",
            "UKB-ICBPmeta750k_DBPsummaryResults.txt.gz",
            "METAANALYSIS_DIAGRAM_SE1.txt.gz",
            "CIMBA_BRCA1_BCAC_TN_meta_summary_level_statistics.txt.gz")

s_list <- fread('/home/ivm/drugs/data/genetic_correlations/sumstats_for_rg.tsv')

s_list_1 <- s_list %>%
  filter(V1 %in% left.1) %>% 
  mutate(V7 = "rsid")

left.2 <- c("meta_v3_onco_euro_overall_ChrAll_1_release.txt.gz",
            "lifegen_phase2_bothpl_alldr_2017_09_18.tsv.gz",
            "Meta-analysis_Locke_et_al+UKBiobank_2018_bmi.txt.gz")

s_list_2 <- s_list %>% 
  filter(V1 %in% left.2)

fwrite(bind_rows(s_list_1, s_list_2), '/home/ivm/drugs/data/genetic_correlations/sumstats_for_rg_round2.tsv', sep ="\t", na = "NA", quote = F, col.names = F)

readarray stats < /home/ivm/drugs/data/genetic_correlations/sumstats_for_rg_round2.tsv
for row in "${stats[11]}";do      
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

echo ${f}

if [ "${nca}" =  "NA" ]
then
  python2 /usr/local/ldsc/munge_sumstats.py \
  --sumstats /home/ivm/drugs/data/genetic_correlations/liftover/sumstats/${f} \
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
  --sumstats /home/ivm/drugs/data/genetic_correlations/liftover/sumstats/${f} \
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