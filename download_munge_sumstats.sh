#!/bin/bash

readarray -t traits < /home/cordioli/drugs/data/adherence_phenolist_breast_canc.txt
job_id=cfc341f1-2fba-4d31-b07a-a4e52a57f27b

conda activate ldsc

for t in "${traits[@]}"
do

gsutil -m cp gs://dsge-cromwell/saige/${job_id}/call-test_combine/shard-*/sub.test_combine/*/call-combine/${t}.gz /home/cordioli/drugs/results/
gsutil -m cp gs://dsge-cromwell/saige/${job_id}/call-test_combine/shard-*/sub.test_combine/*/call-combine/attempt-*/${t}.gz /home/cordioli/drugs/results/

gsutil -m cp gs://dsge-cromwell/saige/${job_id}/call-test_combine/shard-*/sub.test_combine/*/call-combine/${t}.pheweb.gz /home/cordioli/drugs/results/
gsutil -m cp gs://dsge-cromwell/saige/${job_id}/call-test_combine/shard-*/sub.test_combine/*/call-combine/attempt-*/${t}.pheweb.gz /home/cordioli/drugs/results/

# unzip and change SNP id column
zcat /home/cordioli/drugs/results/${t}.gz | awk 'BEGIN{FS=OFS="\t"} {gsub("chr", "", $3); gsub("_[ATGC]*_[ATGC]$", "", $3)} 1' > /home/cordioli/drugs/results/${t}

# separate header
head -n 1 /home/cordioli/drugs/results/${t} > /home/cordioli/drugs/results/header
tail -n +2 /home/cordioli/drugs/results/${t} > tmp && mv tmp /home/cordioli/drugs/results/${t}

# sort and merge with RSID
sort -k 3 /home/cordioli/drugs/results/${t} > tmp && mv tmp /home/cordioli/drugs/results/${t}
join -1 3 -2 2 -t $'\t' -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 2.1 /home/cordioli/drugs/results/${t} /home/cordioli/finngen_rsids.tsv > tmp && mv tmp /home/cordioli/drugs/results/${t}

# re-attach header
awk '{print $0, "\tRSID"}' /home/cordioli/drugs/results/header > tmp && mv tmp /home/cordioli/drugs/results/header
cat /home/cordioli/drugs/results/header /home/cordioli/drugs/results/${t} > tmp && mv tmp /home/cordioli/drugs/results/${t}
rm /home/cordioli/drugs/results/header

# remove columns 3,4 (rsid, SNPID)
cut -f-2,5-21 /home/cordioli/drugs/results/${t} > tmp && mv tmp /home/cordioli/drugs/results/${t}

# Munge summary stats
/home/cordioli/ldsc/munge_sumstats.py \
--sumstats /home/cordioli/drugs/results/${t} \
--out /home/cordioli/drugs/results/${t} \
--snp RSID \
--a2 Allele1 \
--a1 Allele2 \
--p p.value \
--signed-sumstats BETA,0 \
--merge-alleles /home/cordioli/ldsc/w_hm3.snplist
done