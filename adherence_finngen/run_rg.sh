sumstats_list=$(ls -d /home/ivm/drugs/sumstats_for_rg/processed/*gz | tr '\n' ',')
sumstats_list=${sumstats_list::-1}

mkdir /home/ivm/drugs/results/R8_20220112/rg
mkdir /home/ivm/drugs/results/R8_20220112/h2

cd /home/ivm/drugs/results/R8_20220112/munged/
for t in *.sumstats.gz
do
# # # rg # # # 
python2 /usr/local/ldsc/ldsc.py \
--rg /home/ivm/drugs/results/R8_20220112/munged/${t},${sumstats_list} \
--ref-ld-chr /home/ivm/eur_w_ld_chr/ \
--w-ld-chr /home/ivm/eur_w_ld_chr/ \
--out /home/ivm/drugs/results/R8_20220112/rg/${t}_RGS

# # # h2 # # #
python2 /usr/local/ldsc/ldsc.py \
--h2 /home/ivm/drugs/results/R8_20220112/munged/${t} \
--ref-ld-chr /home/ivm/eur_w_ld_chr/ \
--w-ld-chr /home/ivm/eur_w_ld_chr/ \
--out /home/ivm/drugs/results/R8_20220112/h2/${t}_H2
done


# # # RG stopping vs traits # # #
mkdir /home/ivm/drugs/results/R8_20220112/rg_stop
mkdir /home/ivm/drugs/results/R8_20220112/h2_stop

cd /home/ivm/drugs/results/R8_20220112/munged_stop/
for t in *.sumstats.gz
do
# # # rg # # # 
python2 /usr/local/ldsc/ldsc.py \
--rg /home/ivm/drugs/results/R8_20220112/munged_stop/${t},${sumstats_list} \
--ref-ld-chr /home/ivm/eur_w_ld_chr/ \
--w-ld-chr /home/ivm/eur_w_ld_chr/ \
--out /home/ivm/drugs/results/R8_20220112/rg_stop/${t}_RGS

# # # h2 # # #
python2 /usr/local/ldsc/ldsc.py \
--h2 /home/ivm/drugs/results/R8_20220112/munged_stop/${t} \
--ref-ld-chr /home/ivm/eur_w_ld_chr/ \
--w-ld-chr /home/ivm/eur_w_ld_chr/ \
--out /home/ivm/drugs/results/R8_20220112/h2_stop/${t}_H2
done



# # # RG adherence - stop / adherence - adherence
sumstats_stop=$(ls -d /home/ivm/drugs/results/R8_20220112/munged_stop/*gz | tr '\n' ',')
sumstats_stop=${sumstats_stop::-1}

sumstats_adh=$(ls -d /home/ivm/drugs/results/R8_20220112/munged/*gz | tr '\n' ',')
sumstats_adh=${sumstats_adh::-1}

mkdir /home/ivm/drugs/results/R8_20220112/rg_adh_stop
mkdir /home/ivm/drugs/results/R8_20220112/rg_adh_adh

cd /home/ivm/drugs/results/R8_20220112/munged/
for t in *.sumstats.gz
do
# # # rg adherence - stop # # # 
python2 /usr/local/ldsc/ldsc.py \
--rg /home/ivm/drugs/results/R8_20220112/munged/${t},${sumstats_stop} \
--ref-ld-chr /home/ivm/eur_w_ld_chr/ \
--w-ld-chr /home/ivm/eur_w_ld_chr/ \
--out /home/ivm/drugs/results/R8_20220112/rg_adh_stop/${t}_RGS

# # # rg adherence - adherence # # # 
python2 /usr/local/ldsc/ldsc.py \
--rg /home/ivm/drugs/results/R8_20220112/munged/${t},${sumstats_adh} \
--ref-ld-chr /home/ivm/eur_w_ld_chr/ \
--w-ld-chr /home/ivm/eur_w_ld_chr/ \
--out /home/ivm/drugs/results/R8_20220112/rg_adh_adh/${t}_RGS
done




# # RG statins FINNGEN - UKB
# python2 /usr/local/ldsc/ldsc.py \
# --rg /home/ivm/drugs/results/R8_20220112/munged/finngen_R7_statins.af_0.001.info_0.6.txt.gz.rsid.txt.gz.sumstats.gz,ukbb_statins.AF0.001.gzsumstats.gz \
# --ref-ld-chr /home/ivm/eur_w_ld_chr/ \
# --w-ld-chr /home/ivm/eur_w_ld_chr/ \
# --out statins_finngen_ukb_RG
# 
# # RG statins UKB - all
# python2 /usr/local/ldsc/ldsc.py \
# --rg ukbb_statins.AF0.001.gz.sumstats.gz,${sumstats_list} \
# --ref-ld-chr /home/ivm/eur_w_ld_chr/ \
# --w-ld-chr /home/ivm/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/20210520_UKB/rg/statins_RGS
# 
# python2 /usr/local/ldsc/ldsc.py \
# --h2 ukbb_statins.AF0.001.gz.sumstats.gz \
# --ref-ld-chr /home/ivm/eur_w_ld_chr/ \
# --w-ld-chr /home/ivm/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/20210520_UKB/h2/statins_H2
# 
# # RG statins meta - all
# python2 /usr/local/ldsc/ldsc.py \
# --rg meta_statins.sumstats.gz,${sumstats_list} \
# --ref-ld-chr /home/ivm/eur_w_ld_chr/ \
# --w-ld-chr /home/ivm/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/20210522_meta/meta_statins_RGS
# 
# python2 /usr/local/ldsc/ldsc.py \
# --h2 meta_statins.sumstats.gz \
# --ref-ld-chr /home/ivm/eur_w_ld_chr/ \
# --w-ld-chr /home/ivm/eur_w_ld_chr/ \
# --out /home/cordioli/drugs/results/20210522_meta/meta_statins_H2