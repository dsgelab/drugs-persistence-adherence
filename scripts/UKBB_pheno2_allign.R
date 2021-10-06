doac <- fread("~/data/doac_pheno.tsv")
ap <- fread("~/data/ap_pheno.tsv")
bc <- fread("~/data/bc_pheno.tsv")
bp <- fread("~/data/bp_pheno.tsv")
gla <- fread("~/data/gla_pheno.tsv")
stat <- fread("~/data/st_pheno.tsv")

cov <- fread("~/data/all_pops_non_eur_pruned_within_pop_pc_covs.tsv")
sam <- fread("~/genotype/bgen.samples.txt")
bgenX_sam <- fread("~/data/samples_chrX_missing.txt")




doac_phen <- doac %>%
   select(eid,age_first_purch,isgood,adherence,adherence_std)  %>%
   merge(cov,by.x="eid",by.y="s") %>%
   filter(eid %in% sam$V1,
          !eid %in% bgenX_sam$V1) %>%
   mutate(adherence_std = as.numeric(scale(adherence)))
 
ap_phen <- ap %>%
   select(eid,age_first_purch,isgood,adherence,adherence_std)  %>%
   merge(cov,by.x="eid",by.y="s") %>%
   filter(eid %in% sam$V1,
          !eid %in% bgenX_sam$V1) %>%
   mutate(adherence_std = as.numeric(scale(adherence)))

bc_phen <- bc %>%
   select(eid,age_first_purch,isgood,adherence,adherence_std)  %>%
   merge(cov,by.x="eid",by.y="s")%>%
   filter(eid %in% sam$V1,
          !eid %in% bgenX_sam$V1) %>%
   mutate(adherence_std = as.numeric(scale(adherence)))

bp_phen <- bp %>%
   select(eid,age_first_purch,isgood,adherence,adherence_std)  %>%
   merge(cov,by.x="eid",by.y="s") %>%
   filter(eid %in% sam$V1,
          !eid %in% bgenX_sam$V1) %>%
   mutate(adherence_std = as.numeric(scale(adherence)))
 
gla_phen <- gla %>%
   select(eid,age_first_purch,isgood,adherence,adherence_std)  %>%
   merge(cov,by.x="eid",by.y="s") %>%
   filter(eid %in% sam$V1,
          !eid %in% bgenX_sam$V1) %>%
   mutate(adherence_std = as.numeric(scale(adherence))) 

stat_phen <- stat %>%
   select(eid,age_first_purch,isgood,adherence,adherence_std)  %>%
   merge(cov,by.x="eid",by.y="s") %>%
   filter(eid %in% sam$V1,
          !eid %in% bgenX_sam$V1) %>%
   mutate(adherence_std = as.numeric(scale(adherence)))

id_pheno_list <- unique(c(doac_phen$eid,
                          ap_phen$eid,
                          bc_phen$eid,
                          bp_phen$eid,
                          gla_phen$eid,
                          stat_phen$eid))


abc <- id_pheno_list[id_pheno_list %in% bgenX_sam$V1]


fwrite(doac_phen,"~/data/doac_pheno_final.tsv",quote=F,sep="\t")
fwrite(ap_phen,"~/data/antiplat_pheno_final.tsv",quote=F,sep="\t")
fwrite(bc_phen,"~/data/breastcanc_pheno_final.tsv",quote=F,sep="\t")
fwrite(bp_phen,"~/data/bloodpres_pheno_final.tsv",quote=F,sep="\t")
fwrite(gla_phen,"~/data/glaucoma_pheno_final.tsv",quote=F,sep="\t")
fwrite(stat_phen,"~/data/statins_pheno_final.tsv",quote=F,sep="\t")



