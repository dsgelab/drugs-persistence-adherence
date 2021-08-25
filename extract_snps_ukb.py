#!/usr/bin/env python
# coding: utf-8

##Mattia's script for extracting FinnGen genotypes
import hail as hl
hl.init()

vcf_list = ['gs://finngen-production-library-red/finngen_R7/genotype_2.0/data/finngen_R7_chr%s.vcf.gz' % chrom for chrom in [4, 5, 12]]
mt = hl.import_vcf(vcf_list, force_bgz=True, reference_genome='GRCh38')

loci = ["12:21178615", "12:21130388", "4:88131171", "5:75359673", "5:75355259"]
loci_to_extract = [hl.parse_locus('chr'+s, reference_genome='GRCh38') for s in loci]

mt_f = mt.filter_rows(hl.literal(loci_to_extract).contains(mt.locus))
mt_f = mt_f.key_rows_by('locus', 'alleles', 'rsid')

mt_f.GT.n_alt_alleles().export('gs://mattia/pharmacogen_statins_25052021_GT.tsv')
mt_f.DS.export('gs://mattia/pharmacogen_statins_25052021_DS.tsv')
mt_f.rows().info.INFO.export('gs://mattia/pharmacogen_statins_25052021_INFO.tsv')
