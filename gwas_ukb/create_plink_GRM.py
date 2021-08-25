import hail as hl
hl.init(min_block_size=128)

mfi_table_path = "gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.mfi.ht"
samples_table_path = "gs://ukb31063/ukb31063.neale_gwas_samples.both_sexes.ht"
r2 = 0.05
window = 1e7

samples = hl.read_table(samples_table_path)

mfi = hl.read_table(mfi_table_path)
mfi_filt = mfi.filter((mfi.info > 0.95) & (mfi.maf > 0.01)).key_by('locus', 'alleles')

mfi_filt.describe()

# Import bgens, importing only qc'ed variants
mt = hl.import_bgen("gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,"
                    "14,15,16,17,18,19,20,21,22}_v3.bgen",
                    sample_file="gs://ukb31063/ukb31063.autosomes.sample",
                    entry_fields=['GT'],
                    variants=mfi_filt)

# Filter samples
mt = mt.filter_cols(hl.is_defined(samples[mt.col_key]))

# Remove chr8 inversion and HLA
mt = mt.filter_rows(~hl.parse_locus_interval('8:8055789-11980649').contains(mt.locus) &
                    ~hl.parse_locus_interval('6:28477797-33448354').contains(mt.locus))

ht = hl.ld_prune(mt.GT, r2=float(r2), bp_window_size=int(float(window)))
mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
mt = mt.filter_rows(hl.rand_bool(0.55))

hl.export_plink(mt, "gs://mattia-ukb/plink_grm/ukb.pruned.r0.05.w1e7.GRM")
