import hail as hl
import argparse


def main(args):

    hl.init(default_reference='GRCh37',
            min_block_size=128)

    chrom_list = args.chr.split(',')

    for chrom in chrom_list:

        print("\nCHROM: {}\n".format(chrom))

        gwas_variants_path = "gs://ukb31063/ukb31063.neale_gwas_variants.ht"
        gwas_variants_mfi_path = "gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.mfi.ht"
        gwas_samples_path = "gs://ukb31063/ukb31063.neale_gwas_samples.both_sexes.ht"

        bgen_path = "gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen".format(chrom)
        bgen_samples_path = "gs://ukb31063/ukb31063.autosomes.sample"

        covariates_table_path = "gs://ukb31063/ukb31063.neale_gwas_covariates.both_sexes.ht"
        pheno_file_path = "gs://mattia-ukb/UKB_statins_summarized_after_enrolment.txt"

        out_dir = "gs://mattia-ukb/gwas_adherence/"
        out_prefix = "ukbb_statins_chr{}".format(chrom)

        # import bgen, filter only Neal Lab QC'ed variants and gwas samples
        gwas_variants = hl.read_table(gwas_variants_path)
        gwas_samples = hl.read_table(gwas_samples_path)

        mt = hl.import_bgen(bgen_path,
                            sample_file=bgen_samples_path,
                            entry_fields=['GT', 'dosage'])

        mt = mt.filter_rows(hl.is_defined(gwas_variants[mt.row_key]))
        mt = mt.filter_cols(hl.is_defined(gwas_samples[mt.col_key]))

        # annotate variants with mfi
        mfi = hl.read_table(gwas_variants_mfi_path).key_by('locus', 'alleles')
        mt = mt.annotate_rows(mfi=mfi[mt.row_key])

        # annotate phenotype and covariates
        cov = hl.read_table(covariates_table_path)
        mt = mt.annotate_cols(cov=cov[mt.col_key])

        pheno = hl.import_table(pheno_file_path,
                                delimiter='\t',
                                types={'eid': hl.tstr,
                                       'adherence_std': hl.tfloat64}).key_by('eid')
        mt = mt.annotate_cols(pheno=pheno[mt.col_key])

        mt = mt.filter_cols(hl.is_defined(mt.pheno.adherence_std))

        mt = hl.variant_qc(mt)

        gwas = hl.linear_regression_rows(
            y=mt.pheno.adherence_std,
            x=mt.dosage,
            covariates=[1.0,
                        mt.cov.isFemale,
                        mt.cov.age,
                        mt.cov.PC1,
                        mt.cov.PC2,
                        mt.cov.PC3,
                        mt.cov.PC4,
                        mt.cov.PC5,
                        mt.cov.PC6,
                        mt.cov.PC7,
                        mt.cov.PC8,
                        mt.cov.PC9,
                        mt.cov.PC10
                        ],
            pass_through=[mt.variant_qc.AF, mt.mfi.info, mt.rsid])

        res = gwas.annotate(CHR=gwas.locus.contig,
                            POS=gwas.locus.position,
                            A0=gwas.alleles[0],
                            A1=gwas.alleles[1],
                            AF_A1=gwas.AF[1]).key_by()

        res.select('rsid',
                   'CHR',
                   'POS',
                   'A0',
                   'A1',
                   'AF_A1',
                   'info',
                   'n',
                   'sum_x',
                   'y_transpose_x',
                   'beta',
                   'standard_error',
                   't_stat',
                   'p_value').export(out_dir + out_prefix + '.tsv', header=True, delimiter='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--chr', help='Comma-separated list of chromosomes to run (default 21 for testing)', default='21')
    args = parser.parse_args()

    main(args)
