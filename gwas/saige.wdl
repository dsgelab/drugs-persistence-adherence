import "saige_sub.wdl" as sub
import "saige_summary.wdl" as summary

task null {

    String pheno
    File bedfile
    File bimfile = sub(bedfile, "\\.bed$", ".bim")
    File famfile = sub(bedfile, "\\.bed$", ".fam")
    String prefix = basename(bedfile, ".bed") + "-" + pheno
    File phenofile
    String covariates
    String covariates_age_first = covariates + ",age_first_" + pheno
    String sampleID
    String traitType
    String docker
    String loco
    Float traceCVcutoff
    Float ratioCVcutoff
    Int minCovariateCount
    Int cpu = 32
    Boolean invNormalize=false

    command {

        # continuous traits don't have this file and optional outputs are not currently supported
        echo "placeholder" > ${prefix}"_30markers.SAIGE.results.txt"

        step1_fitNULLGLMM.R \
            --plinkFile=${sub(bedfile, "\\.bed$", "")} \
            --phenoFile=${phenofile} \
            --phenoCol=${pheno} \
            --covarColList=${covariates_age_first} \
            --sampleIDColinphenoFile=${sampleID} \
            --traitType=${traitType} \
            --outputPrefix=${prefix} \
            --nThreads=${cpu} \
            --LOCO=${loco} \
            --traceCVcutoff ${traceCVcutoff} \
            --ratioCVcutoff ${ratioCVcutoff} \
            --minCovariateCount ${minCovariateCount} \
            ${true='--invNormalize=TRUE' false=' ' invNormalize}
    }

    output {

        File modelfile = prefix + ".rda"
        File varianceratio = prefix + ".varianceRatio.txt"
        File randommarkers = prefix + "_30markers.SAIGE.results.txt"
    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "7 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow saige {

    String docker
    File phenolistfile
    String traitType
    Array[String] phenos = read_lines(phenolistfile)
    String loco
    String analysisType

    scatter (pheno in phenos) {

        call null {
            input: docker=docker, pheno=pheno, traitType=traitType, loco=loco
        }

        call sub.test_combine {
            input: docker=docker, pheno=pheno, traitType=traitType, nullfile=null.modelfile, loco=loco, analysisType=analysisType
        }

        call summary.summary{
            input: 
                pheno=pheno,
                input_file=test_combine.out
        }
    }
}
