task read_bim {

    String docker
    File bimfile
    
    command {

        wc -l ${bimfile}
        wc -l ${bimfile} > n_variants.txt

        
    }

    output {

        File out = "n_variants.txt"
    }

    runtime {

        docker: "${docker}"
        cpu: "1"
        memory: "7 GB"
        disks: "local-disk 20 HDD"
        zones: "us-central1-a"
        preemptible: 2
        noAddress: true
    }
}

workflow test {

    String docker
    File bimfile
    
    call read_bim{

        input: 
            bimfile=bimfile,
            docker=docker
    
    }

}