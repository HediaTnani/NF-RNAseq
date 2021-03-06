#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process star_idx {
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    path genome
    path gtf
    
    output:
    path "index_dir/", emit: index

    script:
    """
    mkdir index_dir
    
    STAR --runThreadN $params.threads \\
      --runMode genomeGenerate \\
      --genomeDir index_dir/ \\
      --genomeFastaFiles $genome \\
      --genomeSAindexNbases $params.genomeSAindexNbases \\
      --sjdbGTFfile $gtf
    """
}

process star_alignment {
    publishDir "$params.outdir/mapped-reads/", mode: 'copy', overwrite: true  //, pattern: "*.bam"  
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    // (trim_fq, IDX.out, gtf)
    tuple val(sample), path(reads) 
    path indexDir
    path gtf

    output:
    tuple val(sample), path("*.bam"), emit: bam

    script:
    def data = params.paired ? "${reads[0]} ${reads[1]}" : "${reads}"
    """
    STAR  \\
        --readFilesIn ${data} \\
        --runThreadN $task.cpus \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix $sample. \\
        --genomeDir ${indexDir}
    """
}


