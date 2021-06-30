#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process featurecounts {
  tag "$sample"
  label 'high'
  publishDir "$params.outdir/counts", mode: 'copy', overwrite: true
  container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"


  input:
    tuple val(sample), path(bam)
    path gtf


  output:
    tuple val(sample), path("${sample}.counts.txt"), emit: counts 
    path "${sample}.counts.txt.summary", emit: summary

  script:
"""
featureCounts \\
        -T $task.cpus \\
        -t exon -g gene_id \\
        -a ${gtf} \\
        -o ${sample}.counts.txt \\
        ./*.bam
"""
}

