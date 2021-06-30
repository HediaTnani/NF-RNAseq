#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process finalcount {
  tag "$sample"
  label 'high'
  publishDir "$params.outdir/mergedCounts", mode: 'copy', overwrite: true
  
 

  input:
    tuple val(sample), path("${sample}.counts.txt")


  output:
    path("${sample}.counts.tsv"), emit: mergedcounts

  script:
"""

sed '1d' ${sample}.counts.txt | cut -f1,7 > ${sample}.counts.tsv
"""
}


