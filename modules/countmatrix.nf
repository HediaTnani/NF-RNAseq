#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process countmatrix {
  label 'high'
  publishDir "$params.outdir/countmatrix", mode: 'copy', overwrite: true
  
input:
    path("*counts.tsv")
output:
    path("countmatrix.tsv"), emit: countmatr


  script:
"""
Rscript --vanilla  $launchDir/bin/countmatrix.R
"""
}


