#!/usr/bin/env nextflow
// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process deseq {
  label 'high'
  publishDir "$params.outdir/output_deseq2", mode: 'copy', overwrite: true
  container "hediatnani/nf-renv:d57d4d6"

input:
    path("countmatrix.tsv")
    path metadata
output:
    path("*_results.txt"), emit: degs

script:
"""
Rscript --vanilla  $launchDir/bin/deseq.R ${metadata}
"""
}

