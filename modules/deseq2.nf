#!/usr/bin/env nextflow
// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process deseq {
  label 'high'
  publishDir "$params.outdir/output_deseq2", mode: 'copy', overwrite: true
  container "quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0"

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

