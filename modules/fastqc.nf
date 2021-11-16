#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

process fastqc {
  publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
  label 'low'
  container 'quay.io/biocontainers/fastqc:0.11.9--0'
  
  input:
  tuple val(sample), path(reads)

  output:
  path("*_fastqc.{zip,html}"), emit: fastqc_out

  script:
  def data = params.paired ? "${reads[0]} ${reads[1]}" : "${reads}"
  """
  fastqc ${data}
  """
}
 
