#!/usr/bin/env nextflow

// author: Hédia Tnani (Pasteur Institute of Tunis)

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

println """\
      LIST OF PARAMETERS
================================
            GENERAL
Data-folder      : $params.datadir
Results-folder   : $params.outdir
================================
      INPUT & REFERENCES 
Input-files      : $params.reads
Reference genome : $params.genome
GTF-file         : $params.gtf
================================
          TRIMMOMATIC
Sliding window   : $params.slidingwindow
Average quality  : $params.avgqual
================================
             STAR
Threads          : $params.threads
Length-reads     : $params.lengthreads
SAindexNbases    : $params.genomeSAindexNbases
================================
              METADATA
Metadata         : $params.metadata
================================ 
"""


// We start by creating a channel 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

genome = file(params.genome)
gtf = file(params.gtf)
metadata = file(params.metadata)

// We include the modules

include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "${launchDir}/modules/fastqc" 
include { trimmomatic } from "${launchDir}/modules/trimmomatic"
include { star_idx; star_alignment } from "${launchDir}/modules/star"
include { multiqc } from "${launchDir}/modules/multiqc"
include {featurecounts} from "${launchDir}/modules/featurecounts2" 
include {finalcount} from "${launchDir}/modules/merge"
include {countmatrix} from "${launchDir}/modules/countmatrix"
//include {deseq} from "${launchDir}/modules/deseq2"

// We run the workflow  
workflow {
  // QC on raw reads
  read_pairs_ch.view()
  fastqc_raw(read_pairs_ch) 
	
  // Trimming & QC
  trimmomatic(read_pairs_ch)
  fastqc_trim(trimmomatic.out.trim_fq)
	
  // Mapping
  star_idx(genome, gtf)
  star_alignment(trimmomatic.out.trim_fq, star_idx.out.index, gtf)
  

  // Multi QC on all results
  multiqc((fastqc_raw.out.fastqc_out).mix(fastqc_trim.out.fastqc_out).collect())
  
  // Sorting BAM
  featurecounts(star_alignment.out.bam,gtf)
  star_alignment.out.bam.view()

  // Clean counts
  finalcount(featurecounts.out.counts)
  
  // countmatrix
  countmatrix(finalcount.out.mergedcounts.collect())
  
 // deseq 
 //deseq(countmatrix.out.countmatr, metadata)
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Time to complete workflow execution: $workflow.duration"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: $workflow.errorMessage"
}

