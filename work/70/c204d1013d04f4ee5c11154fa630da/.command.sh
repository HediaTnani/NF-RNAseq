#!/bin/bash -ue
trimmomatic PE -threads 2         mouse_cns_E18_rep1_1.fastq.gz mouse_cns_E18_rep1_2.fastq.gz mouse_cns_E18_rep11_P.fq mouse_cns_E18_rep11_U.fq mouse_cns_E18_rep12_P.fq mouse_cns_E18_rep12_U.fq         SLIDINGWINDOW:4:15 AVGQUAL:30
