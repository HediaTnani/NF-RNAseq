#!/bin/bash -ue
trimmomatic PE -threads 2         mouse_cns_E14_rep2_1.fastq.gz mouse_cns_E14_rep2_2.fastq.gz mouse_cns_E14_rep21_P.fq mouse_cns_E14_rep21_U.fq mouse_cns_E14_rep22_P.fq mouse_cns_E14_rep22_U.fq         SLIDINGWINDOW:4:15 AVGQUAL:30
