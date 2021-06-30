# NF-RNAseq

The pipeline is built using Nextflow [Nextflow](https://www.nextflow.io) which a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers. This workflow uses the Nextflow DSL2 implementation. Singularity containers will be addeded.

## Introduction

**NF-RNAseq** is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. 

## Installation and requirements

# Nextflow

To run this pipeline you need to have Nextflow installed. Nextflow can be used on any POSIX compatible system (Linux, OS X, etc). It requires Bash 3.2 (or later) and Java 8 (or later, up to 15) to be installed.

It only needs two easy steps:

Download the executable package by copying and pasting either one of the following commands in your terminal window: `wget -qO- https://get.nextflow.io | bash`

Or, if you prefer curl: `curl -s https://get.nextflow.io | bash`

This will create the nextflow main executable file in the current directory.

Make the binary executable on your system by running `chmod +x nextflow`.

# Docker

**Run Docker Commands Without Sudo**
It is advisable to keep the settings as is. However, you can bypass typing sudo every time. Adding the user to the docker group grants privileges equivalent to root.

1. First, create the docker group with the command:`sudo groupadd docker`

2. Then, type the following command (making sure to replace [user] with your username): `sudo usermod -aG docker [user]`

3. Enable the new settings with: `su - [user]`

# Singularity

TODO


## Pipeline summary

This pipeline has been run with test data. The steps 

1. QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment[`STAR`](https://github.com/alexdobin/STAR) 
4. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
5. Quantification ([`FeatureCounts`](http://subread.sourceforge.net/)
6. Merge counts (Rscript)

# Pipeline parameters
The Nextflow script defines a pipeline parameter `params.input`. Pipeline parameters enable you to change the input to the workflow at runtime, via the command line or a configuration file, so they are not hard-coded into the script.

Pipeline parameters are declared in the workflow by prepending the prefix params, separated by dot character, to a variable name e.g., params.input. Their value can be specified on the command line by prefixing the parameter name with a double dash character, e.g., --input.

` nextflow run RNAseq.nf --input 'data/reads/*R{1,2}.fq.gz'`

# Modules

Each step is separated in a module with an input and an output.

Nextflow (DSL2) allows the definition of module scripts that can be included and shared across workflow pipelines.

A module file is a Nextflow script containing one or more process definitions that can be imported from another Nextflow script.



## Run the pipeline

To run a Nextflow script use the command `nextflow run <script_name>`.


