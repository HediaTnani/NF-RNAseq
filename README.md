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
**Installing Docker**
First, update your existing list of packages:
`sudo apt update`

Next, install a few prerequisite packages which let apt use packages over HTTPS:
`sudo apt install apt-transport-https ca-certificates curl software-properties-common`

Then add the GPG key for the official Docker repository to your system:
`curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -`

Add the Docker repository to APT sources:
`sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu focal stable"`

Next, update the package database with the Docker packages from the newly added repo:
`sudo apt update`

Make sure you are about to install from the Docker repo instead of the default Ubuntu repo:
`apt-cache policy docker-ce`

Finally, install Docker:
`sudo apt install docker-ce`

Docker should now be installed, the daemon started, and the process enabled to start on boot. Check that it’s running:
`sudo systemctl status docker`


**Run Docker Commands Without Sudo**
It is advisable to keep the settings as is. However, you can bypass typing sudo every time. Adding the user to the docker group grants privileges equivalent to root.

1. First, create the docker group with the command:`sudo groupadd docker`

2. Then, type the following command (making sure to replace [user] with your username): `sudo usermod -aG docker ${USER}`

3. Enable the new settings with: `su - ${USER}`

4. Confirm that your user is now added to the docker group by typing: `id -nG`

5. Declare that username explicitly using: `sudo usermod -aG docker username`



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

# Worflow infrastructure

**Channels**
Channels are the way in which Nextflow sends data around a workflow. Channels connect processes via their inputs and outputs. Channels can store multiple items, such as files (e.g., fastq files) or values. The number of items a channel stores determines how many times a process runs.
The `fromPath` factory method creates a queue channel emitting one or more files matching a file path. 

**Processes**
A process is the way Nextflow execute commands you would run on the command line or custom scripts. It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output channels.

**Modules**

Each step is separated in a module with an input and an output.

Nextflow (DSL2) allows the definition of module scripts that can be included and shared across workflow pipelines.

A module file is a Nextflow script containing one or more process definitions that can be imported from another Nextflow script.

**Workflow**
A workflow id combination of multiple processes.

## Run the pipeline

To run a Nextflow script use the command `nextflow run main.nf -profile docker`.


## Result architecture
# Quality control
![image](https://user-images.githubusercontent.com/59562743/123944223-d19c4a00-d994-11eb-91ec-d4ae1e7c685f.png)


# Trimming 
![image](https://user-images.githubusercontent.com/59562743/123944093-b03b5e00-d994-11eb-9134-c635c1540973.png)

# Mapping 
![image](https://user-images.githubusercontent.com/59562743/123947189-0362e000-d998-11eb-828e-0c3ada2fea64.png)

# Quantification
![image](https://user-images.githubusercontent.com/59562743/123947308-28575300-d998-11eb-9925-a73d52693b4c.png)




