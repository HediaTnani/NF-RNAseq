#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
install.packages(c("tidyverse","ashr"),repos = "http://cran.us.r-project.org")
library(tidyverse)
library(ashr)
library(DESeq2)
cat("Count matrix formatting")
counts=read_tsv("countmatrix.tsv", col_names =T, comment = "#")
#counts <-read.csv("countmatrix.tsv",sep="\t",row.names="Geneid")
cnt1 <- strsplit(colnames(counts)[-1], split= ".", fixed=TRUE)
colnames(counts)<- c("Geneid", as.vector(unlist(map(cnt1, 2))))
colnames(counts) <- gsub("/" , "", colnames(counts))
geneID <-counts[,1]
cat("Read metadata")
metafile <- args[1]
metadata <- read.table(metafile, header=T)
samples<-metadata[,1]
colnames<-colnames(metadata)
metadata<-as.data.frame(metadata[,-c(1)])
colnames(metadata)<-colnames[-1]
cat(" Specifying design and contrast")
if(ncol(metadata) > 1){
  design = as.factor(paste0(colnames(metadata)[1:ncol(metadata)], collapse=" + "))
  if(!(condition %in% colnames(metadata))) {
    stop("ERROR: condition '",levels(condition),"' is not a column of the metadata file, please specify a valid condition")
}
} else {
  design = as.factor(paste0(colnames(metadata)[1], collapse=" + "))
}

cat("Setting up DeseqDataset from matrix")
dds <- DESeqDataSetFromMatrix(countData = counts[,-1],
                              colData = metadata,
                              design = eval(parse(text=paste0("~ ", design))))
#---------
#- Run DGE
#----------
dds <- DESeq(dds)

for (i in 1:ncol(metadata)) {
  pairs<-combn(unique(metadata[,i]),2)
  for (j in 1:ncol(pairs)) {

    #- alpha is fdr threshold for summary display only
    res<-results(dds, contrast = c(colnames(metadata)[i],as.character(pairs[,j])), alpha=0.05)
    resNorm <- lfcShrink(dds, contrast = c(colnames(metadata)[i],as.character(pairs[,j])), res=res, type="ashr")
    write.table(resNorm, file=paste0(colnames(metadata)[i], "_", paste0(as.character(pairs[,j]),collapse="_vs_"), "_results.txt"), sep="\t",quote=F,row.names=T, col.names=NA)
  }
}

