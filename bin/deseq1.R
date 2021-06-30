#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
install.packages(c("tidyverse","ashr"),repos = "http://cran.us.r-project.org")
library(tidyverse)
library(ashr)
library(DESeq2)
print("Count matrix formatting")
counts=read_tsv(args[1], col_names =T, comment = "#")
#counts <-read.csv("countmatrix.tsv",sep="\t",row.names="Geneid")
cnt1 <- strsplit(colnames(counts)[-1], split= ".", fixed=TRUE)
colnames(counts)<- c("Geneid", as.vector(unlist(map(cnt1, 2))))
colnames(counts) <- gsub("/" , "", colnames(counts))
geneID <-counts[,1]
print("Read metadata")
metafile <- args[2]
metadata <- read.table(metafile, header=T)
samples<-metadata[,1]
colnames<-colnames(metadata)
metadata<-as.data.frame(metadata[,-c(1)])
colnames(metadata)<-colnames[-1]
print("Specifying design and contrast")
design = as.factor(paste0(colnames(metadata)[1:ncol(metadata)], collapse=" + "))
print("Setting up DeseqDataset from matrix")
dds <- DESeqDataSetFromMatrix(countData = counts[,-1],
                              colData = metadata,
                              design = eval(parse(text=paste0("~ ", design))))
#---------
#- Run DGE
#----------
dds <- DESeq(dds)
res <- results(dds)
res05 <- results(dds, alpha=0.05)
resOrdered <- res05[order(res05$pvalue),]
write.csv(as.data.frame(resOrdered),
          file="DEGS_results.csv")
