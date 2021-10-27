#!/usr/bin/env Rscript

if (!require("purrr")){
  install.packages("purrr", dep=T)
  suppressPackageStartupMessages(library("purrr"))
}-

if (!require("tidyverse")){
  install.packages("tidyverse", dep=T)
  suppressPackageStartupMessages(library("tidyverse"))
}


library(purrr)
library(tidyverse)


f_files<- list.files(pattern = "counts.tsv", full.names = F)
read_in_feature_counts<- function(file){
 cnt<- read_tsv(file, col_names =T, comment = "#")
 return(cnt) 
}   
    
raw_counts<- purrr::map(f_files, read_in_feature_counts)
raw_counts_df<- purrr::reduce(raw_counts, inner_join) 
write_tsv(raw_counts_df, "countmatrix.tsv")

