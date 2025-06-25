# Name: load GC data
# Author: EY
# Date: 12/9/2022 Updated: 03/06/2023
# Version:4.1.2
# Description: will load in the parsed gene count data
# This is done using DESeq

# set cwd
setwd('/scratch/sunflower_cle_data/gene_count_data/')
library("DESeq2")

# put all of the gene_count files into a character vector in order to load them 
gene_count_files<-dir(pattern="*\\.tab$")

# shorted names 
shortened_names <- gsub("_S\\d+_L004ReadsPerGene.out.tab", "", gene_count_files)


# read the data into a dataframe
data_table<-data.frame(sampleName=shortened_names,fileName=gene_count_files)


# load the DESeq data set
# we read in the data using the DESeqDataSetFromHTSeqCount bc our dataset is in the same format as HTSeq data
dds_set<-DESeqDataSetFromHTSeqCount(sampleTable=data_table,directory='sunflower_cle_data/gene_count_data/',design=~0)

# save the data set 
saveRDS(dds_set,file="/home/cle_RNAseq/analysis/gene_count_sunflower_cle_deseq.Rdata")


