# Name: load GC data
# Author: EY
# Date: 12/9/2022 Updated: 03/06/2023
# Version:4.1.2
# Description: will load in the parsed gene count data
# This is done using DESeq

# set cwd
setwd('sunflower_dev_data/gene_count_data/')
library("DESeq2")

# put all of the gene_count files into a character vector in order to load them 
gene_count_files<-dir(pattern="*\\.tab$")

samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

# read the data into a dataframe
data_table<-data.frame(sampleName=samples,fileName=gene_count_files)


# load the DESeq data set
# we read in the data using the DESeqDataSetFromHTSeqCount bc our dataset is in the same format as HTSeq data
dds_set<-DESeqDataSetFromHTSeqCount(sampleTable=data_table,directory='sunflower_dev_data/gene_count_data/',design=~0)

# save the data set 
saveRDS(dds_set,file="gene_count_sunflower_dev_deseq.Rdata")


