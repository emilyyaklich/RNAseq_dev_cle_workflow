# Name: run_combatseq.R
# Author: EY
# Date: 03/06/2023
# Version:4.2.1
# Description: Will run combatseq normalization

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
source("sunflower/Functions.R")


# read in and process data

setwd('/home/dev_RNAseq/')

# read in the data matrix
summed_counts<-readRDS("/home/dev_RNAseq/sunflower/gene_count_sunflower_dev_deseq.Rdata")
samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)

# create the model (wrt dev_stage)
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)


# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

extracted_row <- subset(summed_counts, rownames(summed_counts) == 'g24641.t2')
counts(extracted_row)


# get a dataframe of the counts
count_matrix <- as.matrix(counts(summed_counts_filt))

# sort by batch (sample group) and group (dev_stage)
batch <- c(1,1,1,1,2,3,2,3,2,3,2,3)
group <- c(1,2,3,4,1,1,2,2,3,3,4,4)

# adjust counds using combat seq...output is a matrix of adjusted counts 
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)

# write to a CSV file
write.csv(as.data.frame((adjusted_counts)), file='sunflower/adjusted_counts_combatseq.csv')

# plot MDS of adjusted counts...can compare this with previous plot pre-combat seq
glimmaMDS(adjusted_counts, group=metadata)
