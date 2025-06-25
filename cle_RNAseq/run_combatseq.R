# Name: run_combatseq.R
# Author: EY
# Date: 03/06/2023
# Version:4.1.2
# Description: Will run combatseq normalization


setwd('/home/cle_RNAseq/analysis')


library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
source("Functions.R")

# read in and process data

# read in the data matrix
summed_counts<-readRDS("/home/cle_RNAseq/analysis/gene_count_sunflower_cle_deseq.Rdata")
samples<-colnames(summed_counts)
samples
cle_id <- paste0("CLE", gsub(".*-(\\d+)_.*", "\\1", samples))
cle_id <- as.factor(cle_id)
metadata<-data.frame(samples, cle_id)

# create the factors of interest
metadata$cle_id<-factor(metadata$cle_id)
write.csv(as.data.frame((metadata)), file='metadata.csv')

# create the model (wrt dev_stage)
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+cle_id)

# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# get a dataframe of the counts
count_matrix <- as.matrix(counts(summed_counts_filt))
write.csv(as.data.frame((count_matrix)), file='raw_cle_counts.csv')


reordered <- as.numeric(gsub("CLE", "", metadata$cle_id))

# Sort the dataframe based on numeric part
metadata <- metadata[order(reordered), ]



# sort by batch (sample group) and group (CLE-ID)
batch <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
group <- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10)

# adjust counds using combat seq...output is a matrix of adjusted counts 
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)

# write to a CSV file
write.csv(as.data.frame((adjusted_counts)), file='adjusted_counts_combatseq.csv')

# plot MDS of adjusted counts...can compare this with previous plot pre-combat seq
glimmaMDS(adjusted_counts, group=metadata)
glimmaMDS(count_matrix, group=metadata)





