# Name: plot_MDS
# Author: EY
# Date: 12/24/2022 (edited 08/21/2023)
# Version:4.1.2
# Description: will read the data matrix from load_GC_data_and_sum_reps into
# DESeq2, read in metadata,pre-filter and plotMDS

library(readxl)
library(DESeq2)
library(dplyr)
library(Glimma)
library(sva)
install.packages("Glimma")

setwd('/home/dev_RNAseq/')

# read in the data matrix
summed_counts<-readRDS("/home/dev_RNAseq/sunflower/gene_count_sunflower_dev_deseq.Rdata")
dim(summed_counts)


summed_counts


samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
           "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
           "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)


metadata<-data.frame(samples, dev_stage)
write.csv(as.data.frame(metadata), file='sunflower/metadata.csv')

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)



# create the model
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)
summed_counts$samples
# pre-filter for reads where at least 3 samples (in ANY STAGE) have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

counts(summed_counts)>=1

# will load the plot...need to save within the html
glimmaMDS(summed_counts_filt, groups=metadata)

# plot PCA 
#png("plots/pca_raw_data.png", res=215, width = 1200, height=1000)
#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#DESeq2::plotPCA((summed_counts_filt),labels=FALSE,col=dev_stage)
#legend("topright",inset=c(-0.4,0),legend=unique(dev_stage), fill=dev_stage)
#dev.off()

dds_set
