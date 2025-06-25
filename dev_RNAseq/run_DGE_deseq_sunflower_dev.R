# Name: run_DGE_sunflower_dev.R
# Author: EY
# Date: 08/29/2023
# Version:4.1.2
# Description: Will run the differential gene expression on the infloresence stages accounting for variation with combatseq

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)

# read in and process data
setwd('/home/dev_RNAseq/')

source("Functions.R")

packageVersion("DESeq2")



adjusted_counts<-read.csv('sunflower/adjusted_counts_combatseq.csv', row.names=1)
metadata<-read.csv('sunflower/metadata.csv', row.names=1)
# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)


# create the model with adjusted counts
adjusted_counts_deseq<-DESeqDataSetFromMatrix((adjusted_counts),colData = metadata, design=~0+dev_stage)

# run the DGE analysis
DESeq_dataset_results_combatseq<-DESeq(adjusted_counts_deseq,parallel=TRUE)



normalized_counts<-counts(DESeq_dataset_results_combatseq, normalized=TRUE)
write.csv(normalized_counts, "sunflower/deseq_results/normalized_counts_size_factor.csv")




resultsNames(DESeq_dataset_results_combatseq)

# save this because it will be read in to WGCNA (WGCNA takes all genes into account, but uses the normalization of DESeq)
saveRDS(DESeq_dataset_results_combatseq, file='sunflower/deseq_results/deseq_dataset_results_pairwise_combatseq.RData')

# set up the pairwise contrasts (10v20, 20v30, 30v35)
result_10D_v_20D_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","20D","10D"),alpha=0.05,parallel=TRUE)
result_20D_v_30D_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","30D","20D"),alpha=0.05,parallel=TRUE)
result_30D_v_35D_combatseq<-results(DESeq_dataset_results_combatseq,contrast=c("dev_stage","35D","30D"),alpha=0.05,parallel=TRUE)

# write to CSV file
write.csv(as.data.frame(result_10D_v_20D_combatseq), file='sunflower/deseq_results/pairwise/result_10D_v_20D.csv')
write.csv(as.data.frame(result_20D_v_30D_combatseq), file='sunflower/deseq_results/pairwise/result_20D_v_30D.csv')
write.csv(as.data.frame(result_30D_v_35D_combatseq), file='sunflower/deseq_results/pairwise/result_30D_v_35D.csv')


# do the individual contrasts
# look at results for each treatment
results_10D<-results(DESeq_dataset_results_combatseq, alpha=0.05,name='dev_stage10D')
results_20D<-results(DESeq_dataset_results_combatseq, alpha=0.05,name='dev_stage20D')
results_30D<-results(DESeq_dataset_results_combatseq, alpha=0.05,name='dev_stage30D')
results_35D<-results(DESeq_dataset_results_combatseq, alpha=0.05,name='dev_stage35D')

write.csv(as.data.frame(results_10D), file='sunflower/deseq_results/indiv_res/result_10D.csv')
write.csv(as.data.frame(results_20D), file='sunflower/deseq_results/indiv_res/result_20D.csv')
write.csv(as.data.frame(results_30D), file='sunflower/deseq_results/indiv_res/result_30D.csv')
write.csv(as.data.frame(results_35D), file='sunflower/deseq_results/indiv_res/result_35D.csv')
