# Name: analyze DGE deseq CLE
# Author: EY (based off of code written by E. Dittmar)
# Date: 08/29/2023
# Version:4.1.2
# Description: Will analyze the output from the DESeq DGE with combatseq
# need Functions.R written by ED


setwd('/home/analysis')

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("Functions.R")

# now analyze 
# read in the data for 5h
DEData_pairwise_cs_5h<-ImportCSVs('deseq_results/x5h/',0.05)
# filter out significant results
mydataSig_pairwise_cs_5h<-lapply(DEData_pairwise_cs_5h,SigDEdf,PvaluesCol=7,CritP=0.05)

SigOverlapGraph_pairwise_cs_5h<-lapply(mydataSig_pairwise_cs_5h, function(x) {x$Gene})
# create an upset plot of DE expression by pairwise dev_stage
png("plots/upset_cle_5h.png", res=215, width = 2000, height=1500)
upset(fromList(SigOverlapGraph_pairwise_cs_5h),order.by="freq",nsets=13,nintersects=150, text.scale = 1.5)
dev.off()
mydataSig_pairwise_cs_5h$result_AsterCLV3_5h[1]

# see which genes overlap 
SigOverlap_pairwise_5h<-GeneSets_four(mydataSig_pairwise_cs_5h$result_AtCLV3p_5h[1], mydataSig_pairwise_cs_5h$result_AtCLV3p_S5L_5h[1],mydataSig_pairwise_cs_5h$result_AtCLV3p_5h[1],mydataSig_pairwise_cs_5h$result_AtCLV3p_S5L_5h[1])
names(SigOverlap_pairwise_5h)
lapply(SigOverlap_pairwise_5h,function(x) {length(x$Gene)})


# read in the data for 24h
DEData_pairwise_cs_24h<-ImportCSVs('deseq_results/x24h/',0.05)
# filter out significant results
mydataSig_pairwise_cs_24h<-lapply(DEData_pairwise_cs_24h,SigDEdf,PvaluesCol=7,CritP=0.05)

SigOverlapGraph_pairwise_cs_24h<-lapply(mydataSig_pairwise_cs_24h, function(x) {x$Gene})

# create an upset plot of DE expression by pairwise dev_stage
png("plots/upset_cle_24h.png", res=215, width = 2000, height=1500)
upset(fromList(SigOverlapGraph_pairwise_cs_24h),order.by="freq",nsets=13,nintersects=150, text.scale = 1.5)
dev.off()


# see which genes overlap 
SigOverlap_pairwise_24h<-GeneSets_four(mydataSig_pairwise_cs_24h$result_AtCLV3p_24h[1], mydataSig_pairwise_cs_24h$result_AtCLV3p_S5L_24h[1],mydataSig_pairwise_cs_24h$result_AtCLV3p_24h[1],mydataSig_pairwise_cs_24h$result_AtCLV3p_S5L_24h[1])
names(SigOverlap_pairwise_24h)
lapply(SigOverlap_pairwise_24h,function(x) {length(x$Gene)})


# see where WUS is
# transcript 1, 5h
DEData_pairwise_cs_5h$result_AtCLV3p_5h%>% filter_all(any_vars(. %in% c("g51546.t1")))
DEData_pairwise_cs_5h$result_AtCLV3p_S5L_5h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 
DEData_pairwise_cs_5h$result_AsterCLV3_5h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 
DEData_pairwise_cs_5h$result_AsterCLV3_L5S_5h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 


# transcript 2, 5h
DEData_pairwise_cs_5h$result_AtCLV3p_5h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 
DEData_pairwise_cs_5h$result_AtCLV3p_S5L_5h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 
DEData_pairwise_cs_5h$result_AsterCLV3_5h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 
DEData_pairwise_cs_5h$result_AsterCLV3_L5S_5h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 


# transcript 1, 24h
DEData_pairwise_cs_24h$result_AtCLV3p_24h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 
DEData_pairwise_cs_24h$result_AtCLV3p_S5L_24h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 
DEData_pairwise_cs_24h$result_AsterCLV3_24h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 
DEData_pairwise_cs_24h$result_AsterCLV3_L5S_24h%>% filter_all(any_vars(. %in% c("g51546.t1"))) 


# transcript 2, 24h
DEData_pairwise_cs_24h$result_AtCLV3p_24h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 
DEData_pairwise_cs_24h$result_AtCLV3p_S5L_24h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 
DEData_pairwise_cs_24h$result_AsterCLV3_24h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 
DEData_pairwise_cs_24h$result_AsterCLV3_L5S_24h%>% filter_all(any_vars(. %in% c("g62225.t1"))) 



