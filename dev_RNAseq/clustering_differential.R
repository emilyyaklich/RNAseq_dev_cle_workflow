# Name: clustering differential
# Author: EY 
# Date: 09/03/2024 
# Version:4.1.2
# Description: Will cluster the DE output


##### setup ####

setwd('/home/dev_RNAseq/')

# load packages
library(tidyverse)
library(tidyr)
library(DESeq2)
source("sunflower/Functions.R")
library(DEGreport)
library(factoextra)
library(cluster)
library(fpc)


# read in metadata
metadata<-read.csv('sunflower/metadata.csv', row.names=1)

# read in deseq results (output from run_DGE_deseq_sunflower_dev.R)
# note this is reading in all of the data and is NOT filtered for genes that are
# differentially expressed 
deseq <- readRDS('sunflower/deseq_results/deseq_dataset_results_pairwise_combatseq.RData')



# transform data into a matrix
vsd <- vst(deseq,blind=FALSE)
vsd_matrix<-assay(vsd)
write.csv(vsd_matrix, "sunflower/deseq_results/normalized_counts_vst.csv")

# now analyze 
# read in the data
DEData_pairwise_cs<-ImportCSVs('sunflower/deseq_results/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)



# Initialize an empty vector to store gene names
DE_genes <- c()

# Loop through each dataframe in the list
for (df in mydataSig_pairwise_cs) {
  # Filter rows where padj < 0.05 and extract the gene names
  filtered_genes <- df$Gene[df$padj < 0.05]
  # Append the gene names to the vector
  DE_genes <- c(DE_genes, filtered_genes)
}

# Get unique gene names
DE_genes <- unique(DE_genes)



subset_matrix <- vsd_matrix[rownames(vsd_matrix) %in% DE_genes, ]



average_expression_matrix <- cbind(
  D10 = rowMeans(subset_matrix[, c(1, 5, 6)]), 
  D20 = rowMeans(subset_matrix[, c(2, 7, 8)]), 
  D30 = rowMeans(subset_matrix[, c(3, 9, 10)]), 
  D35 = rowMeans(subset_matrix[, c(4, 11, 12)]))

scaled_expression_matrix <- t(scale(t(average_expression_matrix)))


gene_dist <- dist(scaled_expression_matrix, method="euclidean")

gene_hclust <- hclust(gene_dist, method = "complete")



png("sunflower/plots/clustering_differential/hclust_tree_all_cuts.png", width=2700, height=2100, res=300)
plot(gene_hclust, labels = FALSE)
abline(h = 3.46, col = "red", lwd = 2)
abline(h = 3.3, col = "purple", lwd = 2)
abline(h = 3.0, col = "brown", lwd = 2)
abline(h = 2.7, col = "blue", lwd = 2)# add horizontal line to illustrate cutting dendrogram
abline(h = 2.4, col = "salmon", lwd = 2)
abline(h = 2.0, col = "darkgreen", lwd = 2)
abline(h = 1.5, col = "darkorange", lwd = 2)
abline(h = 1.0, col = "yellow", lwd = 2)
abline(h = 0.1, col = "magenta", lwd = 2)
dev.off()



# plot tree heights
png("sunflower/plots/clustering_differential/hclust_treeheights.png", width=2700, height=2100, res=300)
heights<-sort(gene_hclust$height, decreasing = TRUE)
plot(heights, type="p", xlab = "Number of Clusters", ylab= "Tree Cut Position")
dev.off()



# plot clusters for differing cut heights

cut_heights <- c(3.46, 3.3, 3.0, 2.7, 2.4,2.0, 1.5, 1.0, 0.1)

average_line_colors <- c(
  "3.46" = "red",       
  "3.3"  = "purple",    
  "3.0"  = "blue",     
  "2.7"  = "brown",      
  "2.4"  = "salmon",    
  "2.0"  = "darkgreen", 
  "1.5"  = "darkorange", 
  "1.0" = "yellow",
  "0.1" = "magenta"
)
output_dir <- "sunflower/plots/clustering_differential/"

# loop through each height to cut the tree and plot
for (i in seq_along(cut_heights)) {
  h <- cut_heights[i]
  
  # cut the dendrogram at the specified height
  gene_cluster <- cutree(gene_hclust, h = h)
  gene_cluster_df <- enframe(gene_cluster)
  
  # rename columns
  names(gene_cluster_df) <- c("Gene", "cluster")
  
  
  # check which clusters contain specific genes (WUS and CLV3)
  specific_genes <- c("g51546.t1", "g23024.t1")  # WUS and CLV3
  clusters_with_genes <- gene_cluster_df %>%
    filter(Gene %in% specific_genes) %>%
    select(Gene, cluster)
  
  # print the clusters for specific genes
  wus_cluster <- clusters_with_genes %>% filter(Gene == "g51546.t1") %>% pull(cluster)
  clv_cluster <- clusters_with_genes %>% filter(Gene == "g23024.t1") %>% pull(cluster)
  print(paste("Cut Height:", h))
  print(paste("WUS: Cluster", wus_cluster))
  print(paste("CLV3: Cluster", clv_cluster))
  

  # prepare average expression data
  average_expression_df <- as.data.frame(scaled_expression_matrix)
  average_expression_df$Gene <- rownames(scaled_expression_matrix)
  

  df_cluster <- average_expression_df %>% 
    inner_join(gene_cluster_df, by = "Gene")
  
  # write clustering data to CSV
  write.csv(df_cluster, file = paste0("sunflower/deseq_results/clustering/", "clustering_data_", h, ".csv"), row.names = FALSE)
  
  # reshape data for plotting
  df_long <- df_cluster %>%
    pivot_longer(cols = starts_with(c("D")), names_to = "samples", values_to = "Expression")
  
  cluster_gene_count <- df_long %>%
    group_by(cluster) %>%
    summarise(number_of_genes = n_distinct(Gene))
  

  facet_labels <- cluster_gene_count %>%
    mutate(label = paste("Cluster", cluster, ":", number_of_genes)) %>%
    pull(label, name = cluster)
  

  facet_labels <- setNames(facet_labels, cluster_gene_count$cluster)
  
  

  df_long$samples <- as.numeric(gsub("D", "", df_long$samples))
  
  

  plot_file <- paste0(output_dir, "hclust_clusters_", h, "_cut.png")
  

  p <- ggplot(df_long, aes(x = samples, y = Expression, group = Gene)) +
    geom_line() +  # Set color for all lines
    geom_line(stat = "summary", fun = "mean", size = 1.5, color = average_line_colors[as.character(format(h, nsmall=1))], aes(group = 1)) +  # Average line
    facet_wrap(~ cluster, labeller = labeller(cluster = facet_labels)) +
    labs(
      title = paste("Clustered Expression Data (Cut Height =", h, ")"),
      x = "Developmental Stage",
      y = "Scaled Expression"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1) + scale_y_continuous(limits = c(min(df_long$Expression, na.rm = TRUE), max(df_long$Expression, na.rm = TRUE)))
)
  

  ggsave(plot_file, plot = p, width = 2700/300, height = 2100/300, dpi = 300)
}


# just show two clusters CLV3 and WUS is in for the lowest cut

clusters_of_interest <- c("558", "126")  

# subset the data to include only the clusters of interest
df_subset <- df_long %>%
  filter(cluster %in% clusters_of_interest)

# update the cluster_gene_count and facet_labels for the selected clusters
cluster_gene_count_subset <- df_subset %>%
  group_by(cluster) %>%
  summarise(number_of_genes = n_distinct(Gene))

facet_labels_subset <- cluster_gene_count_subset %>%
  mutate(label = paste("Cluster", cluster, ":", number_of_genes)) %>%
  pull(label, name = cluster)

facet_labels_subset <- setNames(facet_labels_subset, cluster_gene_count_subset$cluster)


df_subset$samples <- as.numeric(gsub("D", "", df_subset$samples))


png("sunflower/plots/clustering_differential/hclust_clusters_0_1_cut_subset.png", width=2700, height=2100, res=300)
ggplot(df_subset, aes(x = samples, y = Expression, group = Gene)) +
  geom_line() +
  geom_line(stat = "summary", fun = "mean", color = "darkorange", size = 1.5, aes(group = 1)) +
  facet_wrap(~ cluster, labeller = labeller(cluster = facet_labels_subset)) +
  scale_x_continuous(
    breaks = c(10, 20, 30, 35),
    labels = c("D10", "D20", "D30", "D35")
  ) +
  labs(
    x = "Developmental Stage",
    y = "Scaled Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



