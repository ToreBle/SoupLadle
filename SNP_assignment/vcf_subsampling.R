library(tidyverse)
library(ggpubr)
library(clue)
library(vcfR)
library(proxy)

setwd("/data/SNP_Multiplexing/output/KH30_unfiltered/")

#Read-in the VCF of the merged SNPs from all bulkRNA-seq samples
VCF_ref <- read.vcfR("KH30_merged_bulkRNA.vcf")
gt_ref <- extract.gt(VCF_ref, element = 'GT')

#VCF_Soup <- read.vcfR("./output_souporcell_SNP_only/cluster_genotypes.vcf")
#gt_Soup <- extract.gt(VCF_Soup, element = 'GT')

calculatePrecisionRecall <- function(df_PBMC, tmp) {
  # Inner join the dataframes
  df <- inner_join(df_PBMC, tmp, by = "cell")
  df <- df %>% dplyr::select(donor_id, Clustering)
  
  # Define true labels and set levels
  true_labels <- df$Clustering
  true_labels <- factor(true_labels, levels = c("KH_75", "KH_76", "KH_77", "KH_78", "KH_79", "doublet"))
  
  # Create a copy of the dataframe for predictions
  predicted_df <- df
  predicted_df$Clustering <- NULL
  
  # Initialize variables to store precision and recall for each model
  models <- colnames(predicted_df)
  precisions <- vector("numeric", length(models))
  recalls <- vector("numeric", length(models))
  
  # Calculate precision and recall for each model
  for (i in 1:length(models)) {
    model_predictions <- predicted_df[, i]
    TP <- sum(true_labels == model_predictions)
    FP <- sum(model_predictions %in% unique(true_labels)) - TP
    FN <- sum(true_labels %in% unique(model_predictions)) - TP
    
    precisions[i] <- TP / (TP + FP)
    recalls[i] <- TP / (TP + FN)
  }
  
  # Create a summary data frame
  results_df <- data.frame(Method = models, Precision = precisions, Recall = recalls)
  
  # Add precision and recall to tmp dataframe
  tmp$Precision <- results_df$Precision
  tmp$Recall <- results_df$Recall
  
  return(tmp)
}

#Prepare meta data
PBMC <- readRDS("KH30_PBMC_Multiplex.rds")
df_PBMC <- PBMC@meta.data
df_PBMC$cell <- rownames(df_PBMC)
df_PBMC <- df_PBMC %>% mutate(Clustering = recode(Cluster,
                                                  "CD14+ Monocytes" = "KH_75",
                                                  "CD16+ Granulo" = "KH_77",
                                                  "CD4+ T cells" = "KH_79",
                                                  "CD8+ T cells" = "KH_76",
                                                  "CD19+ B cells" = "KH_78",
                                                  "Mixed" = "doublet",
                                                  "Unassigned" = "Unassigned"))
df_PBMC <- df_PBMC %>% select(cell, Clustering)

#Initialize the seeds
for (seed in seq(10, 100, by = 10)) {
  set.seed(seed)
  print(seed)
  
  #Define the fraction of SNPs to be considered
  fractions <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1)
  result_df <- data.frame()
  Post <- data.frame()
  
  for (fraction in fractions) {
    sampled_df <- data.frame(gt_ref) %>%
      sample_frac(fraction, replace = FALSE)
    # Extract "chr" and "pos" values for all row names
    result_df <- data.frame(chr = character(0), pos = numeric(0))  # Initialize an empty result dataframe
    
    # Loop through row names and extract values
    for (row_name in rownames(sampled_df)) {
      # Split the row name into "chr" and "pos"
      split_string <- strsplit(row_name, "_")
      chr_value <- split_string[[1]][1]
      pos_value <- as.numeric(split_string[[1]][2])
      
      # Append to the result dataframe
      result_df <- rbind(result_df, data.frame(chr = chr_value, pos = pos_value))
    }
    
    #order and check the results
    result_df <- result_df[order(result_df$chr, result_df$pos), ]
    print(fraction)
    print(nrow(result_df))
    
    write.table(result_df, file = paste("./SNP_assignment/SNP_fraction_PBMC_bulk_RNA_", fraction, "_", seed, ".tab", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  }}