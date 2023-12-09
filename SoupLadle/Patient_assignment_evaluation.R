library(vcfR)
library(tidyverse)
library(ggpubr)
library(clue)
library(vcfR)
library(proxy)

setwd("/Users/tbleckwehl/Documents/SNP_demultiplexing/Publication/PBMC/")

VCF_ref <- read.vcfR("/Users/tbleckwehl/Documents/SNP_demultiplexing/Publication/PBMC_V3/KH30_merged_bulkRNA.vcf")
gt_ref <- extract.gt(VCF_ref, element = 'GT')

VCF_Soup <- read.vcfR("/Users/tbleckwehl/Documents/SNP_demultiplexing/Publication/PBMC_V3/cluster_genotypes.vcf")
gt_Soup <- extract.gt(VCF_Soup, element = 'GT')



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
set.seed(89)
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

  result_df <- result_df[order(result_df$chr, result_df$pos), ]
  print(fraction)
  print(nrow(result_df))
  
  #write.table(result_df, file = paste("SNP_fraction_PBMC_bulk_RNA_", fraction, ".tab", sep = ""), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  ##Added for Post
  
  shared_SNP_loci <- intersect(rownames(gt_Soup), rownames(sampled_df))
  
  # Subset both matrices based on common row names
  
  sampled_df <- filter_rows(sampled_df)
  
  SOUP <- gt_Soup[shared_SNP_loci, , drop = FALSE]
  REF <- sampled_df[shared_SNP_loci, , drop = FALSE]
  
  SOUPs <- replace_genotype_values(SOUP)
  REFs <- replace_genotype_values(REF)
  
  dist_matrix <- dissimilarity_matrix_between_matrices(t(SOUPs), t(REFs))
  
  lsap_result <- solve_LSAP(dist_matrix)
  
  for (i in 1:length(lsap_result)){
    print(colnames(dist_matrix)[lsap_result[i]])
    print(rownames(dist_matrix)[i])}
  
  
  
  assignment_df <- data.frame()
  
  for (i in 1:length(lsap_result)) {
    tmp <- data.frame("Reference" = colnames(dist_matrix)[lsap_result[i]],
                      "Soup" = rownames(dist_matrix)[i])
    
    # Create a new row in the result dataframe
    assignment_df <- rbind(assignment_df , tmp)
  }
  
  # Print the resulting dataframe
  print(assignment_df)
  assignment_df$Soup <- c(0:4)
  
  clusters <- read.table("/Users/tbleckwehl/Documents/SNP_demultiplexing/Publication/PBMC_V3/clusters.tsv", sep="\t", h=T, row.names = 1)
  
  colnames(clusters)[colnames(clusters) == "assignment"] <- "SouporCell_Assignment"
  bulk_assignments <- assignment_df$Reference[match(clusters$SouporCell_Assignment, assignment_df$Soup)]
  clusters$donor_id <- ifelse(is.na(bulk_assignments), clusters$status, bulk_assignments)
  clusters$cell <- rownames(clusters)
  
  clusters$Fraction <- fraction
  clusters$Method <- "SNP_fishing"
  clusters <- subset(clusters, clusters$donor_id!="unassigned")
  
  clusters <- calculatePrecisionRecall(df_PBMC, clusters)
  Post <- rbind(Post, clusters)
}

#Improve
# Usage:
# result_tmp <- calculatePrecisionRecall(df_PBMC, tmp)



fractions <- c(1, 0.8, 0.6, 0.4, 0.2,0.1)
Prior <- data.frame()
for (fraction in fractions) {
  tmp <- read.table(paste("./", fraction, "/vartrix_bulkRNA/donor_ids.tsv", sep=""), header=T)
  print(head(tmp))
  tmp$Fraction <- fraction
  tmp$Method <- "Vatrix - Vireo"
  tmp <- subset(tmp, tmp$donor_id!="unassigned")

  tmp <- calculatePrecisionRecall(df_PBMC, tmp)
  Prior <- rbind(Prior, tmp)
  
  tmp <- read.table(paste("./", fraction, "/vireo_bulkRNA/donor_ids.tsv", sep=""), header=T)
  print(head(tmp))
  tmp$Fraction <- fraction
  tmp$Method <- "cellSNP - Vireo"
  tmp <- subset(tmp, tmp$donor_id!="unassigned")
  
  tmp <- calculatePrecisionRecall(df_PBMC, tmp)
  
  Prior <- rbind(Prior, tmp)
}

sum_Prior <- Prior %>% group_by(Fraction, Precision, Method) %>% dplyr::count()
sum_Post <- Post %>% group_by(Fraction, Precision, Method) %>% dplyr::count()
SNP_assign <- rbind(sum_Prior,sum_Post)
SNP_assign$SNP <- nrow(gt_ref)*SNP_assign$Fraction

SNP_assign$Percentage <- SNP_assign$n/ncol(PBMC)


pdf("/Users/tbleckwehl/Documents/SNP_demultiplexing/Publication/Version4/Framework_evaluation_Prior_vs_Post.pdf", height = 6, width = 5)
ggplot(SNP_assign, aes(x=Fraction, y=n, color=Method))+geom_line(linetype = "dashed", aes(color=Method))+geom_point(aes(size=Precision))+
theme(panel.background = element_blank(),
      legend.title = element_text(size = 9, colour = "black"),
      axis.text = element_text(size = 9, colour = "black"),
      axis.title = element_text(size = 10),
      legend.position = "right")+scale_x_reverse()+
  labs(x="Fraction of bulk RNA-seq SNPs for the assignment", y="# of assigned cells")

ggplot(SNP_assign %>% filter(Fraction>0.1), aes(x=SNP, y=Percentage, color=Method))+geom_line(linetype = "dashed", aes(color=Method))+geom_point(aes(size=Precision))+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 9, colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 10),
        legend.position = "right")+scale_x_reverse()+scale_y_continuous(labels = scales::percent)+scale_size_continuous(labels = scales::percent)+
  labs(x="SNPs from bulk RNA-seq SNPs used for the assignment", y="Percentage of assigned cells")+scale_color_manual(values=c("#B2182B", "#2166AC", "#67A9CF"))


ggplot(SNP_assign %>% filter(Fraction>0.1), aes(x=SNP, y=Percentage, color=Method))+geom_line(linetype = "dashed", aes(color=Method))+geom_point(aes(size=Precision))+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 9, colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 10),
        legend.position = "right")+scale_x_reverse()+scale_y_continuous(labels = scales::percent)+scale_size_continuous(labels = scales::percent)+
  labs(x="SNPs from bulk RNA-seq SNPs used for the assignment", y="Percentage of assigned cells")+scale_color_manual(values=c("#B2182B", "#2166AC", "#67A9CF"))+
  ylim(0.6,1)

dev.off()
