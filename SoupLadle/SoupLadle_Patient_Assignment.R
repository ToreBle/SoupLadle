library(ggplot2)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(clue)
library(vcfR)
library(proxy)
library(viridis)
library(RColorBrewer)


assign_clusters_to_samples <- function(vcf_Soup, vcf_bulk) {

  # Filter SNPs that are not common in all samples
  filter_rows <- function(dataframe) {
    # Check for rows where more than one non-NA value exists
    discard_rows <- apply(dataframe, 1, function(row) {
      # Filter out NA values
      non_na_values <- na.omit(row)
      # Check if more than one non-NA value exists
      length(unique(non_na_values)) <= 1
    })
    
    # Subset the dataframe to keep rows that don't meet the discard condition
    filtered_df <- dataframe[!discard_rows | rowSums(!is.na(dataframe)) <= 1, ]
    
    return(filtered_df)
  }
  
  #Replace the variant with a binary string
  replace_genotype_values <- function(matrix) {
    matrix[matrix == "0/0"] <- "100"
    matrix[matrix == "0/1"] <- "010"
    matrix[matrix == "1/1"] <- "001"
    matrix[matrix == "0/2"] <- "011"
    matrix[matrix == "0/3"] <- "011"
    matrix[is.na(matrix)] <- "000"
    
    return(matrix)
  }
  
  # Calculate a dissimilarity matrix for Hamming distances between the SNPs of the reference
  # and the deconvoluted single-cell clusters
  dissimilarity_matrix_between_matrices <- function(matrix1, matrix2) {
    if (nrow(matrix1) != nrow(matrix2)) {
      stop("Matrices must have the same number of rows")
    }
    
    num_rows <- nrow(matrix1)
    dissimilarity_matrix <- matrix(0, nrow = num_rows, ncol = num_rows)
    
    # Get row names from the original matrices
    rownames_matrix1 <- rownames(matrix1)
    rownames_matrix2 <- rownames(matrix2)
    
    for (i in 1:num_rows) {
      row1 <- matrix1[i, , drop = FALSE]
      
      for (j in 1:num_rows) {
        row2 <- matrix2[j, , drop = FALSE]
        
        # Ensure both rows have the same length
        if (length(row1) != length(row2)) {
          stop("Rows must have the same length")
        }
        
        # Calculate Hamming distance for the current pair of rows
        dissimilarity_matrix[i, j] <- sum(row1 != row2)
      }
    }
    # Set row and column names for the dissimilarity matrix
    rownames(dissimilarity_matrix) <- rownames_matrix1
    colnames(dissimilarity_matrix) <- rownames_matrix2
    
    return(dissimilarity_matrix)
  }
  
  
  read_extract_filter_vcf <- function(vcf_file) {
    cat(paste("Processing ", vcf_file, ":\n", sep=""))
    VCF <- read.vcfR(vcf_file)
    gt <- extract.gt(VCF, element = 'GT')
    filtered_gt <- filter_rows(gt)
    cat(paste("... Completed\n", sep=""))
    return(filtered_gt)
  }
  
  gt_Soup <- read_extract_filter_vcf(vcf_Soup)
  gt_ref <- read_extract_filter_vcf(vcf_bulk)
  
  # Subset both matrices based on common row names
  shared_SNP_loci <- intersect(rownames(gt_Soup), rownames(gt_ref))
  SOUP <- gt_Soup[shared_SNP_loci, , drop = FALSE]
  REF <- gt_ref[shared_SNP_loci, , drop = FALSE]
  
  # Binary string
  SOUPs <- replace_genotype_values(SOUP)
  REFs <- replace_genotype_values(REF)
  
  # Dissimilarity matrix
  dist_matrix <- dissimilarity_matrix_between_matrices(t(SOUPs), t(REFs))
  
  # Hungarian algorithm
  lsap_result <- solve_LSAP(dist_matrix)
  
  assignment_df <- data.frame(
    "Reference" = colnames(dist_matrix)[lsap_result],
    "Cluster" = rownames(dist_matrix))
  
  return(assignment_df)

}

# Usage
assignments <- assign_clusters_to_samples("cluster_genotypes.vcf", "KH30_merged_bulkRNA.vcf")


#Add meta data column
clusters <- read.table("clusters.tsv", sep="\t", h=T, row.names = 1)
colnames(clusters)[colnames(clusters) == "assignment"] <- "SouporCell_Assignment"
bulk_assignments <- assignments$Reference[match(clusters$SouporCell_Assignment, assignments$Cluster)]
clusters$SoupLadle <- ifelse(is.na(bulk_assignments), clusters$status, bulk_assignments)
#Export as meta_data for Seurat / Scanpy and others

#Read PBMC data set for evaluation
PBMC <- readRDS("KH30_PBMC_Multiplex.rds")


pdf("SNP_fishing_Posterior_Assignment.pdf", height = 6, width = 6)
Patient_colors <- c("#195E83", "#DFC8A2", "#E07B39", "#69BDD2", "#80391E", "#606060", "gray")
clusters <- clusters %>%
  mutate(SoupLadle  = recode(SoupLadle,
                               "KH_75" = "bulk_SNP_Patient_1",
                               "KH_77" = "bulk_SNP_Patient_2",
                               "KH_79" = "bulk_SNP_Patient_3",
                               "KH_76" = "bulk_SNP_Patient_4",
                               "KH_78" = "bulk_SNP_Patient_5"),
  )
PBMC <- AddMetaData(PBMC, metadata = clusters)


DimPlot(PBMC, group.by = c("SoupLadle"))+ coord_equal()+scale_color_manual(values = Patient_colors)
VlnPlot(PBMC, features = c("log_prob_singleton"), group.by = "SoupLadle")+scale_fill_manual(values = Patient_colors)
VlnPlot(PBMC, features = c("log_prob_singleton"), group.by = "SoupLadle", pt.size = "")+scale_fill_manual(values = Patient_colors)
FeaturePlot(PBMC, features = c("doublet_posterior"))+ coord_equal()
dev.off()

#The Doublet score can be used for filtering before any downstream analysis. This are SNP-proven doublets.
#However, low quality nuclei/cells and intera-patient doublets might persist and can be filtered out in the normal QC steps.