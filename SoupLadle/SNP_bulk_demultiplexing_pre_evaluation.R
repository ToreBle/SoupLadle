library(GenomicFeatures)
library(RCAS)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)

  args <- commandArgs(trailingOnly = TRUE)
#Adapt directly from vcf file to read

  transcript_features <- data.frame()
  sites <- read.table(paste(args[2], "/SNP_evaluation/bulk_RNA_isec/sites.txt", sep=""))
  print(head(sites))
  sites$end = sites$V2+1
  sites <- sites[,c(1,2,6,3,4,5)]
  colnames(sites) <- c("chromosome", "start", "end", "REF", "ALT", "Sample")
  
  #Translate the binary Sample information to a column
  sites$Binary <- as.character(sites$Sample)
  
  # Determine the maximum length of binary numbers
  max_length <- max(nchar(sites$Binary))
  
  # Pad the binary numbers with leading zeros to ensure equal length
  sites$Binary <- sprintf("%0*s", max_length, sites$Binary)
  
  # Create new columns for each binary digit position
  for (i in 1:max_length) {
    col_name <- paste0("Patient_", i)
    sites[[col_name]] <- as.integer(sapply(sites$Binary, function(x) substr(x, i, i) == "1"))
  }
  
  #Determine the SNPs per Patient
  sites$SNP_type <- sapply(sites$Sample, function(x) sum(as.numeric(strsplit(as.character(x), "")[[1]])))
  
  #sites <- sites %>%
  #  rename(
  #    Patient_1 = Patient_1,
  #    Patient_3 = Patient_2,
  #    Patient_5 = Patient_3,
  #    Patient_2 = Patient_4,
  #    Patient_4 = Patient_5)
  
  #Pivot the dataframe by the Patients
  df_long <- pivot_longer(sites, cols = starts_with("Patient_"), names_to = "Samples", values_to = "Value") %>%
    filter(Value == "1") %>% group_by(SNP_type, Samples) %>% dplyr::count()

  df_Patient <- sites %>% dplyr::select(Patient_1, Patient_2, Patient_3, Patient_4, Patient_5)
  df_char <- data.frame(lapply(df_Patient, as.character), stringsAsFactors = FALSE)
  
  
  calculate_co_occurrence_percentage <- function(patient1, patient2) {
    # Check if both vectors have the same length
    if (length(patient1) != length(patient2)) {
      stop("Vectors 'patient1' and 'patient2' must have the same length.")
    }
    
    # Calculate the total number of elements
    total_elements <- length(patient1)
    
    # Calculate the intersection, i.e., the number of elements where both are equal to 1
    intersection <- sum(patient1 == 1 & patient2 == 1)
    
    # Calculate the percentage of elements where both are equal to 1
    percentage <- (intersection / total_elements) * 100
    
    return(percentage)
  }
  
  
  # Initialize a matrix to store co-occurrence percentages
  num_columns <- ncol(df_Patient)
  co_occurrence_matrix <- matrix(0, nrow = num_columns, ncol = num_columns, dimnames = list(names(df_Patient), names(df_Patient)))
  
  # Calculate co-occurrence percentages for each column pair
  for (i in 1:num_columns) {
    for (j in 1:num_columns) {
      co_occurrence_matrix[i, j] <- calculate_co_occurrence_percentage(df_char[, i], df_char[, j])
    }
  }
  
  co_occurrence_df <- as.data.frame(as.table(co_occurrence_matrix))
  
  
pdf(paste(args[2], "/SNP_evaluation/SNP_fishing_bulk_SNP_evaluation_", args[3],".pdf", sep=""), height = 8, width = 12)
  

  colors <- brewer.pal(7, name = "RdYlBu")
  m = make_comb_mat(df_Patient)
  UpSet(m, pt_size = unit(5, "mm"), lwd = 3, set_order = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5"), comb_order = order(comb_size(m)), 
        comb_col = colors[comb_degree(m)])

  p1 <- ggplot(co_occurrence_df, aes(Var1, Var2, fill = Freq)) +
    geom_tile() + scale_fill_distiller(palette = "RdBu", direction = -1, labels = scales::percent_format(scale = 1)) +
    labs(title = "Co-occurrent SNPs", x = "", y = "", fill="Shared\nSNPs") +
    theme_minimal() + coord_equal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 12),
          axis.text = element_text(size = 10, colour = "black"),
          legend.position = "top") +
    geom_text(aes(label = sprintf("%.1f%%", Freq)), vjust = 0.5, size = 2)
  
  p2 <- ggplot(df_long, aes(y = Samples, x = n, fill = factor(SNP_type))) +
    geom_bar(position="stack", stat = "identity") +
    theme(panel.background = element_blank(),
          legend.title = element_text(size = 9, colour = "black"),
          axis.text = element_text(size = 9, colour = "black"),
          axis.title = element_text(size = 10),
          legend.position = "top")+
    labs(x="Called SNPs", y="", fill="SNP\nper\nPatient")+scale_fill_brewer(palette="Blues", direction = -1)
  wrap_plots(p1, p2)


  
  #Import the GTF file and get the additional features with RCAS
  gff <- importGtf(filePath = args[1])
  txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
  
    
  sites <- sites[sites$chromosome %in% unique(data.frame(gff[,1])$seqnames), ]
  queryRegions <- makeGRangesFromDataFrame(sites, keep.extra.columns=TRUE)
  transcript_features <- data.frame()
  for(i in 1:length(unique(queryRegions$SNP_type))){
    tmp <- subset(queryRegions, SNP_type ==i)
    tmp <- as.data.frame(summarizeQueryRegions(queryRegions = tmp, txdbFeatures = txdbFeatures))
    tmp$features <- rownames(tmp)
    tmp$SNP_type <- as.character(unique(queryRegions$SNP_type)[i])
    transcript_features = rbind(transcript_features , tmp)
  }



  transcript_features <- transcript_features %>% filter(!features %in% c("NoFeatures", "promoters", "transcripts", "exons"))
  transcript_features$features <- factor(transcript_features$features, levels = c("fiveUTRs", "cds", "introns", "threeUTRs"))

  transcript_features <- transcript_features %>% filter(features %in% c("fiveUTRs", "cds", "introns", "threeUTRs"))
  transcript_features <- transcript_features %>% mutate(features = recode(features,
                                                                        "fiveUTRs" = "5'-UTR",
                                                                        "cds" = "CDS",
                                                                        "introns" = "Intron",
                                                                        "threeUTRs" = "3'-UTR"))
  transcript_features$features <- factor(transcript_features$features, levels = c("5'-UTR", "CDS", "Intron", "3'-UTR"))

  p3 <- ggplot(transcript_features, aes(x= SNP_type, y = count, fill=SNP_type)) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(cols = vars(features))+theme(panel.background = element_blank(), legend.position = "", axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10))+
  labs(x="SNPs per Patient", y="Patient SNPs", fill="")+scale_fill_brewer(palette="Blues", direction = -1)



  Patient_colors <- c("#195E83", "#DFC8A2", "#E07B39", "#69BDD2", "#80391E", "#606060", "gray")
  pca_res <- prcomp(df_Patient, scale. = TRUE)
  res <- pca_res$rotation
  p4 <- ggplot(data.frame(res), aes(x=PC1, y=PC2, color=rownames(res)))+geom_point(size=5) + coord_equal()+
  theme(panel.background = element_blank(), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "right")+coord_fixed(ratio = 1)+
  scale_color_manual(values = Patient_colors)+
  labs(color="")
  wrap_plots(p3, p4)
  dev.off()
