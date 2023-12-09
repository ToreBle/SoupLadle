library(GenomicFeatures)
library(RCAS)
library(Seurat)
library(ggplot2)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

#Import the GTF file and get the additional features with RCAS
gff <- importGtf(filePath = args[1])
txdbFeatures <- getTxdbFeaturesFromGRanges(gff)


pdf(paste(args[2], "/SNP_evaluation/Transcriptional_features_bulk_WES_", args[3],".pdf", sep=""), height = 3, width = 4)
transcript_features <- data.frame()
for(data_type in c("bulk_RNA", "WES")){
  print(data_type)
  sites <- read.table(paste(args[2], "/SNP_evaluation/", data_type, "_isec/sites.txt", sep=""))
  print(data_type)
  print(head(sites))
  sites$end = sites$V2+1
  #sites <- sites[,c(1,6,2,3,4,5)] check
  sites <- sites[,c(1,2,6,3,4,5)]
  colnames(sites) <- c("chromosome", "start", "end", "REF", "ALT", "Sample")
  sites$SNP_type <- sapply(sites$Sample, function(x) sum(as.numeric(strsplit(as.character(x), "")[[1]])))
  
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
  
 #For PBMC

  df_long <- pivot_longer(sites, cols = starts_with("Patient_"), names_to = "Samples", values_to = "Value") %>%
    filter(Value == "1") %>% group_by(SNP_type, Samples) %>% dplyr::count()
  
  p <- ggplot(df_long, aes(x = Samples, y = n, fill = factor(SNP_type))) +
    geom_bar(position="stack", stat = "identity") +
    theme(panel.background = element_blank(),
          legend.title = element_text(size = 9, colour = "black"),
          axis.text = element_text(size = 9, colour = "black"),
          axis.title = element_text(size = 10),
          legend.position = "top")+
    labs(x="", y="Called SNPs", fill="SNP\nper\nPatient", title=data_type)+scale_fill_brewer(palette="Blues", direction = -1)
  
  print(p)
  queryRegions <- makeGRangesFromDataFrame(sites, keep.extra.columns=TRUE)
  
  
  for(i in 1:length(unique(queryRegions$SNP_type))){
    tmp <- subset(queryRegions, SNP_type ==i)
    tmp <- as.data.frame(summarizeQueryRegions(queryRegions = tmp, txdbFeatures = txdbFeatures))
    tmp$features <- rownames(tmp)
    tmp$data_type <- data_type
    tmp$SNP_type <- as.character(unique(queryRegions$SNP_type)[i])
    transcript_features = rbind(transcript_features , tmp)
  }
}


transcript_features <- transcript_features %>% filter(!features %in% c("NoFeatures", "promoters", "transcripts", "exons"))
transcript_features$features <- factor(transcript_features$features, levels = c("fiveUTRs", "cds", "introns", "threeUTRs"))
transcript_features <- transcript_features %>% mutate(features = recode(features,
                          "fiveUTRs" = "5'-UTR",
                          "cds" = "CDS",
                          "introns" = "Intron",
                          "threeUTRs" = "3'-UTR"))
transcript_features <- transcript_features %>% mutate(data_type = recode(data_type,
                                                                        "WES" = "WES",
                                                                        "bulk_RNA" = "bulk RNA",))

ggplot(transcript_features, aes(x= data_type, y = count, fill=SNP_type)) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(cols = vars(features))+theme(panel.background = element_blank(), legend.title = element_text(size = 8, colour = "black"), axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 8))+
  labs(x="", y="Patient SNPs", fill="SNP\nper\nPatient")+scale_fill_brewer(palette="Blues", direction = -1)
dev.off()





