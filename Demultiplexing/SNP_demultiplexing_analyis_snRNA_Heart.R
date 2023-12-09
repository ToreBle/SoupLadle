library(ggplot2)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

setwd("/data/SNP_Multiplexing/output/KH108_unfiltered/")

# CMO assignment from CellRanger
CMO <- read_csv("/data/SNP_Multiplexing/data/multiplex/KH108/outs/multi/multiplexing_analysis/assignment_confidence_table.csv")
CMO <- CMO %>% mutate(Assignment = ifelse(grepl("Blank", Assignment), "Unassigned", Assignment))
CMO$Method <- "CellPlex"
CMO$CMO_assignment <- CMO$Assignment[match(CMO$Barcode, CMO$Barcode)]

# Hashtag assignment from CellRanger
HTO <- read_csv("/data/SNP_Multiplexing/data/multiplex/KH108_HTO/outs/multi/multiplexing_analysis/assignment_confidence_table.csv")
HTO <- HTO %>% mutate(Assignment = ifelse(grepl("Blank", Assignment), "Unassigned", Assignment))
HTO$Method <- "Hashtag"
CMO$HTO_assignment <- HTO$Assignment[match(CMO$Barcode, HTO$Barcode)]


CMO <- CMO %>%
  mutate(CMO_assignment = recode(CMO_assignment,
                                 "CMO308" = "Patient_38",
                                 "CMO309" = "Patient_41",
                                 "CMO310" = "Patient_42",
                                 "CMO311" = "Patient_70",
                                 "CMO312" = "Patient_78"),
         HTO_assignment = recode(HTO_assignment,
                                 "HTO_A0451" = "Patient_38",
                                 "HTO_A0452" = "Patient_41",
                                 "HTO_A0453" = "Patient_42",
                                 "HTO_A0454" = "Patient_70",
                                 "HTO_A0455" = "Patient_78"))


# Function to assign the clusters with the CMO demultiplexed data
assign_patient_to_cluster <- function(data_frame, column_name) {
  vireo_singlet <- subset(data_frame, !data_frame[[column_name]] %in% c("Unassigned", "Multiplet"))
  cluster_assignment <- data.frame(table(vireo_singlet[[column_name]], vireo_singlet$CMO_assignment))
  cluster_assignment <- subset(cluster_assignment, Var2 != "Unassigned" & Var2 != "Multiplet")
  cluster_table <- cluster_assignment %>% group_by(Var1) %>% top_n(1, Freq)
  #Multiple best hits: cluster_table$Var2 = paste(cluster_table$Var1, "_", cluster_table$Var2, sep="")
  
  for (i in 0:nrow(cluster_table)) {
    data_frame[[column_name]][data_frame[[column_name]] == i] <- as.character(subset(cluster_table, Var1 == i)$Var2)
  }
  return(data_frame)
}

# Assignment of SouporCell Clusters
process_souporcell_clusters <- function(file_path, method_name, target_column, assign_to_cluster = FALSE) {
  soup_data <- read_tsv(file_path) %>%
    mutate(assignment = ifelse(grepl("\\d+/\\d+", assignment), status, assignment)) %>%
    mutate(Method = method_name)
  
  CMO[[target_column]] <- soup_data$assignment[match(CMO$Barcode, soup_data$barcode)]
  CMO[[target_column]] <- recode(CMO[[target_column]], !!!c("unassigned" = "Unassigned", "doublet" = "Multiplet"))
  
  
  if (assign_to_cluster) {
    assign_patient_to_cluster(CMO, target_column)
  } else {
    return(CMO)
  }
}
#  assign_patient_to_cluster(CMO, target_column)
#}


CMO <- process_souporcell_clusters("./output_souporcell_SNP_only/clusters.tsv", "SouporCell - SNP", "Souporcell", assign_to_cluster = TRUE)
CMO <- process_souporcell_clusters("./output_souporcell_minimap/clusters.tsv", "SouporCell - minimap2", "Soup_Minimap2", assign_to_cluster = TRUE)

# Assignment of Vireo Clusters
process_vireo_data <- function(file_path, method_name, target_column, assign_to_cluster = FALSE) {
  vireo_data <- read_tsv(file_path) %>%
    mutate(Method = method_name)
  
  CMO[[target_column]] <- vireo_data$donor_id[match(CMO$Barcode, vireo_data$cell)]
  CMO[[target_column]] <- gsub("^donor\\s*", "", CMO[[target_column]])
  CMO[[target_column]] <- recode(CMO[[target_column]], !!!c("unassigned" = "Unassigned", "doublet" = "Multiplet"))
  
  if (assign_to_cluster) {
    assign_patient_to_cluster(CMO, target_column)
  } else {
    return(CMO)
  }
}

# Update the CMO data frame with each dataset and assign to cluster for SNP only based methods
CMO <- process_vireo_data("./vireo_WES/donor_ids.tsv", "Vireo - WES", "Vireo_WES")
CMO <- process_vireo_data("./vireo_bulkRNA/donor_ids.tsv", "Vireo - bulkRNA", "Vireo_bulkRNA")
CMO <- process_vireo_data("./vireo_SNP_only/donor_ids.tsv", "Vireo - SNP", "Vireo", assign_to_cluster = TRUE)
CMO <- process_vireo_data("./vartrix_vireo_bulkRNA/donor_ids.tsv", "Vireo - vartrix", "VarTrix_Vireo_bulkRNA", assign_to_cluster = TRUE)

CMO <- CMO %>%
  mutate(Vireo_WES = recode(Vireo_WES,
                       "M36483Ex" = "Patient_38",
                       "M36484Ex" = "Patient_41",
                       "M36485Ex" = "Patient_42",
                       "M36486Ex" = "Patient_70",
                       "M36487Ex" = "Patient_78"),
  Vireo_bulkRNA = recode(Vireo_bulkRNA,
                           "KH113" = "Patient_38",
                           "KH114" = "Patient_41",
                           "KH115" = "Patient_42",
                           "KH116" = "Patient_70",
                           "KH117" = "Patient_78"),
  VarTrix_Vireo_bulkRNA = recode(VarTrix_Vireo_bulkRNA,
                         "KH113" = "Patient_38",
                         "KH114" = "Patient_41",
                         "KH115" = "Patient_42",
                         "KH116" = "Patient_70",
                         "KH117" = "Patient_78")
  )

meta.data <- CMO %>% select(Barcode, Assignment_Probability, Method, CMO_assignment, HTO_assignment,
                            Souporcell, Soup_Minimap2, Vireo,
                            Vireo_WES, Vireo_bulkRNA, VarTrix_Vireo_bulkRNA)

meta.data <- CMO %>% select(Barcode, Assignment_Probability, Method, CMO_assignment, HTO_assignment,
                            Souporcell, VarTrix_Vireo_bulkRNA,
                            Vireo_bulkRNA, Vireo_WES)

meta.data <- data.frame(meta.data, row.names = meta.data$Barcode)
meta.data$Barcode <- NULL

pdf("CellPlex_Heart_Heatmaps.pdf", height = 6, width = 6)
for (i in 3:ncol(meta.data)){
  group_counts <- table(meta.data$CMO_assignment, meta.data[,i])
  group_counts <- data.frame(group_counts)
  print(ggplot(data = group_counts, aes(x=Var1, y=Var2, fill=Freq)) + geom_tile()+coord_fixed()+scale_fill_gradient(low="white", high="#900C3F")+
    theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=12, colour = "black"))+RotatedAxis()+
    labs(title=paste("Comparison of CellPlex vs. ", colnames(meta.data[i]), sep=""), y=colnames(meta.data[i]), x="CellPlex Sample", fill="Cells"))
}
dev.off()


Heart <- Read10X("/data/SNP_Multiplexing/data/multiplex/KH108/outs/multi/count/raw_feature_bc_matrix/")
Heart <- CreateSeuratObject(counts = Heart$`Gene Expression`, meta.data = meta.data, project = "Heart", )

cells <- read_tsv("all_cell_barcodes.tsv", col_names = "Barcode")
Heart <- subset(Heart, cells = cells$Barcode)

set.seed(1989)
Heart <- NormalizeData(Heart)
Heart <- FindVariableFeatures(Heart, selection.method = "vst", nfeatures = 4000)
Heart <- ScaleData(Heart)
Heart <- RunPCA(Heart, features = VariableFeatures(object = Heart))
Heart <- FindNeighbors(Heart)
Heart <- RunUMAP(Heart, dims = 1:20)
Heart <- FindClusters(Heart, resolution = 0.2)

new.cluster.ids <- c("Cardiomyocyte", "Cardiomyocyte", "Cardiomyocyte", "Fibroblast", "Cardiomyocyte", "Endothelial",
                     "Immune")
names(new.cluster.ids) <- levels(Heart)
Heart <- RenameIdents(Heart, new.cluster.ids)
Heart$Cluster <- Idents(Heart)


saveRDS(Heart, "KH108_Heart_Multiplex.rds")



####UMAP###

same_levels <- function(column) {
  factor(column, levels = c("Patient_38", "Patient_41", "Patient_42", "Patient_70", "Patient_78", "Multiplet", "Unassigned"))
}

# Apply the function to specific columns in the Heart dataframe
Heart$CMO_assignment <- same_levels(Heart$CMO_assignment)
Heart$HTO_assignment <- same_levels(Heart$HTO_assignment)
Heart$Souporcell <- same_levels(Heart$Souporcell)
Heart$Soup_Minimap2 <- same_levels(Heart$Soup_Minimap2)
Heart$Vireo <- same_levels(Heart$Vireo)
Heart$Vireo_WES <- same_levels(Heart$Vireo_WES)
Heart$Vireo_bulkRNA <- same_levels(Heart$Vireo_bulkRNA)
Heart$VarTrix_Vireo_bulkRNA <- same_levels(Heart$VarTrix_Vireo_bulkRNA)

Heart <- subset(Heart, subset = HTO_assignment != "NA")

pdf("UMAP_Heart_Demultiplexing.pdf", bg = "transparent", height = 6, width = 6)
Patient_colors <- c("#195E83", "#DFC8A2", "#E07B39", "#69BDD2", "#80391E", "#606060", "gray")
#Patient_colors <- c(brewer.pal(n = 5, name = 'OrRd'), "#606060", "gray")
DimPlot(Heart, reduction = "umap", group.by = "CMO_assignment")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(Heart, reduction = "umap", group.by = "HTO_assignment")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(Heart, reduction = "umap", group.by = "Souporcell")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
#DimPlot(Heart, reduction = "umap", group.by = "Soup_Minimap2")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
#DimPlot(Heart, reduction = "umap", group.by = "Vireo")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(Heart, reduction = "umap", group.by = "Vireo_WES")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(Heart, reduction = "umap", group.by = "Vireo_bulkRNA")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(Heart, reduction = "umap", group.by = "VarTrix_Vireo_bulkRNA")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(Heart, reduction = "umap", group.by = "Cluster")+theme(aspect.ratio=1)+scale_color_brewer(palette = "Dark2")
dev.off()

####

pdf("QC_plots_Heart.pdf", width = 12, height = 6)
DefaultAssay(Heart) <- "RNA"
Heart[["percent.mt"]] <- PercentageFeatureSet(Heart, pattern = "^MT-")
p1 <- VlnPlot(Heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = Patient_colors)
p2 <- DotPlot(Heart, features = c("TNNT2", "MYH6", "DCN", "FBLN2", "PECAM1", "VWF", "CD163", "CD86"), group.by = "Cluster")+RotatedAxis()+labs(x="")+coord_flip()
wrap_plots(p1, p2, widths = c(3,1))
dev.off()

#####

pdf("Demultiplexing_Method_Comparision_Heart.pdf", height = 4, width = 12)
df <- Heart@meta.data
df$Barcode <- row.names(Heart@meta.data)
df <- df %>% select(Barcode, Cluster, CMO_assignment, HTO_assignment, Souporcell, VarTrix_Vireo_bulkRNA, Vireo_bulkRNA, Vireo_WES) %>% filter(!is.na(HTO_assignment))

df_Method <- df %>%
  pivot_longer(cols = c(CMO_assignment:Vireo_WES), names_to = "Method", values_to = "Patient") %>%
  group_by(Cluster, Method, Patient) %>% count()

df_Method$Patient <- factor(df_Method$Patient, levels = c("Patient_38", "Patient_41", "Patient_42", "Patient_70", "Patient_78", "Multiplet", "Unassigned"))
#df_Method$Method <- factor(df_Method$Method, levels = c("CMO_assignment","HTO_assignment",  "Souporcell", "Soup_Minimap2", "Vireo", "Vireo_bulkRNA", "Vireo_WES"))
#Patient_colors <- c("#195E83", "#DFC8A2", "#E07B39", "#69BDD2", "#80391E", "#606060", "gray")


ggplot(df_Method, aes(x = Cluster, y = n, fill = Patient)) + geom_bar(stat = "identity", position="fill")+facet_grid(.~Method, scales = "free")+RotatedAxis()+
  theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=10, colour = "black"))+RotatedAxis()+
  labs(y="Proportion of single cells", x="")+scale_fill_manual(values = Patient_colors)


df <- df %>% dplyr::select(CMO_assignment:Vireo_WES)

# Convert dataframe columns to character vectors
df_char <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

# Create a common set of levels from all columns
all_levels <- unique(unlist(df_char))

# Initialize a matrix to store co-occurrence percentages
num_columns <- ncol(df)
co_occurrence_matrix <- matrix(0, nrow = num_columns, ncol = num_columns, dimnames = list(names(df), names(df)))


# Function to calculate co-occurrence percentage of patients in two columns
calculate_co_occurrence_percentage <- function(patient1, patient2, levels) {
  intersection <- sum(patient1 == patient2)
  total_elements <- length(patient1)
  percentage <- (intersection / total_elements) * 100
  return(percentage)
}


# Calculate co-occurrence percentages for each column pair
for (i in 1:num_columns) {
  for (j in 1:num_columns) {
    co_occurrence_matrix[i, j] <- calculate_co_occurrence_percentage(df_char[, i], df_char[, j], all_levels)
  }
}

co_occurrence_df <- as.data.frame(as.table(co_occurrence_matrix))


ggplot(co_occurrence_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() + scale_fill_distiller(palette = "RdBu", direction = -1, labels = scales::percent_format(scale = 1)) +
  labs(title = "Co-occurrent assignments", x = "", y = "", fill="matched\nassignments") +
  theme_minimal() +coord_equal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 10),
        legend.position = "right") +
  geom_text(aes(label = sprintf("%.1f%%", Freq)), vjust = 0.5, size = 2)


# Only for Assigned Cells
calculate_co_occurrence_patients <- function(patient1, patient2, levels) {
  # Check if values start with "Pat"
  start_with_pat <- grepl("^Pat", patient1) & grepl("^Pat", patient2)
  
  # Filter values that start with "Pat"
  patient1_pat <- patient1[start_with_pat]
  patient2_pat <- patient2[start_with_pat]
  
  # Calculate intersection of values that start with "Pat"
  intersection <- sum(patient1_pat == patient2_pat)
  total_elements <- length(patient1_pat)
  
  # Handle the case where there are no values that start with "Pat"
  if (total_elements == 0) {
    return(0)
  }
  
  percentage <- (intersection / total_elements) * 100 #Version1
  #percentage <- intersection #Version2
  
  return(percentage)
}



# Calculate co-occurrence percentages for each column pair
for (i in 1:num_columns) {
  for (j in 1:num_columns) {
    co_occurrence_matrix[i, j] <- calculate_co_occurrence_patients(df_char[, i], df_char[, j], all_levels)
  }
}

co_occurrence_df <- as.data.frame(as.table(co_occurrence_matrix))

#Version1
ggplot(co_occurrence_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() + scale_fill_distiller(palette = "Reds", direction = 1, labels = scales::percent_format(scale = 1)) +
  labs(title = "Overlap of assigned cells only", x = "", y = "", fill="matched\nassignments") +
  theme_minimal() + coord_equal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 10),
        legend.position = "right") +
  geom_text(aes(label = sprintf("%.1f%%", Freq)), vjust = 0.5, size = 2)
dev.off()


pdf("Percentage_of_assigned_cells_Heart.pdf", height = 4, width = 4)
#Percentage of assigned single-cells
df <- df %>% dplyr::select(CMO_assignment:Vireo_WES)
percentage_df <- df %>%
  summarise_all(~ sum(grepl("^Pat", .)) / n() * 100)
percentage_df_long <- percentage_df %>%
  gather(key = "Method", value = "Percentage")


ggplot(percentage_df_long %>% mutate(Method = fct_reorder(Method, Percentage)), aes(x = Method, y = Percentage, fill=Method)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Blues", direction = 1)+
  theme(panel.background = element_blank(),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.position = "") +
  labs(x = "", y = "Patient assigned cells [%]") +
  coord_flip()
dev.off()

