library(ggplot2)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

setwd("/data/SNP_Multiplexing/output/KH30_unfiltered/")

CMO <- read_csv("/data/SNP_Multiplexing/data/multiplex/KH30/outs/multi/multiplexing_analysis/assignment_confidence_table.csv")
CMO <- CMO %>% mutate(Assignment = ifelse(grepl("Blank", Assignment), "Unassigned", Assignment))
CMO$Method <- "CellPlex"
CMO$CMO_assignment <- CMO$Assignment[match(CMO$Barcodes, CMO$Barcodes)]


CMO <- CMO %>%
  mutate(CMO_assignment = recode(CMO_assignment,
                                 "PBMC1" = "Patient_1",
                                 "PBMC2" = "Patient_2",
                                 "PBMC3" = "Patient_3",
                                 "PBMC4" = "Patient_4",
                                 "PBMC5" = "Patient_5"))


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



process_souporcell_clusters <- function(file_path, method_name, target_column, assign_to_cluster = FALSE) {
  soup_data <- read_tsv(file_path) %>%
    mutate(assignment = ifelse(grepl("\\d+/\\d+", assignment), status, assignment)) %>%
    mutate(Method = method_name)
  
  CMO[[target_column]] <- soup_data$assignment[match(CMO$Barcodes, soup_data$barcode)]
  CMO[[target_column]] <- recode(CMO[[target_column]], !!!c("unassigned" = "Unassigned", "doublet" = "Multiplet"))
  
  if (assign_to_cluster) {
    assign_patient_to_cluster(CMO, target_column)
  } else {
    return(CMO)
  }
}

CMO <- process_souporcell_clusters("./output_souporcell_SNP_only/clusters.tsv", "SouporCell - SNP", "Souporcell", assign_to_cluster = TRUE)
CMO <- process_souporcell_clusters("./output_souporcell_minimap/clusters.tsv", "SouporCell - minimap2", "Soup_Minimap2", assign_to_cluster = TRUE)

process_vireo_data <- function(file_path, method_name, target_column, assign_to_cluster = FALSE) {
  vireo_data <- read_tsv(file_path) %>%
    mutate(Method = method_name)
  
  CMO[[target_column]] <- vireo_data$donor_id[match(CMO$Barcodes, vireo_data$cell)]
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
CMO <- process_vireo_data("./vireo_vartrix/donor_ids.tsv", "Vireo - vartrix", "VarTrix_Vireo", assign_to_cluster = TRUE)


CMO <- CMO %>%
  mutate(Vireo_WES = recode(Vireo_WES,
                            "M34287" = "Patient_1",
                            "M34288" = "Patient_2",
                            "M34289" = "Patient_3",
                            "M34290" = "Patient_4",
                            "M34291" = "Patient_5"),
         Vireo_bulkRNA = recode(Vireo_bulkRNA,
                                "KH_75" = "Patient_1",
                                "KH_77" = "Patient_2",
                                "KH_79" = "Patient_3",
                                "KH_76" = "Patient_4",
                                "KH_78" = "Patient_5"),
  )



meta.data <- CMO %>% select(Barcodes, Method,CMO_assignment,
                            Souporcell, Soup_Minimap2,
                            Vireo, VarTrix_Vireo,
                            Vireo_WES, Vireo_bulkRNA)
meta.data <- data.frame(meta.data, row.names = meta.data$Barcodes)
meta.data$Barcodes <- NULL

pdf("CellPlex_PBMC_Heatmaps.pdf", height = 6, width = 6)
for (i in 3:ncol(meta.data)){
  group_counts <- table(meta.data$CMO_assignment, meta.data[,i])
  group_counts <- data.frame(group_counts)
  print(ggplot(data = group_counts, aes(x=Var1, y=Var2, fill=Freq)) + geom_tile()+coord_fixed()+scale_fill_gradient(low="white", high="darkred")+
          theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=12, colour = "black"))+RotatedAxis()+
          labs(title=paste("Comparison of CellPlex vs. ", colnames(meta.data[i]), sep=""), y=colnames(meta.data[i]), x="CellPlex Sample", fill="Cells"))
}
dev.off()


# Process single-cell RNA-seq data by loading the 10X matrix
PBMC <- Read10X("/data/SNP_Multiplexing/data/multiplex/KH30/outs/multi/count/raw_feature_bc_matrix/")
PBMC <- CreateSeuratObject(counts = PBMC$`Gene Expression`, meta.data = meta.data, project = "PBMC", )

cells <- read_tsv("all_cell_barcodes.tsv", col_names = "Barcode")
PBMC <- subset(PBMC, cells = cells$Barcode)

set.seed(1989)
PBMC <- NormalizeData(PBMC)
PBMC <- FindVariableFeatures(PBMC, selection.method = "vst", nfeatures = 4000)
PBMC <- ScaleData(PBMC)
PBMC <- RunPCA(PBMC, features = VariableFeatures(object = PBMC))
PBMC <- FindNeighbors(PBMC)
PBMC <- RunUMAP(PBMC, dims = 1:20)
PBMC <- FindClusters(PBMC, resolution = 0.1)


# Cluster assignment based on marker gene expression
new.cluster.ids <- c("CD4+ T cells", "CD14+ Monocytes", "CD16+ Granulo", "CD19+ B cells", "CD8+ T cells", "CD8+ T cells",
                     "Mixed", "Mixed", "Mixed", "Mixed", "Mixed", "Mixed")
names(new.cluster.ids) <- levels(PBMC)
PBMC <- RenameIdents(PBMC, new.cluster.ids)
PBMC$Cluster <- Idents(PBMC)


#single-cell QC and marker gene expression
pdf("QC_plots_single_cell_PBMC.pdf", width = 12, height = 6)
PBMC[["percent.mt"]] <- PercentageFeatureSet(PBMC, pattern = "^MT-")
p1 <- VlnPlot(PBMC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = Patient_colors)
p2 <- DotPlot(PBMC, features = c("CD3G", "CD4", "CD14", "MS4A7", "FCGR3A", "CD19", "MS4A1", "CD8A"), group.by = "Cluster")+RotatedAxis()+labs(x="")+coord_flip()
wrap_plots(p1, p2, widths = c(3,1))
dev.off()

saveRDS(PBMC, "KH30_PBMC_Multiplex.rds")

###

# Function to simplify factor level assignments
same_levels <- function(column) {
  factor(column, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Multiplet", "Unassigned"))
}

# Apply the function to specific columns in the PBMC dataframe
PBMC$CMO_assignment <- same_levels(PBMC$CMO_assignment)
PBMC$Souporcell <- same_levels(PBMC$Souporcell)
PBMC$Soup_Minimap2 <- same_levels(PBMC$Soup_Minimap2)
PBMC$Vireo <- same_levels(PBMC$Vireo)
PBMC$VarTrix_Vireo <- same_levels(PBMC$VarTrix_Vireo)
PBMC$Vireo_bulkRNA <- same_levels(PBMC$Vireo_bulkRNA)
PBMC$Vireo_WES <- same_levels(PBMC$Vireo_WES)

#PBMC$CMO_assignment <- factor(PBMC$CMO_assignment, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Unassigned", "Multiplet"))
#PBMC$Soup_SNP <- factor(PBMC$Soup_SNP, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Unassigned", "Multiplet"))

PBMC$Cluster <- factor(PBMC$Cluster, levels = c("CD14+ Monocytes", "CD16+ Granulo", "CD4+ T cells", "CD8+ T cells", "CD19+ B cells", "Mixed"))
Patient_colors <- c("#195E83", "#DFC8A2", "#E07B39", "#69BDD2", "#80391E", "#606060", "gray")




#UMAP of the different demultiplexing approaches
pdf("UMAP_Demultiplexing_PBMC.pdf", bg = "transparent", height = 6, width = 6)
DimPlot(PBMC, reduction = "umap", group.by = "CMO_assignment")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "Souporcell")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "Soup_Minimap2")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "Vireo")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "VarTrix_Vireo")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "Vireo_bulkRNA")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "Vireo_WES")+theme(aspect.ratio=1)+scale_color_manual(values = Patient_colors)
DimPlot(PBMC, reduction = "umap", group.by = "Cluster")+theme(aspect.ratio=1, )+scale_color_manual(values = Patient_colors)
dev.off()


pdf("UMAP_PBMC_Demultiplexing2.pdf", height = 6, width = 18)
p1 <- DimPlot(PBMC, reduction = "umap", group.by = "Cluster")+theme(aspect.ratio=1, legend.position = "bottom")+scale_color_manual(values = Patient_colors)+guides(color=guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=5)))
p2 <- DimPlot(PBMC, reduction = "umap", group.by = "CMO_assignment")+theme(aspect.ratio=1, legend.position = "bottom")+scale_color_manual(values = Patient_colors)+guides(color=guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=5)))
p3 <- DimPlot(PBMC, reduction = "umap", group.by = "Souporcell")+theme(aspect.ratio=1, legend.position = "bottom")+scale_color_manual(values = Patient_colors)+guides(color=guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=5)))
ggarrange(p1,ggplot()+theme_void(),p2,ggplot()+theme_void(),p3, nrow =1, widths = c(1, -0.5, 1, -0.5, 1))
dev.off()



pdf("Cluster_and_Patient_Plots_PBMC.pdf", height = 4, width = 12)
df <- PBMC@meta.data
df$Barcode <- row.names(PBMC@meta.data)

df <- df %>% dplyr::select(Barcode, Cluster, CMO_assignment:Vireo_bulkRNA)
df$Cluster <- factor(df$Cluster, levels = c("CD14+ Monocytes", "CD16+ Granulo", "CD4+ T cells", "CD8+ T cells", "CD19+ B cells", "Mixed"))


#Per Cluster Plot
df_Cluster <- df %>%
  pivot_longer(cols = c(CMO_assignment:Vireo_bulkRNA), names_to = "Name", values_to = "Patient") %>%
  group_by(Cluster) %>% add_count(Name, name = "Clustering")
df_Cluster$Patient <- factor(df_Cluster$Patient, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Multiplet", "Unassigned"))
ggplot(df_Cluster, aes(y=Clustering, x=Name, fill=Patient)) + geom_bar(stat = "identity")+facet_grid(.~Cluster, scales = "free")+RotatedAxis()+
  theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=12, colour = "black"))+RotatedAxis()+
  labs(y="Cells", x="Method", fill="Sample")+scale_fill_manual(values = Patient_colors)

#Per Patient plot
df_Patient <- df %>%
  pivot_longer(cols = c(CMO_assignment:Vireo_bulkRNA), names_to = "Name", values_to = "Value") %>%
  group_by(Value) %>% add_count(Name, name = "Patient")
df_Patient$Value <- factor(df_Patient$Value, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Multiplet", "Unassigned"))
ggplot(df_Patient, aes(y=Patient, x=Name, fill=Value)) + geom_bar(stat = "identity")+facet_grid(.~Value, scales = "free")+RotatedAxis()+
  theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=10, colour = "black"), legend.position = "")+RotatedAxis()+
  labs(y="Cells", x="Method", fill="Sample")+scale_fill_manual(values = Patient_colors)

df_Method <- df %>%
  pivot_longer(cols = c(CMO_assignment:Vireo_bulkRNA), names_to = "Method", values_to = "Patient") %>%
  group_by(Cluster, Method, Patient) %>% dplyr::count()

df_Method$Patient <- factor(df_Method$Patient, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Multiplet", "Unassigned"))
df_Method$Method <- factor(df_Method$Method, levels = c("CMO_assignment", "Souporcell", "Soup_Minimap2", "Vireo", "VarTrix_Vireo", "Vireo_bulkRNA", "Vireo_WES"))

ggplot(df_Method, aes(x=Patient, y=n, fill=Cluster)) + geom_bar(stat = "identity")+facet_grid(.~Method, scales = "free")+RotatedAxis()+
  theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=10, colour = "black"))+RotatedAxis()+
  labs(y="Cells", x="")+scale_fill_manual(values = Patient_colors)

ggplot(df_Method, aes(x = Cluster, y = n, fill = Patient)) + geom_bar(stat = "identity", position="fill")+facet_grid(.~Method, scales = "free")+RotatedAxis()+
  theme(panel.background = element_blank(), axis.title=element_text(size=14), axis.text=element_text(size=10, colour = "black"))+RotatedAxis()+
  labs(y="Proportion of single cells", x="")+scale_fill_manual(values = Patient_colors)

dev.off()


###### Evaluation of the assignments


pdf("Assignment_evaluation_PBMC.pdf")

# Precision - Recall


df <- df %>% mutate(Clustering = recode(Cluster,
                                        "CD14+ Monocytes" = "Patient_1",
                                        "CD16+ Granulo" = "Patient_2",
                                        "CD4+ T cells" = "Patient_3",
                                        "CD8+ T cells" = "Patient_4",
                                        "CD19+ B cells" = "Patient_5",
                                        "Mixed" = "Multiplet",
                                        "Unassigned" = "Unassigned"))

df <- df %>% dplyr::select(CMO_assignment:Vireo_bulkRNA, Clustering)

# Clustering as ground truth
true_labels <- df$Clustering
true_labels <- factor(true_labels, levels = c("Patient_1", "Patient_2", "Patient_3", "Patient_4", "Patient_5", "Multiplet", "Unassigned"))

# Predict the different method performances
predicted_df <- df
predicted_df$Clustering <- NULL

# Initialize variables to store precision and recall for each model
models <- colnames(predicted_df)
precisions <- vector("numeric", length(models))
recalls <- vector("numeric", length(models))

# Calculate precision and recall for each model
for (i in 1:length(models)) {
  model_predictions <- predicted_df[, i]
  TP <- sum(true_labels == model_predictions) #& true_labels == unique(true_labels))
  FP <- sum(model_predictions %in% unique(true_labels)) - TP
  FN <- sum(true_labels %in% unique(model_predictions)) - TP
  
  precisions[i] <- TP / (TP + FP)
  recalls[i] <- TP / (TP + FN)
}

# Create a summary data frame
results_df <- data.frame(Method = models, Precision = precisions, Recall = recalls)


Category <- c("CMO_assignment" = "Labeling", "Souporcell" = "Deconvolution", "Soup_Minimap2" = "Deconvolution",
              "Vireo" = "Deconvolution", "VarTrix_Vireo" = "Deconvolution",
              "Vireo_WES" = "Assignment", "Vireo_bulkRNA" = "Assignment")
results_df <- results_df %>%
  mutate(Call = Category[Method])
#results_df$Call <-  c("Labeling", "Deconvolution", "Deconvolution", "Assignment", "Assignment", "Deconvolution", "Deconvolution")
results_df$Call <- factor(results_df$Call, levels = c("Labeling", "Deconvolution", "Assignment"))
results_df$Method <- factor(results_df$Method, levels = c("CMO_assignment", "Souporcell", "Soup_Minimap2", "Vireo", "VarTrix_Vireo", "Vireo_bulkRNA", "Vireo_WES"))

ggplot(results_df, aes(x=Precision, Recall, color = Method)) +
  geom_point(size=5, aes(shape=Call)) + coord_equal()+
  theme(panel.background = element_blank(), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.position = "right")+coord_fixed(ratio = 1)+scale_color_brewer(palette = "Paired", direction = -1)

####

# Comparision of the Recall for the different methods
df <- df %>% dplyr::select(CMO_assignment:Vireo_bulkRNA, Clustering)

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

ggplot(co_occurrence_df %>% filter(Var1=="Clustering"), aes(x=Var2, y = Freq/100, fill=Var2)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Reds", direction = 1) +
  scale_y_continuous(labels = scales::percent)+
  labs(title = "Overlap of assigned cells only", x = "", y = "", fill="matched\nassignments") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 10),
        legend.position = "right") +
  geom_text(aes(label = sprintf("%.1f%%", Freq)), vjust = 3.5, size = 2)

#Version2
ggplot(co_occurrence_df, aes(Var1, Var2, fill = Freq)) +
  geom_tile() + scale_fill_distiller(palette = "RdBu", direction = -1) +
  labs(title = "Overlap of assigned cells only", x = "", y = "", fill="matched\nassignments") +
  theme_minimal() + coord_equal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        axis.title = element_text(size = 10),
        legend.position = "right")+
  geom_text(aes(label = Freq), vjust = 0.5, size = 2)

dev.off()


