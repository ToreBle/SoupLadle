library(Seurat)
library(ggplot2)
library(tidyverse)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

sites <- read.table(paste(args[1], "sites.txt", sep=""))
sites$end = sites$V2+1
sites <- sites[,c(1,2,6,3,4,5)]
colnames(sites) <- c("chromosome", "start", "end", "REF", "ALT", "Sample")
sites$SNP_type <- paste(sites$REF, sites$ALT, sep="->")

sites <- sites %>%
  mutate(
    Method = case_when(
      Sample == 1 ~ "VarTrix",
      Sample == 10 ~ "CellSNP",
      Sample == 11 ~ "Common",
      TRUE ~ as.character(Sample))) %>%
  mutate(
    Shared = case_when(
      Sample == 1 ~ "VarTrix",
      Sample == 10 ~ "CellSNP",
      Sample == 11 ~ "Common",
      TRUE ~ as.character(Sample))) %>%
  bind_rows(
    sites %>% filter(Sample == 11) %>%
      mutate(Shared = "VarTrix", Method = "Common"),
    sites %>% filter(Sample == 11) %>%
      mutate(Shared = "CellSNP", Method = "Common")) %>%
  filter(Shared != "Common") %>%
  mutate(SNPs = ifelse(grepl("^[A-Z]->[A-Z]$", SNP_type), SNP_type, "others"))

summary <- sites %>%
  dplyr::count(Method, Shared, SNPs)

pdf(paste("SNP_calling_cellSNP_Vartrix_", args[2],".pdf", sep = ""), height = 12, width = 6)
ggplot(summary, aes(x=n, y = SNPs, fill=Method)) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(rows = vars(Shared), scales = "free")+theme(panel.background = element_blank(), legend.title = element_text(size = 8, colour = "black"), axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 8))+
  labs(x="Called SNPs", y="Single-nucleotid Polymorphism (SNP)")+NoLegend()+scale_fill_brewer(palette="Reds", direction = 1)

ggplot(summary, aes(x=log2(n), y = SNPs, fill=Method)) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(rows = vars(Shared), scales = "free")+theme(panel.background = element_blank(), legend.title = element_text(size = 8, colour = "black"), axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 8))+
  labs(x="log2 (Called SNPs)", y="Single-nucleotid Polymorphism (SNP)")+NoLegend()+scale_fill_brewer(palette="Reds", direction = 1)

summary2 <- sites %>%
  dplyr::count(Method, Shared)
summary2$Method <- factor(summary2$Method, levels = c("VarTrix","CellSNP", "Common"))

ggplot(summary2, aes(x=Shared, y=n, fill=Method)) + 
  geom_bar(position="stack", stat="identity")+
  theme(panel.background = element_blank(), legend.title = element_text(size = 9, colour = "black"), legend.position = "top", axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 12))+
  labs(x="", y="Called SNPs")+scale_fill_brewer(palette="Blues", direction = 1)

#summary2$SNP <- "SNP"
#ggplot(summary2, aes(x=SNP, y = n, fill=Method)) + 
#  geom_bar(position="stack", stat="identity")+
#  theme(panel.background = element_blank(), legend.title = element_text(size = 9, colour = "black"), legend.position = "top", axis.text = element_text(size = 9, colour = "black"), axis.title = element_text(size = 10))+
#  labs(x="", y="# of Called SNPs")+scale_fill_brewer(palette="Blues", direction = 1)+scale_y_continuous(labels = unit_format(unit = "M", scale = 1e-6))

dev.off()
