# combination OR recluster by ATAC+RNA signal 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
DefaultAssay(ORN)<-"raw_RNA"
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

# show a typical combination;
# select a beautiful track to show :
library(ggplot2)
library(gggenes)
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
library(RColorBrewer)
library(UpSetR)

# ONLY IN CLUSTER CO-EXP 

obj<-subset(ORN,idents="p4:12");
obj_features<- unique(dotplot_data[dotplot_data$id%in%"p4:12",]$features.plot)
## Upset plot for 4 coreceptor barcode 
# 1: Observation with Or2 

ORN_count<-obj@assays$RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%c(obj_features,"Or2")),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or2 = names(which(ORN_matrix[4,]>0)), 
        Or43a_b = names(which(ORN_matrix[2,]>0)), 
        Or151 = names(which(ORN_matrix[1,]>0)), 
        Or152 = names(which(ORN_matrix[3,]>0)))
pdf("./Or151-152-43A-upsetR_Observation_withOrco.pdf", width=6, height=4)
upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()

listInput <- list(
        #Or2 = names(which(ORN_matrix[4,]>0)), 
        Or43a_b = names(which(ORN_matrix[2,]>0)), 
        Or151 = names(which(ORN_matrix[1,]>0)), 
        Or152 = names(which(ORN_matrix[3,]>0)))
pdf("./Or151-152-43A-upsetR_Observation_withoutOrco.pdf", width=6, height=4)
upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()

# All cell exp OR151 152,43a
## Upset plot for 4 coreceptor barcode 
# 1: Observation with Or2 

ORN_count<-ORN@assays$RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%c(obj_features,"Or2")),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or2 = names(which(ORN_matrix[4,]>0)), 
        Or43a_b = names(which(ORN_matrix[2,]>0)), 
        Or151 = names(which(ORN_matrix[1,]>0)), 
        Or152 = names(which(ORN_matrix[3,]>0)))
pdf("./All_cell-Or151-152-43A-upsetR_Observation_withOrco.pdf", width=6, height=4)
upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()

listInput <- list(
        #Or2 = names(which(ORN_matrix[4,]>0)), 
        Or43a_b = names(which(ORN_matrix[2,]>0)), 
        Or151 = names(which(ORN_matrix[1,]>0)), 
        Or152 = names(which(ORN_matrix[3,]>0)))
pdf("./All_cell-Or151-152-43A-upsetR_Observation_withoutOrco.pdf", width=6, height=4)
upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()
