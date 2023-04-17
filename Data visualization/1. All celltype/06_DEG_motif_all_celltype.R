library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(RColorBrewer)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
set.seed(1234)
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")

# find the DEG among all celltypes and get the specifical TF in ORNs

