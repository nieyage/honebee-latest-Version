## create Seurat by gex matrix 
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(DropletUtils)
library(Matrix)
library(DoubletFinder)
set.seed(1234)
bin20_obj<- readRDS("/md01/nieyg/project/honeybee/data/Spatial_data/bin20.seurat.rds")
bin50_obj<- readRDS("/md01/nieyg/project/honeybee/data/Spatial_data/bin50.seurat.rds")
bin20_obj@meta.data[1:4,]


