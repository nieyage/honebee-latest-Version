library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(VennDiagram)
library(RColorBrewer)
set.seed(1234)
# remove nopower cluster,then plot UMAP
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )
# load DEP 
tau_data<-read.csv("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_tau-0.9_cluster_specfic_data_peak_ORN.csv",row.names=1)

# load DEG 
DEG <-read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEG_tau_cor_withoutGRN.csv",row.names=1)

# for cluster p2:21
cluster<- "p2:21"
peak_in_cluster<- tau_data[tau_data$cluster==cluster,]$peak;
# motif enrichment by TOBIAS 

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)

# Get a list of motif position frequency matrices from the JASPAR database
obj<- subset(ORN,idents=cluster)
pfm <-  readJASPARMatrix("/data/R04/liyh526/project/Honeybee/03_CisBP2JASPAR/honeybee_jaspar/CisBP-honeybee.jaspar", matrixClass=c("PFM", "PWM", "PWMProb"))

DefaultAssay(obj) <- 'peaks_ORN_subcluster'
# add motif information
obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = pfm
)

# test enrichment
enriched.motifs <- FindMotifs(
  object = obj,
  features = peak_in_cluster
)
MotifPlot(
  object = obj,
  motifs = head(rownames(enriched.motifs))
)

obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor
)

DefaultAssay(obj) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = obj,
  features = head(rownames(enriched.motifs)),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

# exp specfic in cluster 
gene_in_cluster<- DEG[DEG$cluster==cluster,]$gene;
TF_list_info<- read.table("./05_ORN_cluster2/07_DEG_and_DEP/TF_Information.txt",sep="\t")
ID_trans<- read.table("/md01/nieyg/project/honeybee/joint/IDtrans/ID.list")
TF_list_info$gene_symbol<- ID_trans[match(TF_list_info$V6,ID_trans$V1),2]
TF_list<- unique(TF_list_info$gene_symbol)[unique(TF_list_info$gene_symbol)%in% rownames(ORN)]
TF_specific_in_cluster<- gene_in_cluster[gene_in_cluster%in% TF_list];

# promoter scan results 


