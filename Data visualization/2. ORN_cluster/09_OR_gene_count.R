library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)
set.seed(1234);
library(pheatmap)
# remove nopower cluster,then plot UMAP
ORN<- readRDS("./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_order_by_tree_recall_peak.rds")
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )
# the OR gene exp order plot 
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))

ORN_count<-ORN@assays$raw_RNA@counts
ORN_count<-ORN_count[,which(colnames(ORN_count)%in%colnames(ORN))]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
OR_gene<- all_receptor_gene[which(all_receptor_gene%in%rownames(ORN_matrix))]
ORN_matrix<-ORN_matrix[OR_gene,]

OR_count_data<- data.frame(OR=rownames(ORN_matrix),count=rowSums(ORN_matrix))
OR_count_data<- OR_count_data[order(OR_count_data$count,decreasing=T),]
OR_count_data$OR<- factor(OR_count_data$OR,levels=OR_count_data$OR)

pdf("./05_ORN_cluster/03_coexp_cluster/OR_total_count_order.pdf",width=20,height=4)
ggplot(OR_count_data[-1,],aes(x=OR,y=log10(count)))+geom_point()+theme_light()+
theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
dev.off()


