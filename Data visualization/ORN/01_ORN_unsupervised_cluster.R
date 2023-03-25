#Step1: ORN cluster by OR gene 
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype_Neuron_subcluster.rds")
DefaultAssay(honeybee) <- "raw_RNA"
#Step1: extract Orco highly exp cluster#####
Idents(honeybee)<- honeybee$Annotation_subcluster
ORN <- subset(honeybee,idents=as.character("Orco+Neuron"))
ORN.list <- SplitObject(ORN, split.by = "orig.ident");
for (i in 1:length(ORN.list)) {
  DefaultAssay(ORN.list[[i]]) <- "raw_RNA"
  ORN.list[[i]] <- SCTransform(ORN.list[[i]], verbose = FALSE,return.only.var.genes = FALSE,variable.features.n=5000,verbose = FALSE)
}
ORN.features <- SelectIntegrationFeatures(object.list = ORN.list, nfeatures = 3000);
ORN.list <- PrepSCTIntegration(object.list = ORN.list, anchor.features = ORN.features)
#integrate RNA using rpca
Neuron_list <- lapply(
  X = ORN.list,
  FUN = ScaleData,
  features = ORN.features,
  verbose = FALSE
)
Neuron_list <- lapply(
  X = Neuron_list,
  FUN = RunPCA,
  features = ORN.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list = Neuron_list,
  normalization.method = "SCT",
  anchor.features = ORN.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)
ORN <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA_ORN",
  dims = 1:30
)

#run LSI on new seurat object with integrated RNA assay
DefaultAssay(ORN) <- "ATAC"
ORN <- RunTFIDF(ORN)
ORN <- FindTopFeatures(ORN, min.cutoff = "q25")
ORN <- RunSVD(ORN,n=25)
#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
ORN_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI_ORN",
  reductions = ORN@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
ORN@reductions$integratedLSI_ORN <- ORN_atac@reductions$integratedLSI_ORN
#####done integrate ATAC and RNA ################
# RNA analysis
library(dplyr)
#first cluster for ORN
# vst
DefaultAssay(ORN) <- "integratedRNA_ORN"
ORN <- FindVariableFeatures(ORN, selection.method = "vst",features=500)
top200 <- head(VariableFeatures(ORN),200)
hvf.info <- HVFInfo(object = ORN,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/01_first_cluster/Find_var_RNA.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
DefaultAssay(ORN) <- "integratedRNA_ORN"
ORN <- ScaleData(ORN,features=rownames(ORN))
ORN <- RunPCA(ORN,features=top200)
ORN <- RunUMAP(ORN,dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
ORN <- RunTSNE(ORN,dims = 1:50, min.dist = 0.001,check_duplicates = FALSE, reduction.name = 'tsne.rna', reduction.key = 'rnaTSNE_')
ORN <- FindNeighbors(object = ORN, reduction = 'pca', dims = 1:50)
ORN <- FindClusters(object = ORN, verbose = FALSE,  resolution =4, algorithm = 3)
table(ORN$seurat_clusters)

Idents(ORN)<-ORN$orig.ident
ORN$orig.ident<-factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )

pdf("./05_ORN_cluster/01_first_cluster/raw_ORN_cluster.pdf",width=6,height=5)
DimPlot(ORN, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
DimPlot(ORN, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "umap.rna",group.by = "seurat_clusters")
#DimPlot(ORN, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "wnn.umap",group.by = "seurat_clusters")
dev.off()

# plot the dotplot by all receptor gene 
DefaultAssay(ORN) <- "raw_RNA"
Idents(ORN)<-ORN$seurat_clusters
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
pdf("./05_ORN_cluster/01_first_cluster/first_raw_OR_dotplot_rawRNA-New-filtration-standard.pdf",width=30, height=12)
p1<-DotPlot(ORN, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p1;
dev.off()

DefaultAssay(ORN) <- "RNA"
Idents(ORN)<-ORN$seurat_clusters
pdf('./05_ORN_cluster/01_first_cluster/first_coreceptor_VlnPlot-New-filtration-standard.pdf',width=15, height=10)
print( VlnPlot(ORN, features = Orco, ncol = 1, pt.size = 0) )
dev.off()

saveRDS(ORN,"./05_ORN_cluster/01_first_cluster/ORN_integrated_antenna_first_cluster.rds")

#Step2: remove cluster without Orco 
ORN <- readRDS("./05_ORN_cluster/01_first_cluster/ORN_integrated_antenna_first_cluster.rds")
all_cluster<-levels(ORN$seurat_clusters)
rm_cluster<-c("44")
onecluster <- subset(ORN,idents=setdiff(all_cluster,rm_cluster))
onecluster.list <- SplitObject(onecluster, split.by = "orig.ident")
for (i in 1:length(onecluster.list)) {
  DefaultAssay(onecluster.list[[i]]) <- "raw_RNA"
  onecluster.list[[i]] <- SCTransform(onecluster.list[[i]],return.only.var.genes = FALSE,variable.features.n=5000,verbose = FALSE)
}
onecluster.features <- SelectIntegrationFeatures(object.list = onecluster.list, nfeatures = 2000)
onecluster.list <- PrepSCTIntegration(object.list = onecluster.list, anchor.features = onecluster.features)

#integrate RNA using rpca
onecluster.list <- lapply(
  X = onecluster.list,
  FUN = ScaleData,
  features = onecluster.features,
  verbose = FALSE
)
Neuron_list <- lapply(
  X = onecluster.list,
  FUN = RunPCA,
  features = onecluster.features,
  verbose = FALSE
)
integration_anchors <- FindIntegrationAnchors(
  object.list = Neuron_list,
  normalization.method = "SCT",
  anchor.features = onecluster.features,
  dims = 1:30,
  reduction = "rpca",
  k.anchor = 20,
)

all_features <- lapply(Neuron_list, row.names) %>% Reduce(intersect,.) 
onecluster <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA_onecluster",
  features.to.integrate = all_features,
  dims = 1:30
)

#run LSI on new seurat object with integrated RNA assay
DefaultAssay(onecluster) <- "ATAC"
onecluster <- RunTFIDF(onecluster)
onecluster <- FindTopFeatures(onecluster, min.cutoff = "q15")
onecluster <- RunSVD(onecluster)
#integrate embeddings and output new object to prevent overwriting integrated RNA
#I can't quite remember why I did this the first time, but it was the work around that I needed
onecluster_atac <- IntegrateEmbeddings(
  anchorset = integration_anchors,
  new.reduction.name = "integratedLSI_onecluster",
  reductions = onecluster@reductions$lsi
)

#copy integrated LSI from duplicate seurat object to original object
onecluster@reductions$integratedLSI_onecluster <- onecluster_atac@reductions$integratedLSI_onecluster

DefaultAssay(onecluster) <- "integratedRNA_onecluster"
onecluster <- FindVariableFeatures(onecluster, selection.method = "vst",features=500)
top300 <- head(VariableFeatures(onecluster),300)
hvf.info <- HVFInfo(object = onecluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[301:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/02_second_cluster/Find_var_RNA_top300.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()

# RNA analysis
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
onecluster <- ScaleData(onecluster, features =rownames(onecluster),verbose = FALSE)
onecluster <- RunPCA(onecluster,features=top300) 
pdf("./05_ORN_cluster/02_second_cluster/ElbowPlot.pdf")
ElbowPlot(onecluster,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
onecluster <- FindNeighbors(object = onecluster, reduction = 'pca', dims = 1:50)
onecluster <- FindClusters( object = onecluster, verbose = FALSE, resolution =4,algorithm = 3)
table(onecluster$seurat_clusters)
onecluster <- RunTSNE(
  object = onecluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:40
)

pdf("./05_ORN_cluster/02_second_cluster/second_ORN_cluster_WNN.pdf",width=6,height=5)
DimPlot(onecluster, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(onecluster, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
# plot the dotplot by all receptor gene 
DefaultAssay(onecluster) <- "raw_RNA"
Idents(onecluster)<-onecluster$seurat_clusters
pdf("./05_ORN_cluster/02_second_cluster/second_raw_OR_dotplot_rawRNA.pdf",width=25, height=12)
p<-DotPlot(onecluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
pdf('./05_ORN_cluster/02_second_cluster/second_coreceptor_VlnPlot_WNN.pdf',width=15, height=10)
print( VlnPlot(onecluster, features =Orco, ncol = 1, pt.size = 0) )
dev.off()
saveRDS(onecluster,"./05_ORN_cluster/02_second_cluster/ORN_integrated_antenna_withOr2_second.rds")

# Step3: distinguish the OR cluster 
 
dotplot_data<- p$data
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled > 1){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- c("23","28","32","33","44")

# each multi OR cluster 
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster/02_second_cluster/distinguish_multi_OR_for_secondcluster.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(onecluster,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
# add max_exp OR label for each cell
ORN_count<-obj@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
  for (i in 1:length(colnames(ORN_matrix))){
    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
}
  }
barcode_label<-barcode_label[barcode_label$label!="NA",]
barcode_label<-barcode_label[order(barcode_label$label),]
##cell cosine simility heatmap 
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col <- myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
names(col)<-unique(barcode_label_pheatmap$OR)
ann_colors= list(OR = col)
p1<-pheatmap(trans_dist,
       cluster_cols = TRUE,
       cluster_rows = TRUE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
#calculate the cosine simility within group and between groups;
rownames(barcode_label)<-barcode_label$barcode
within_group<-c()
between_group<-c()
for (i in 1:nrow(trans_dist)){
  for (j in 1:ncol(trans_dist)){
    if(i!=j){
      if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
        within_group<-c(within_group,trans_dist[i,j]);
      }
      else{between_group<-c(between_group,trans_dist[i,j])}
    }
  }
}
# calculate the FC 
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame(cluster,log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)

# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
              geom_density(alpha=.25) + theme_classic() 
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#method = "t.test"
# plot OR heatmap 
DefaultAssay(obj)<-"raw_RNA";
obj<-ScaleData(obj,features=all_receptor_gene);
#p4<-dittoHeatmap(obj,obj_features,slot ="scale.data",cluster_cols=T,scaled.to.max = TRUE)

top_right<-plot_grid(p2,p3,labels = c("B","C"),ncol = 1)
#right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
last<-plot_grid(p1$gtable, top_right, labels = c('A', ''),rel_widths = c(1, 1),label_size = 12, ncol = 2)
title <- ggdraw() + 
  draw_label(
    paste("Cluster",cluster,":","log2FC=",log2FC),
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
print(add_title)
}
dev.off()

# Step4: select the cluster to subcluster 
# first subcluster 
# multiple_stop_cluster log2FC < 0.05 pvalue < 0.01
> log2FCdata
   cluster       log2FC        pvalue
1        0  0.361562447  0.000000e+00
2        1  0.057765515  1.342092e-25
3        2  0.010681568  2.994228e-01
4        4  0.054940626  1.506329e-14
5        7  0.027203784  1.778686e-03
6        8  0.006426137  1.610359e-04
7        9  0.086226159  9.504924e-17
8       10  0.053994417  2.338018e-16
9       11  0.064322690  7.200988e-08
10      12  0.094915834  1.221055e-22
11      19  0.045666546  3.854494e-03
12      24 -0.055565720  1.686493e-02
13      27  0.424511061 3.524018e-145
14      29  0.276256828 1.847488e-118
15      30  0.380640442  2.218520e-66
16      37  0.028155310  6.751134e-01



multiple_stop_cluster <- as.character(c(24))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))

first_subcluster <- subset(onecluster,idents=need2subcluster)

# vst
DefaultAssay(first_subcluster) <- "RNA"
first_subcluster <- FindVariableFeatures(first_subcluster, selection.method = "vst",features=500)
top300 <- head(VariableFeatures(first_subcluster),300)
hvf.info <- HVFInfo(object = first_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[301:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/02_second_cluster/01_first_subcluster/Find_var_RNA_top300.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
# RNA analysis
DefaultAssay(first_subcluster) <- "integratedRNA_onecluster"
first_subcluster <- RunPCA(first_subcluster,features=top300) 
#build a tSNE visualization
first_subcluster <- FindNeighbors(object = first_subcluster, reduction = 'pca', dims = 1:50)
first_subcluster <- RunTSNE(
  object = first_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)
first_subcluster <- FindClusters( object = first_subcluster, verbose = FALSE, resolution =3,algorithm = 3)
table(first_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/02_second_cluster/01_first_subcluster/second_ORN_cluster_WNN_first_subcluster.pdf",width=6,height=5)
DimPlot(first_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(first_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()

#ORfeatures<-dotplot_data[dotplot_data$id%in%as.character(multiple_cluster),]$features.plot #dotplot_data the second cluster results
DefaultAssay(first_subcluster) <- "raw_RNA"
Idents(first_subcluster)<-first_subcluster$seurat_clusters
pdf("./05_ORN_cluster/02_second_cluster/01_first_subcluster/second_raw_OR_dotplot_rawRNA_first_subcluster.pdf",width=25, height=12)
p<-DotPlot(first_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(first_subcluster,"./05_ORN_cluster/02_second_cluster/01_first_subcluster/second_multiple_classes_first_subcluster.rds");

# Step4: select the cluster to subcluster 
# first_subcluster distinguish OR pipeline 

p<-DotPlot(first_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled > 1){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- c("27")

# distinguish_multi_OR

log2FCdata<-data.frame()
pdf("./05_ORN_cluster/02_second_cluster/01_first_subcluster/first_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(first_subcluster,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
# add max_exp OR label for each cell
ORN_count<-obj@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
  for (i in 1:length(colnames(ORN_matrix))){
    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
}
  }

barcode_label<-barcode_label[barcode_label$label!="NA",]
barcode_label<-barcode_label[order(barcode_label$label),]

##cell cosine simility heatmap 
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col<-myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
names(col)<-unique(barcode_label_pheatmap$OR)
ann_colors= list(OR = col)
p1<-pheatmap(trans_dist,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
#calculate the cosine simility within group and between groups;
rownames(barcode_label)<-barcode_label$barcode
within_group<-c()
between_group<-c()
for (i in 1:nrow(trans_dist)){
  for (j in 1:ncol(trans_dist)){
    if(i!=j){
      if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
        within_group<-c(within_group,trans_dist[i,j]);
      }
      else{between_group<-c(between_group,trans_dist[i,j])}
    }
  }
}
# calculate the FC 
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame(cluster,log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)
# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
              geom_density(alpha=.25) + theme_classic() 
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#method = "t.test"
# plot OR heatmap 
DefaultAssay(obj)<-"raw_RNA";
obj<-ScaleData(obj,features=all_receptor_gene);
#p4<-dittoHeatmap(obj,obj_features,slot ="scale.data",cluster_cols=T,scaled.to.max = TRUE)
top_right<-plot_grid(p2,p3,labels = c("B","C"),ncol = 1)
#right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
last<-plot_grid(p1$gtable, top_right, labels = c('A', ''),rel_widths = c(1, 1),label_size = 12, ncol = 2)
title <- ggdraw() + 
  draw_label(
    paste("Cluster",cluster,":","log2FC=",log2FC),
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
print(add_title)
}
dev.off()

> log2FCdata
   cluster       log2FC        pvalue
1        0  0.040918141  5.317597e-07
2        2  0.171947152 5.438821e-184
3        5  0.073465131  2.569478e-29
4        7  0.014313027  1.767196e-01
5        9  0.064773819  6.890063e-19
6       11  0.060175748  1.499502e-07
7       12  0.034070277  7.066169e-07
8       13  0.273728525 1.648461e-258
9       15  0.004260757  4.954401e-02
10      19 -0.012378076  1.928976e-01
11      20  0.383284281 2.443784e-147
12      22  0.243415498 8.608567e-104
13      25  0.219771920  2.244714e-24
14      28  0.364610027 2.459936e-121
15      33 -0.059067642  1.823119e-01


#second recluster 

multiple_stop_cluster <- as.character(c(19,20,22,25,33))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))
second_subcluster <- subset(first_subcluster,idents=need2subcluster)

# vst
DefaultAssay(second_subcluster) <- "RNA"
second_subcluster <- FindVariableFeatures(second_subcluster, selection.method = "vst",features=500)
top200 <- head(VariableFeatures(second_subcluster),200)
hvf.info <- HVFInfo(object = second_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/02_second_cluster/02_second_subcluster/Find_var_RNA_top200.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
# RNA analysis
DefaultAssay(second_subcluster) <- "integratedRNA_onecluster"
second_subcluster <- RunPCA(second_subcluster,features=top200) 
#build a tSNE visualization
second_subcluster <- FindNeighbors(object = second_subcluster, reduction = 'pca', dims = 1:50)
second_subcluster <- RunTSNE(
  object = second_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)
second_subcluster <- FindClusters( object = second_subcluster, verbose = FALSE, resolution =3,algorithm = 3)
table(second_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/02_second_cluster/02_second_subcluster/second_ORN_cluster_WNN_second_subcluster.pdf",width=6,height=5)
DimPlot(second_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(second_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()

#ORfeatures<-dotplot_data[dotplot_data$id%in%as.character(multiple_cluster),]$features.plot #dotplot_data the second cluster results
DefaultAssay(second_subcluster) <- "raw_RNA"
Idents(second_subcluster)<-second_subcluster$seurat_clusters
pdf("./05_ORN_cluster/02_second_cluster/02_second_subcluster/second_raw_OR_dotplot_rawRNA_second_subcluster.pdf",width=25, height=12)
p<-DotPlot(second_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(second_subcluster,"./05_ORN_cluster/02_second_cluster/02_second_subcluster/second_multiple_classes_second_subcluster.rds");

# Step4: select the cluster to subcluster 
# second_subcluster distinguish OR pipeline 

p<-DotPlot(second_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled > 1){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])

one_classes <- c("34")

# distinguish_multi_OR

log2FCdata<-data.frame()
pdf("./05_ORN_cluster/02_second_cluster/02_second_subcluster/second_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(second_subcluster,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
# add max_exp OR label for each cell
ORN_count<-obj@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
  for (i in 1:length(colnames(ORN_matrix))){
    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
}
  }

barcode_label<-barcode_label[barcode_label$label!="NA",]
barcode_label<-barcode_label[order(barcode_label$label),]

##cell cosine simility heatmap 
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col<-myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
names(col)<-unique(barcode_label_pheatmap$OR)
ann_colors= list(OR = col)
p1<-pheatmap(trans_dist,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
#calculate the cosine simility within group and between groups;
rownames(barcode_label)<-barcode_label$barcode
within_group<-c()
between_group<-c()
for (i in 1:nrow(trans_dist)){
  for (j in 1:ncol(trans_dist)){
    if(i!=j){
      if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
        within_group<-c(within_group,trans_dist[i,j]);
      }
      else{between_group<-c(between_group,trans_dist[i,j])}
    }
  }
}
# calculate the FC 
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame(cluster,log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)
# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
              geom_density(alpha=.25) + theme_classic() 
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#method = "t.test"
# plot OR heatmap 
DefaultAssay(obj)<-"raw_RNA";
obj<-ScaleData(obj,features=all_receptor_gene);
#p4<-dittoHeatmap(obj,obj_features,slot ="scale.data",cluster_cols=T,scaled.to.max = TRUE)
top_right<-plot_grid(p2,p3,labels = c("B","C"),ncol = 1)
#right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
last<-plot_grid(p1$gtable, top_right, labels = c('A', ''),rel_widths = c(1, 1),label_size = 12, ncol = 2)
title <- ggdraw() + 
  draw_label(
    paste("Cluster",cluster,":","log2FC=",log2FC),
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
print(add_title)
}
dev.off()

> log2FCdata
   cluster       log2FC        pvalue
1        0 -0.002756693  8.130950e-01
2        2 -0.014691058  1.242722e-01
3        3  0.036546437  9.384850e-04
4        4  0.200763810 8.290952e-221
5        8  0.074571865  1.802188e-14
6       11  0.122188125  4.676465e-38
7       15 -0.054723321  4.193709e-05
8       16  0.282845296 2.494778e-278
9       17  0.073773401  1.215635e-07
10      18  0.025448196  4.641565e-04


multiple_stop_cluster <- as.character(c(15,17))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))
third_subcluster <- subset(second_subcluster,idents=need2subcluster)

# vst
DefaultAssay(third_subcluster) <- "RNA"
third_subcluster <- FindVariableFeatures(third_subcluster, selection.method = "vst",features=500)
top200 <- head(VariableFeatures(third_subcluster),300)
hvf.info <- HVFInfo(object = third_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/02_second_cluster/03_third_subcluster/Find_var_RNA_top200.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
# RNA analysis
DefaultAssay(third_subcluster) <- "integratedRNA_onecluster"
third_subcluster <- RunPCA(third_subcluster,features=top200) 
#build a tSNE visualization
third_subcluster <- FindNeighbors(object = third_subcluster, reduction = 'pca', dims = 1:50)
third_subcluster <- RunTSNE(
  object = third_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)
third_subcluster <- FindClusters( object = third_subcluster, verbose = FALSE, resolution =3,algorithm = 3)
table(third_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/02_second_cluster/03_third_subcluster/second_ORN_cluster_WNN_third_subcluster.pdf",width=6,height=5)
DimPlot(third_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(third_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.01,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()

#ORfeatures<-dotplot_data[dotplot_data$id%in%as.character(multiple_cluster),]$features.plot #dotplot_data the second cluster results
DefaultAssay(third_subcluster) <- "raw_RNA"
Idents(third_subcluster)<-third_subcluster$seurat_clusters
pdf("./05_ORN_cluster/02_second_cluster/03_third_subcluster/second_raw_OR_dotplot_rawRNA_third_subcluster.pdf",width=25, height=12)
p<-DotPlot(third_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(third_subcluster,"./05_ORN_cluster/02_second_cluster/03_third_subcluster/second_multiple_classes_third_subcluster.rds");

# Step4: select the cluster to subcluster 
# third_subcluster distinguish OR pipeline 

p<-DotPlot(third_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled > 1){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])

one_classes <- c("19","22","23","30","31")

# distinguish_multi_OR

log2FCdata<-data.frame()
pdf("./05_ORN_cluster/02_second_cluster/03_third_subcluster/third_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(third_subcluster,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
# add max_exp OR label for each cell
ORN_count<-obj@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
  for (i in 1:length(colnames(ORN_matrix))){
    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
}
  }

barcode_label<-barcode_label[barcode_label$label!="NA",]
barcode_label<-barcode_label[order(barcode_label$label),]

##cell cosine simility heatmap 
embeddings <- Embeddings(object = obj, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col<-myUmapcolors[1:length(unique(barcode_label_pheatmap$OR))]
names(col)<-unique(barcode_label_pheatmap$OR)
ann_colors= list(OR = col)
p1<-pheatmap(trans_dist,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
#calculate the cosine simility within group and between groups;
rownames(barcode_label)<-barcode_label$barcode
within_group<-c()
between_group<-c()
for (i in 1:nrow(trans_dist)){
  for (j in 1:ncol(trans_dist)){
    if(i!=j){
      if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
        within_group<-c(within_group,trans_dist[i,j]);
      }
      else{between_group<-c(between_group,trans_dist[i,j])}
    }
  }
}
# calculate the FC 
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame(cluster,log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)
# plot density line 
# manage data
type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
var<-c(within_group,between_group)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("within-OR","between-OR"))
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
              geom_density(alpha=.25) + theme_classic() 
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#method = "t.test"
# plot OR heatmap 
DefaultAssay(obj)<-"raw_RNA";
obj<-ScaleData(obj,features=all_receptor_gene);
#p4<-dittoHeatmap(obj,obj_features,slot ="scale.data",cluster_cols=T,scaled.to.max = TRUE)
top_right<-plot_grid(p2,p3,labels = c("B","C"),ncol = 1)
#right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" ","D"))
last<-plot_grid(p1$gtable, top_right, labels = c('A', ''),rel_widths = c(1, 1),label_size = 12, ncol = 2)
title <- ggdraw() + 
  draw_label(
    paste("Cluster",cluster,":","log2FC=",log2FC),
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
print(add_title)
}
dev.off()

> log2FCdata
