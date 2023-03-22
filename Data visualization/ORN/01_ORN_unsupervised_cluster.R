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
  ORN.list[[i]] <- SCTransform(ORN.list[[i]], verbose = FALSE)
}
for (i in seq_len(length(ORN.list))) {
  DefaultAssay(ORN.list[[i]]) <- "SCT"
}
ORN.features <- SelectIntegrationFeatures(object.list = ORN.list, nfeatures = 3000);
ORN.list <- PrepSCTIntegration(object.list = ORN.list, anchor.features = ORN.features)
#integrate RNA using rpca
Neuron_list <- lapply(
  X = ORN.list,
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
DefaultAssay(ORN) <- "RNA"
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
IR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="IR",]$gene_name)
all_receptor_gene <- c(OR_gene,IR_gene,GR_gene)
pdf("./05_ORN_cluster/01_first_cluster/first_raw_OR_dotplot_rawRNA-New-filtration-standard.pdf",width=30, height=12)
p1<-DotPlot(ORN, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p1;
dev.off()

DefaultAssay(ORN) <- "RNA"
Idents(ORN)<-ORN$seurat_clusters
pdf('./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/first_coreceptor_VlnPlot-New-filtration-standard.pdf',width=15, height=10)
print( VlnPlot(ORN, features = c("Or2","LOC100577715","LOC100576282","LOC107963999","LOC102656838"), ncol = 1, pt.size = 0) )
dev.off()

saveRDS(ORN,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/ORN_integrated_antenna_gtf_first.rds")



#Step2: remove cluster without OR2 
all_cluster<-levels(ORN$seurat_clusters)
rm_cluster<-c("31")
onecluster <- subset(ORN,idents=setdiff(all_cluster,rm_cluster))
onecluster.list <- SplitObject(onecluster, split.by = "orig.ident")
for (i in 1:length(onecluster.list)) {
  DefaultAssay(onecluster.list[[i]]) <- "raw_RNA"
  onecluster.list[[i]] <- SCTransform(onecluster.list[[i]], verbose = FALSE)
}
for (i in seq_len(length(onecluster.list))) {
  DefaultAssay(onecluster.list[[i]]) <- "SCT"
}
onecluster.features <- SelectIntegrationFeatures(object.list = onecluster.list, nfeatures = 3000)
onecluster.list <- PrepSCTIntegration(object.list = onecluster.list, anchor.features = onecluster.features)

#integrate RNA using rpca

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
onecluster <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT",
  new.assay.name = "integratedRNA_onecluster",
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

# RNA analysis
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
onecluster <- RunPCA(onecluster) 
#build a tSNE visualization
onecluster <- FindNeighbors(object = onecluster, reduction = 'pca', dims = 1:50)
onecluster <- FindClusters( object = onecluster, verbose = FALSE, resolution =3,algorithm = 3)
onecluster <- RunTSNE(
  object = onecluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)
table(onecluster$seurat_clusters)

pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_ORN_cluster_WNN.pdf",width=6,height=5)
DimPlot(onecluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(onecluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
# plot the dotplot by all receptor gene 
DefaultAssay(onecluster) <- "raw_RNA"
Idents(onecluster)<-onecluster$seurat_clusters
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_raw_OR_dotplot_rawRNA.pdf",width=25, height=12)
p<-DotPlot(onecluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
pdf('./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_coreceptor_VlnPlot_WNN.pdf',width=15, height=10)
print( VlnPlot(onecluster, features = c("Or2","LOC100577715","LOC100576282","LOC107963999","LOC102656838"), ncol = 1, pt.size = 0) )
dev.off()
saveRDS(onecluster,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/ORN_integrated_antenna_withOr2_second.rds")

# Step3: distinguish the OR cluster 
DefaultAssay(onecluster) <- "raw_RNA"
Idents(onecluster) <- "seurat_clusters";
all_cluster<-levels(onecluster$seurat_clusters)

dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled>1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled> 0){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>10&&dotplot_data[i,]$avg.exp.scaled>2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),];
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot=="Or2"),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes<-cluster_info[cluster_info$Freq>1,1]
one_classes<-cluster_info[cluster_info$Freq==1,1]
zero_classes<-cluster_info[cluster_info$Freq==0,1]

# each multi OR cluster 
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/distinguish_multi_OR_for_secondcluster.pdf",width=14,height=6)
for (cluster in multiple_classes_order_id){
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

####OR log2FC
## Random select two ORs to calculate the withingroup and between group FC;
features<-unique(multiple_classes_order_feature)
# add max exp OR label
ORN_count<-onecluster@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
for (i in 1:length(colnames(ORN_matrix))){
    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
}
  }
barcode_label<-barcode_label[barcode_label$label!="NA",]
barcode_label<-barcode_label[order(barcode_label$label),]

#transcriptome distance 
embeddings <- Embeddings(object = onecluster, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]

library(lsa)
trans_dist <- 1-cosine(t(embeddings))
rownames(barcode_label)<-barcode_label$barcode
features<-unique(multiple_classes_order_feature)
data<-data.frame()
# select features pairs 
remaining_gene<-features;
for (gene1 in features){
  print(gene1)
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    print(gene2)
    barcode<-barcode_label[which(barcode_label$label%in% c(gene1,gene2)),1]
    barcode_label_subset<-barcode_label[which(barcode_label$label%in% c(gene1,gene2)),]
    # extract the distance of within and between groups 
    within_group<-c()
    between_group<-c()
    trans_dist_subset<-trans_dist[barcode,barcode];
    for (i in 1:nrow(trans_dist_subset)){
         for (j in 1:ncol(trans_dist_subset)){
           if(i!=j){
             if(barcode_label_subset[rownames(trans_dist_subset)[i],2]==barcode_label_subset[colnames(trans_dist_subset)[j],2]){
               within_group<-c(within_group,trans_dist_subset[i,j]);
             }
             else{between_group<-c(between_group,trans_dist_subset[i,j])}
           }
         }
       }
    # calculate the FC 
    if(!is.null(median(within_group))){
    if(!is.null(median(between_group))){
      log2FC<-log2(median(between_group))-log2(median(within_group));
      test<-wilcox.test(within_group,between_group);
      pvalue<-test$p.value;
      data_subset<-data.frame(gene1,gene2,log2FC,pvalue)
      data<-rbind(data,data_subset)
    }}
    
  }
}

# plot the log2FC distribution 
# plot density line 
# manage data
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/select_log2FC_cutoff_gene_second.pdf",width=10,height=5)
ggplot(data, aes(x=log2FC)) + xlab("log2FC")+
              geom_density(alpha=.25) + theme_classic() 
d <- density(data$log2FC)
d$x[which.min(abs(diff(d$y)))]
hist(data$log2FC,prob=TRUE)
lines(d, col="red", lty=2)
df<-data
km <- kmeans(df$log2FC,centers=5)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=log2FC)) + geom_histogram(aes(fill=clust,y=..count../sum(..count..)),binwidth=0.5, color="grey50")+stat_density(geom="line", color="red")
dev.off()

## define log2FC < 0.5 as cutoff
# For gene FC

last_data<-data.frame()
for(cluster in multiple_classes_order_id){
multi_data<-dotplot_data[dotplot_data$id %in%cluster,]
cluster_features<-multi_data$features.plot
tmp_data<-data[data$gene1%in% cluster_features,]
tmp_data<-tmp_data[tmp_data$gene2%in% cluster_features,]
for(i in 1:nrow(tmp_data)){
  gene1<-tmp_data$gene1[i];
  gene2<-tmp_data$gene2[i];
  if(dotplot_data[which(dotplot_data$features.plot==gene1),5]==dotplot_data[which(dotplot_data$features.plot==gene2),5]){
    data_tmp<-tmp_data[i,];
    data_tmp$cluster<-cluster;
    last_data<-rbind(last_data,data_tmp);
  }
}
}
write.csv(last_data,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/multiOR_pair_log2FC_second.csv")

# Step4: select the cluster to subcluster 
# first subcluster 

multiple_cluster <- c(0,4,7,12,14,16,18,19,20,21,34)
zero_classes
first_subcluster <- subset(onecluster,idents=as.character(multiple_cluster,zero_classes))
# RNA analysis
DefaultAssay(first_subcluster) <- "integratedRNA_onecluster"
VF<-Var
first_subcluster <- RunPCA(first_subcluster) 
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

pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_ORN_cluster_WNN_first_subcluster.pdf",width=6,height=5)
DimPlot(first_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(first_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()

ORfeatures<-dotplot_data[dotplot_data$id%in%as.character(multiple_cluster),]$features.plot #dotplot_data the second cluster results
DefaultAssay(first_subcluster) <- "raw_RNA"
Idents(first_subcluster)<-first_subcluster$seurat_clusters
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_raw_OR_dotplot_rawRNA_first_subcluster.pdf",width=25, height=12)
p<-DotPlot(first_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(first_subcluster,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_multiple_classes_first_subcluster.rds");

# Step4: select the cluster to subcluster 
# first_subcluster distinguish OR pipeline 

p<-DotPlot(first_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled>1.5){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled> 0){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>10&&dotplot_data[i,]$avg.exp.scaled>2.4){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot=="Or2"),];

cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes<-cluster_info[cluster_info$Freq>1,1]
one_classes<-cluster_info[cluster_info$Freq==1,1]
zero_classes<-cluster_info[cluster_info$Freq==0,1]

# distinguish_multi_OR
# each multi OR cluster 
library(pheatmap)
library(dittoSeq)
library(cowplot)
#ORN <- RunPCA(ORN) 
log2FCdata<-data.frame()
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/first_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes_order_id){
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

#second recluster 
second_subcluster <- subset(first_subcluster,idents=c("0_2","0_1","4_3","16_2","19_2","21"))
DefaultAssay(second_subcluster) <- "integratedRNA_onecluster";

second_subcluster <- RunPCA(second_subcluster,features=input_feature) 
#build a tSNE visualization
second_subcluster <- FindNeighbors(object = second_subcluster, reduction = 'pca', dims = 1:50)
second_subcluster <- FindClusters( object = second_subcluster, verbose = FALSE, resolution =2,algorithm = 3)
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
table(second_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_ORN_cluster_WNN_second_subcluster.pdf",width=6,height=5)
DimPlot(second_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(second_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()

ORfeatures<-dotplot_data[dotplot_data$id%in%as.character(c("0_2","0_1","4_3","16_2","19_2","21")),]$features.plot #dotplot_data the first_subcluster results
DefaultAssay(second_subcluster) <- "raw_RNA"
Idents(second_subcluster)<-second_subcluster$seurat_clusters
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_raw_OR_dotplot_rawRNA_second_subcluster.pdf",width=25, height=12)
p<-DotPlot(second_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

Idents(second_subcluster)<-second_subcluster$seurat_clusters
second_subcluster <- RenameIdents(second_subcluster, 
  '0'='4_3','1'='0_2_2','2'='0_1_1','3'='0_2_1','4'='21_3','5'='21_1',
  '6'='21_2','7'='19_2','8'='16_2','9'='16_2')

sort<-c('0_1_1',"0_2_1" ,"0_2_2" ,"4_3" ,"16_2" , "19_2","21_1", "21_2","21_3")
second_subcluster$second_subcluster<-Idents(second_subcluster)
second_subcluster$second_subcluster<-factor(second_subcluster$second_subcluster,levels=rev(sort))
Idents(second_subcluster)<-second_subcluster$second_subcluster

DefaultAssay(second_subcluster)<-"raw_RNA"
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_Supervised_some_WNN_normal_OR_dotplot_rawRNA_order_second_subcluster.pdf",width=18, height=9)
p<-DotPlot(second_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
saveRDS(second_subcluster,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_multiple_classes_second_subcluster.rds");

# second_subcluster distinguish OR pipeline 
p<-DotPlot(second_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled>1.5){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled> 2){dotplot_data[i,]$state="Yes"};
  #if(dotplot_data[i,]$pct.exp>10&&dotplot_data[i,]$avg.exp.scaled>2.4){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes<-cluster_info[cluster_info$Freq>1,1]
multiple_classes_data<-dotplot_data[which(dotplot_data$id%in% multiple_classes),];
multiple_classes_order_id<-unique(multiple_classes_data$id);
multiple_classes_order_feature<-multiple_classes_data$features.plot;
# distinguish_multi_OR
# each multi OR cluster 
library(pheatmap)
library(dittoSeq)
library(cowplot)
#ORN <- RunPCA(ORN) 
log2FCdata<-data.frame()
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes_order_id){
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

#third recluster 
third_subcluster <- subset(second_subcluster,idents=c("4_3" ,"16_2" , "19_2"))
DefaultAssay(third_subcluster) <- "integratedRNA_onecluster";

third_subcluster <- RunPCA(third_subcluster,features=input_feature) 
#build a tSNE visualization
third_subcluster <- FindNeighbors(object = third_subcluster, reduction = 'pca', dims = 1:50)
third_subcluster <- FindClusters( object = third_subcluster, verbose = FALSE, resolution =1.8,algorithm = 3)
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
table(third_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/third_ORN_cluster_WNN_third_subcluster.pdf",width=6,height=5)
DimPlot(third_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(third_subcluster, cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.5,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()

ORfeatures<-dotplot_data[dotplot_data$id%in%as.character(c("4_3" ,"16_2" , "19_2")),]$features.plot #dotplot_data the first_subcluster results
DefaultAssay(third_subcluster) <- "raw_RNA"
Idents(third_subcluster)<-third_subcluster$seurat_clusters
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/third_raw_OR_dotplot_rawRNA_third_subcluster.pdf",width=25, height=12)
p<-DotPlot(third_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

Idents(third_subcluster)<-third_subcluster$seurat_clusters
third_subcluster <- RenameIdents(third_subcluster, 
  '0'='4_3_1','1'='4_3_2','2'='19_2','3'='16_2','4'='None')
third_subcluster$third_subcluster<-Idents(third_subcluster)
#third_subcluster$third_subcluster<-factor(third_subcluster$third_subcluster,levels=rev(sort))
Idents(third_subcluster)<-third_subcluster$third_subcluster

DefaultAssay(third_subcluster)<-"raw_RNA"
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_Supervised_some_WNN_normal_OR_dotplot_rawRNA_order_third_subcluster.pdf",width=18, height=9)
p<-DotPlot(third_subcluster, features = unique(ORfeatures)) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
saveRDS(third_subcluster,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_multiple_classes_third_subcluster.rds");



# the last cluster 

onecluster<-readRDS("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/Supervised_WNN_withpowercluster_ORN_second.rds")
DefaultAssay(onecluster) <- "raw_RNA"
Idents(onecluster) <- "seurat_clusters";
all_cluster<-levels(onecluster$seurat_clusters)
rm_cluster<-c("13","35","42","46","47")
ORN <- subset(onecluster,idents=setdiff(all_cluster,rm_cluster))
# All ORN dotplot map #
barcode=colnames(ORN)
seurat_clusters=as.character(ORN$seurat_clusters)
cell_info<-data.frame(barcode,seurat_clusters)
cell_info$subcluster<-"None";
multiple_cluster <- c(0,4,7,12,14,16,18,19,20,21,34)
recluster<-as.character(multiple_cluster)
all_cluster<-unique(cell_info[,2]);
no_recluster<-setdiff(all_cluster,recluster)
for (i in 1:nrow(cell_info)){
  if(cell_info[i,2]%in%no_recluster){
    cell_info[i,3]<-cell_info[i,2]
  }
}

#first recluster
first_subcluster<-readRDS("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_multiple_classes_first_subcluster.rds")
cell_info_first<-data.frame(colnames(first_subcluster),as.character(first_subcluster$first_subcluster))
cell_info_first<-cell_info_first[-which(cell_info_first[,2]%in%c("0_2","0_1","4_3","16_2","19_2","21")),]

#second recluster
second_subcluster<-readRDS("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_multiple_classes_second_subcluster.rds")
cell_info_second<-data.frame(colnames(second_subcluster),as.character(second_subcluster$second_subcluster))
cell_info_second<-cell_info_second[-which(cell_info_second[,2]%in%c("4_3" ,"16_2" , "19_2")),]

third_subcluster<-readRDS("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/second_multiple_classes_third_subcluster.rds")
cell_info_third<-data.frame(colnames(third_subcluster),as.character(third_subcluster$third_subcluster))

colnames(cell_info_first)<-c("barcode","subcluster")
colnames(cell_info_second)<-c("barcode","subcluster")
colnames(cell_info_third)<-c("barcode","subcluster")

cell_info_recluster<-rbind(cell_info_first,cell_info_second,cell_info_third)

# merge cell information 
subcluster<-cell_info_recluster$subcluster[match(cell_info$barcode,cell_info_recluster$barcode)]
for (i in 1:nrow(cell_info)){
  if(!is.na(subcluster[i])){
    cell_info[i,3]<-subcluster[i]
  }
}
rownames(cell_info)<-cell_info$barcode;
cell_info<-cell_info[colnames(ORN),]
ORN$subcluster<-cell_info$subcluster;
saveRDS(ORN,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/Supervised_ORN_cluster_WNN_add_subcluster.rds")

# tree,order,dotplot,umap 
#make tree for dotplot order
ORN<-readRDS("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/Supervised_ORN_cluster_WNN_add_subcluster.rds")
Idents(ORN)<-ORN$subcluster
DefaultAssay(ORN) <- "integratedRNA_onecluster"
ORN2 <- RunPCA(ORN) 
Idents(ORN2)<-ORN2$subcluster
object<-ORN2;
embeddings <- Embeddings(object = object, reduction = "pca")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);


pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/cluster-ORN-tree-cosine_bytop3000features.pdf",width=6,height=6)
#tree <- groupOTU(data.tree)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/cluster-ORN-tree-cosine_bytop3000features-nocircular.pdf",width=8,height=12)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()
#cluster order by tree

m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

DefaultAssay(ORN)<-"raw_RNA"
receptor<-read.csv("/md01/nieyg/project/Neuron/bulk-RNAseq/blast/receptor.csv")
GR<-receptor[grep("gustatory receptor",receptor$V10),]$gene_name
IR<-receptor[grep("ionotropic receptor",receptor$V10),]$gene_name
OR<-read.csv("/md01/nieyg/project/Neuron/bulk-RNAseq/blast/Or-7tm_6-information-new.csv")
OR<-OR[which(OR$gene_name%in% rownames(ORN)),]
OR<-OR$gene_name;
all_receptor_gene<-unique(c(OR,GR,IR));


Idents(ORN)<-factor(ORN$subcluster,levels=cluster_order)

p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled>1.5){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled> 0){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>10&&dotplot_data[i,]$avg.exp.scaled>2.4){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)))
DefaultAssay(ORN)<-"raw_RNA"
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/Supervised_WNN_dotplot-allfeature_orderbytree.pdf",width=22, height=14)
p<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
# change color bar for p2
p2<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') +  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) +
scale_color_gradientn(colours = c('#008080', '#FF00FF'),  name = 'Average\nexpression', oob = scales::squish)
p2
dev.off()


# plot the transcript heatmap and calculate the within variation in cluster;
DefaultAssay(ORN)<-"raw_RNA"
embeddings <- Embeddings(object = ORN, reduction = "pca")[,1:50]
#embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(cluster=ORN$subcluster)
rownames(barcode_label_pheatmap)<-colnames(ORN)
col<-myUmapcolors[1:length(unique(barcode_label_pheatmap$cluster))]
names(col)<-unique(barcode_label_pheatmap$cluster)
ann_colors= list(cluster = col)
library(pheatmap)
pdf("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/trans_dist_heatmap_for_last_ORN.pdf",width=16,height=15)
pheatmap(trans_dist,
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
dev.off()


barcode_label<-data.frame(barcode=colnames(ORN),label=ORN$subcluster)
rownames(barcode_label)<-barcode_label$barcode
# within group trans_dist variation 
within_group_data<-data.frame()
for (cluster in levels(ORN)){
  print(cluster);
  barcode<-barcode_label[which(barcode_label$label==cluster),]$barcode;
  cluster_trans_dist<-trans_dist[barcode,barcode];
  within_group_sd<-sd(as.numeric(cluster_trans_dist));
  within_group_data_subset<-data.frame(cluster,within_group_sd);
  within_group_data <- rbind(within_group_data,within_group_data_subset)
}

write.csv(within_group_data,"./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/within_group_data_sd.csv")
mean(within_group_data$within_group_sd)
#[1] 0.1742372



ORN<-readRDS("./05_ORN_cluster/03_iterative_clustering_just_filterd_gene/Supervised_ORN_cluster_WNN_add_subcluster.rds")
Idents(ORN)<-ORN$subcluster
DefaultAssay(ORN)<-"raw_RNA"
ORN_tau<-ScaleData(ORN)
pdf("./05_ORN_cluster/01_supervised_raw_ORN/modify_tau098_gene_cluster_specfic_heatmap.pdf",width=15,height=15)
DoHeatmap(object = ORN_tau,disp.min=0,disp.max=0.05,features=input_feature,label=TRUE,size = 3.5) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = ORN_tau,features=input_feature,label=TRUE,size = 3.5) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = ORN_tau,features=tau_gene,label=TRUE,size = 3.5) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()

dev.off();




