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
log2FCdata<-data.frame();
DefaultAssay(ORN)<- "integratedRNA_onecluster"
combination_group<- c("p1:0_1","p1:0_0","p1:4_1","p1:4_0")
obj<-subset(ORN,idents=combination_group);
obj_features<- unique(dotplot_data[dotplot_data$id%in%combination_group,]$features.plot)
#####done integrate ATAC and RNA ################
# RNA analysis
DefaultAssay(obj)<-"RNA"
obj<-FindVariableFeatures(obj, selection.method = "vst")
top20 <- head(VariableFeatures(obj),20)

DefaultAssay(obj)<-"raw_RNA"
obj<-ScaleData(obj,rownames(obj))
DefaultAssay(obj)<-"integratedRNA_onecluster"
obj <- RunPCA(obj,features=obj_features) %>% RunUMAP(dims = 1:3, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(obj) <- "peaks_ORN_subcluster"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q25")
obj <- RunSVD(obj,n=25)

obj <- RunUMAP(obj, reduction = 'integratedLSI_onecluster', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
# We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. 
# We use this graph for UMAP visualization and clustering
#obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "integratedLSI_onecluster"), 
  dims.list = list(1:3, 2:4),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", resolution =1, algorithm = 3, verbose = FALSE)
table(obj$seurat_clusters)
table(obj$seurat_clusters,obj$subcluster)
###reorder the level of sample#####
Idents(obj)<-obj$orig.ident
obj$orig.ident<-factor(obj$orig.ident,levels=c("NE","Nurse","Forager"))
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
pdf("./05_ORN_cluster2/05_combination_group_recluster/combination_group_recluster_WNN.pdf.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
###sample
p1 <- DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "umap.rna", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "umap.atac",group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(obj, cols=my47colors[22:30],pt.size = 1.2, reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))& NoLegend()
dev.off()
pdf("./05_ORN_cluster2/05_combination_group_recluster/combination_group_recluster_trans_exp.pdf.pdf",width=15,height=5)
ORN_count<-obj@assays$SCT
barcode_label<-data.frame(barcode=colnames(obj),label=obj$seurat_clusters)
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-barcode_label[colnames(ORN_matrix),]
barcode_label<-barcode_label[order(barcode_label$label),]
# p1 cell cosine simility heatmap 
embeddings <- Embeddings(object = ORN, reduction = "pca")[,1:50]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(label=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col<-brewer.pal(12,"Set3")[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)
#trans_dist<-trans_dist[rownames(barcode_label_pheatmap),rownames(barcode_label_pheatmap)]
p1<-pheatmap(trans_dist,
         cluster_cols = F,
         cluster_rows = F,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = barcode_label_pheatmap,
         annotation_colors = ann_colors,
         annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F)
# calculate the cosine simility within group and between groups;
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
if(!is.null(median(between_group))){
log2FC<-log2(median(between_group))-log2(median(within_group));
test<-wilcox.test(within_group,between_group);
pvalue<-test$p.value;
data_subset<-data.frame("1",log2FC,pvalue)
log2FCdata<-rbind(log2FCdata,data_subset)
}
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
#raw counts heatmap 
# Heat map of expression  value 
    Idents(obj)<-obj$subcluster
    DefaultAssay(obj)<-"SCT"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
    obj_data<-obj_data[rownames(barcode_label),]
    p4<-pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "red"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
    top_right<-plot_grid(p2,p3,labels = c(" "," "),rel_widths = c(2, 1))
    right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" "," "))
    last<-plot_grid(p1$gtable, right, labels = c(' ', ''), label_size = 12, ncol = 2)
    title <- ggdraw() + 
      draw_label(
        paste("combination","log2FC=",log2FC),
        fontface = 'bold',
        x = 0,
        hjust = 0
      )
    add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
    print(add_title)
dev.off()

# Fig4H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(obj)<-"peaks_ORN_subcluster"
cluster_info<- combination_group
pdf("./05_ORN_cluster2/05_combination_group_recluster/combination_group_recluster_trackplot.pdf",width=15,height=5)
  Idents(obj)<-obj$seurat_clusters
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp######
  idents.plot <- Idents(obj)
  # plot region 
  start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
  ranges.show <- paste(seq,start,end,sep="-")
  col<-brewer.pal(12,"Set3")[1:length(cluster_info)]
  p1<-CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 100,
    annotation = TRUE,
    extend.downstream = 100,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 30,
    links=F
  )
  print(p1)
dev.off()
