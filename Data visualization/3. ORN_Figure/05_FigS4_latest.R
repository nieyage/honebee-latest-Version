library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));


obj<- readRDS("./05_ORN_cluster2/05_combination_group_recluster/obj_recluster.rds")
Idents(obj)<-obj$Annotation
identityMapping <- c('C1' = 'C39', 'C2' = 'C40','C3' = 'C41', 'C4' = 'C42')
obj <- RenameIdents(obj, identityMapping)
obj_features<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
Idents(obj)<-obj$Annotation

# Find the DEG among the 4 cluster:
DefaultAssay(obj)<- "raw_RNA"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
write.csv(markers,"./00_Figure/FigS4/FigS4-DEG_markers-DEG_heatmap_need2select_showmarkers.csv")

solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC);
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[6:8])

pdf("./00_Figure/FigS4/FigS4-DEG_markers-DEG_heatmap.pdf",width=5,height=5)
DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

# show select gene exp and ATAC signal
gene<- unique(markers$gene)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  ######Visulize track and RNA exp#####
  idents.plot <- Idents(obj)
  # plot region 
pdf("./00_Figure/FigS4/FigS4-DEG-track.pdf",width=10,height=5)
for(i in gene){
    p1<-CoveragePlot(
    object = obj,
    region = i,
    features=i,
    window = 150,
    expression.assay = "raw_RNA",
    expression.slot = "data",
    extend.upstream = 500,
    annotation = TRUE,
    extend.downstream = 500
  )
   p2<-  p1 + scale_fill_manual(color_for_cluster)
print(p2)
}


dev.off()
library(dittoSeq)
pdf("./00_Figure/FigS4/FigS4_DEG_label_heatmap.pdf",width=5,height=5)

dittoHeatmap(obj, genes = top10$gene,
             annot.by = c("Annotation"),
             scaled.to.max = F,
             treeheight_row = 10,
             heatmap.colors = colorRampPalette(c('#1A5592','white',"#B83D3D"))(50),
             show_rownames=F,
             highlight.features = gene)
dev.off()

# label DEG in heatmap 
library(ComplexHeatmap)
DefaultAssay(obj)<-"raw_RNA"
mat<-GetAssayData(obj,slot='scale.data')
cluster_info<- sort(obj$Annotation)
mat<- as.matrix(mat[unique(top10$gene),names(cluster_info)])
gene<- unique(gene)
gene_pos<- which(rownames(mat)%in% gene)
row_anno<- rowAnnotation(gene=anno_mark(at=gene_pos,label=gene))
pdf("./00_Figure/FigS4/FigS4_DEG_label_heatmap.pdf",width=5,height=5)
Heatmap(mat,
  cluster_rows=F,
  cluster_columns=F,
  show_column_names=F,
  show_row_names =F,
  column_split=cluster_info,
  right_annotation=row_anno
  )
dev.off()


# Find the DEP among the 4 cluster:
DefaultAssay(obj)<- "peaks_ORN_subcluster"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
write.csv(markers,"./00_Figure/FigS4/FigS4-DEG_markers-DEP_heatmap.csv")

sambaNight<- c("#1873CC", "#1798E5" ,"#00BFFF", "#4AC596" ,"#00CC00" ,"#A2E700", "#FFFF00" ,"#FFD200","#FFA500")
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

obj<- ScaleData(obj)


pdf("./00_Figure/FigS4/FigS4-DEP_markers-DEP_heatmap.pdf",width=5,height=5)
DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1,disp.max = 1,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
dev.off();


obj_features<- markers$gene
barcode_label<-data.frame(barcode=colnames(obj),label=Idents(obj))
DefaultAssay(obj)<-"raw_RNA"
obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
barcode_label<- barcode_label[order(barcode_label$label),]
obj_data<- obj_data[rownames(barcode_label),]
C39_data<-obj_data[rownames(barcode_label[barcode_label$label=="C1",]),]
C40_data<-obj_data[rownames(barcode_label[barcode_label$label=="C2",]),]
C41_data<-obj_data[rownames(barcode_label[barcode_label$label=="C3",]),]
C42_data<-obj_data[rownames(barcode_label[barcode_label$label=="C4",]),]

library(pheatmap)

clusterMatrix <- function(input_matrix) {
  # Define the clustering method and other parameters
  clustering_method <- "complete"  # You can change this to other methods like "ward.D", "single", etc.
  # Perform clustering
  p <- pheatmap(
    input_matrix,
    clustering_method = clustering_method,
    cluster_cols = F,
    cluster_rows = T,
  )
  clustered_matrix <- input_matrix[p$tree_row$order,]
  # Return the clustered matrix
  return(clustered_matrix)
}
smooth_column <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 2:(length(col) - 1)) {
    smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
  }
  smoothed[1] <- (col[1] + col[2]) / 2
  smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
  return(smoothed)
}

C39_data_clustered<- clusterMatrix(C39_data)
C39_data_clustered_smoothed_data <- as.data.frame(lapply(C39_data_clustered, smooth_column))
rownames(C39_data_clustered_smoothed_data)<- rownames(C39_data_clustered)

C40_data_clustered<- clusterMatrix(C40_data)
C40_data_clustered_smoothed_data <- as.data.frame(lapply(C40_data_clustered, smooth_column))
rownames(C40_data_clustered_smoothed_data)<- rownames(C40_data_clustered)
C41_data_clustered<- clusterMatrix(C41_data)
C41_data_clustered_smoothed_data <- as.data.frame(lapply(C41_data_clustered, smooth_column))
rownames(C41_data_clustered_smoothed_data)<- rownames(C41_data_clustered)
C42_data_clustered<- clusterMatrix(C42_data)
C42_data_clustered_smoothed_data <- as.data.frame(lapply(C42_data_clustered, smooth_column))
rownames(C42_data_clustered_smoothed_data)<- rownames(C42_data_clustered)


clustered_smoothed_data<- rbind(C39_data_clustered_smoothed_data,
	C40_data_clustered_smoothed_data,
	C41_data_clustered_smoothed_data,
	C42_data_clustered_smoothed_data)

# 对每一列进行平滑
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[6:8])
barcode_label_pheatmap<-data.frame(label=c(rep("C39",nrow(C39_data)),rep("C40",nrow(C40_data)),rep("C41",nrow(C41_data)),rep("C42",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-rownames(clustered_smoothed_data)
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)

smoothed_data <- as.data.frame(lapply(obj_data, smooth_column))
rownames(smoothed_data)<- rownames(obj_data)
barcode_label_pheatmap<-data.frame(label=c(rep("C39",nrow(C39_data)),rep("C40",nrow(C40_data)),rep("C41",nrow(C41_data)),rep("C42",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-rownames(smoothed_data)

clustered_smoothed_data<-  t(scale(t(clustered_smoothed_data),scale=T))

pdf("./00_Figure/FigS4/FigS4_DEG_smoothed_heatmap.pdf",height=8,width=8)
pheatmap(t(clustered_smoothed_data),
             cluster_cols = F,
             cluster_rows = TRUE,
             #color = colorRampPalette(c("white", "#CC0000"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()





