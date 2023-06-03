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
obj <- RunPCA(obj,features=c(top20,obj_features),reduction.name="obj_features_pca") 
obj <- RunUMAP(obj, dims = 1:20, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(obj) <- "peaks_ORN_subcluster"
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
obj_peaks<- c(rownames(obj)[grep("Group15-695....-.......",rownames(obj))],rownames(obj)[grep("Group15-697....-.......",rownames(obj))],rownames(obj)[grep("Group15-697....-.......",rownames(obj))])
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q25")
#obj <- RunSVD(obj,n=2,reduction.name="obj_vf_lsi",)
obj <- RunSVD(obj,n=4,features=obj_peaks,reduction.name="obj_peaks_lsi",)
obj <- RunUMAP(obj, reduction = 'obj_peaks_lsi', dims = 1:3, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
# We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. 
# We use this graph for UMAP visualization and clustering
#obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("obj_features_pca", "obj_peaks_lsi"), 
  dims.list = list(1:20, 1:4),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", resolution =1.7, algorithm = 3, verbose = FALSE)

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

Idents(obj)<-obj$seurat_clusters
#####further annotation########
obj <- RenameIdents(
  object = obj,
  '0' = 'C2',
  '1' = 'C1',
  '2' = 'C4',
  '3' = 'C1',
  '4' = 'C1',
  '5' = 'C3',
  '6' = 'C1' )
obj@meta.data$Annotation<-factor(Idents(obj),levels=c("C1","C2","C3","C4"))

# color set 
color_for_cluster<- c("#4BA9D1",myUmapcolors[6:8])
color_for_group <- c("#476D87","#E95C59")

ORN_count<-obj@assays$SCT
barcode_label<-data.frame(barcode=colnames(obj),label=obj$Annotation)
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
barcode_label<-barcode_label[colnames(ORN_matrix),]
barcode_label<-barcode_label[order(barcode_label$label),]
# p1 cell cosine simility heatmap 
#DefaultAssay(obj)<-"integratedRNA_onecluster"
#obj<- ScaleData(obj)
#obj <- RunPCA(obj,reduction.name="pca") 
embeddings <- Embeddings(object = obj, reduction = "obj_features_pca")[,1:20]
embeddings <- embeddings[barcode_label$barcode,]
trans_dist <- 1-cosine(t(embeddings))
barcode_label_pheatmap<-data.frame(label=barcode_label$label)
rownames(barcode_label_pheatmap)<-barcode_label$barcode
col <- color_for_cluster
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
p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+guides(color='none')+
              geom_density(alpha=.25) + theme_classic()+ scale_fill_manual(values=color_for_group)+ scale_color_manual(values=color_for_group)
# t-test
p3 <- ggboxplot(data, x="type", y="var", color = "type",width=0.6,) +stat_compare_means()+guides(color = "none")+ scale_color_manual(values=color_for_group)
#raw counts heatmap 
# Heat map of expression  value 
    Idents(obj)<-obj$Annotation
    DefaultAssay(obj)<-"raw_RNA"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
    obj_data<-obj_data[rownames(barcode_label),]
    p4<-pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "#CC0000"))(100),
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

pdf("./00_Figure/Fig4/Fig4ABCD0combination_group_recluster_trans_exp.pdf",width=12,height=5)
print(add_title)
dev.off()

# Fig4H coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(obj)<-"peaks_ORN_subcluster"
  Idents(obj)<-obj$Annotation
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
  col<- color_for_cluster
  p1<-CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 200,
    annotation = TRUE,
    extend.downstream = 600,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 30,
    links=F
  )

pdf("./00_Figure/Fig4/Fig4E-combination_group_recluster_trackplot.pdf",width=10,height=10)
#p1<-p1& scale_fill_manual(values=col)
print(p1)
dev.off()


# obj_feature circle plot 
##########################################################
### chord plot
library(chorddiag)
highReceptorRNA = obj@assays$raw_RNA@counts[obj_features,]
highReceptorRNA.data = obj@assays$raw_RNA@data[obj_features,]
highReceptorRNA.expCol=highReceptorRNA[,colSums(highReceptorRNA) != 0]
highReceptorRNA.data.expCol=highReceptorRNA.data[,colSums(highReceptorRNA.data) != 0]
# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = 4, nrow = 4))
colnames(coexp.df) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df) <- rownames(highReceptorRNA.expCol)
coexp.df2 <-coexp.df
coexp.df3 <- coexp.df

# normalized expression
for (i in 1:4) {
  for (j in 1:4) {
    # print(c('i = ', i))
    # print(c('j = ', j))
    if (i==j) { coexp.df[i,j] <- 0 }
    else {
      lapply(1:length(ncol(highReceptorRNA.data.expCol)), function(x){
        print(c('i = ', i))
        print(c('j = ', j))
        iExpressedCells.mx = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1]
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        # print(jExpressedN)
        coexp.df[i,j]<<-jExpressedN
        
        iExpressedCells.mx2 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 2]
        jExpressedN2=NCOL(iExpressedCells.mx2[, iExpressedCells.mx2[j,] >= 2])
        # print(jExpressedN)
        coexp.df2[i,j]<<-jExpressedN2
        
        iExpressedCells.mx3 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1.5]
        jExpressedN3=NCOL(iExpressedCells.mx3[, iExpressedCells.mx3[j,] >= 1.5])
        # print(jExpressedN)
        coexp.df3[i,j]<<-jExpressedN3
      })
    }
  }
}
write.csv(coexp.df,"./05_ORN_cluster2/05_combination_group_recluster/coexp_1.csv")
write.csv(coexp.df2,"./05_ORN_cluster2/05_combination_group_recluster/coexp_2.csv")
write.csv(coexp.df3,"./05_ORN_cluster2/05_combination_group_recluster/coexp_1_5.csv")
print(chorddiag(as.matrix(coexp.df), groupColors = col, groupnamePadding = 10, showTicks = FALSE))
chorddiag(as.matrix(coexp.df2), groupColors = col, groupnamePadding = 20)
chorddiag(as.matrix(coexp.df3), groupColors = col, groupnamePadding = 20)


# upsetR 


# Fig2F Upset plot for 4  
# 1: Observation
ORN_count<-obj@assays$RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
library(UpSetR)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
#pdf("./00_Figure/Fig4/Fig4F-combination_group_recluster-upsetR_Observation.pdf", width=6, height=4)
#upset(data_plot, sets = c("Or63_b","LOC410623","LOC107963999","LOC100578045"), mb.ratio = c(0.5, 0.5), keep.order = TRUE,decreasing = c(FALSE,FALSE))
##upset(fromList(listInput),  mb.ratio = c(0.5, 0.5), sets = obj_features, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
##upset(fromList(listInput),  mb.ratio = c(0.5, 0.5), sets = obj_features,keep.order = TRUE, order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE))
#dev.off()
library(ComplexUpset)
list_name<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
pdf("./00_Figure/Fig4/Fig4F-combination_group_recluster-upsetR_Observation.pdf", width=6, height=4)
upset(
    data_plot, list_name,
    width_ratio=0.3, sort_sets=FALSE,
    sort_intersections=FALSE,
    intersections=list(list_name[1],list_name[1:2],list_name[1:3],list_name[1:4],
      list_name[2],list_name[2:3],list_name[2:4],list_name[3],list_name[3:4],list_name[4]
    ),
    queries=list(
        upset_query(intersect=list_name[1],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=list_name[1:2],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=list_name[1:3],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=list_name[1:4],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=list_name[2],color="#476D87",fill="#476D87"),
        upset_query(intersect=list_name[2:3],color="#476D87",fill="#476D87"),
        upset_query(intersect=list_name[2:4],color="#476D87",fill="#476D87"),
        upset_query(intersect=list_name[3],color="#E95C59",fill="#E95C59"),
        upset_query(intersect=list_name[3:4],color="#E95C59",fill="#E95C59"),
        upset_query(intersect=list_name[4],color="#E59CC4",fill="#E59CC4")
    )
)
dev.off()

saveRDS(obj,"./05_ORN_cluster2/05_combination_group_recluster/obj_recluster.rds")

obj<- readRDS("./05_ORN_cluster2/05_combination_group_recluster/obj_recluster.rds")
# for C1: promoter1 
obj_C1<- subset(obj,idents="C1")
ORN_count<-obj_C1@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
library(UpSetR)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
pdf("./00_Figure/Fig4/Fig4F-promoter1-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
#upset(data, sets=c("Or63_b","LOC410603","LOC107963999","LOC100578045"), order.by=c("degree","freq"),empty.intersections=TRUE,  mb.ratio = c(0.5, 0.5), keep.order = F)
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]#'Outside of known sets'
         ),
    queries=list(
        upset_query(intersect=obj_upset[1],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:2],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:3],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:4],color="#4BA9D1",fill="#4BA9D1")
    )
)
dev.off()

# simulation 
obj_barcode<-colnames(obj_C1)
ORN_count<-obj_C1@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
data<- as.data.frame(ORN_matrix)
data$features<- rownames(data)
data_long<-melt(data, id.vars = c("features"),
               measure.vars = c(colnames(data)[-length(colnames(data))]),
               variable.name = c('barcode'),
               value.name = 'value')
data_long<- data_long[data_long$value>0,]
listInput <- split(data_long$barcode,data_long$features)
total_cell <- length(colnames(ORN_matrix))
gene1<- names(listInput)[1]
gene2<- names(listInput)[2]
gene3<- names(listInput)[3]
gene4<- names(listInput)[4]
gene1_dropout <- 1-length(listInput[[gene1]])/length(colnames(ORN_matrix))
gene2_dropout <- 1-length(listInput[[gene2]])/length(colnames(ORN_matrix))
gene3_dropout <- 1-length(listInput[[gene3]])/length(colnames(ORN_matrix))
gene4_dropout <- 1-length(listInput[[gene4]])/length(colnames(ORN_matrix))
input <- c(
  "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout*gene4_dropout,
  "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout*gene4_dropout,
  "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout*gene4_dropout,
  "gene4" =total_cell*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene3_dropout,
  "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout*gene4_dropout, 
  "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout*gene4_dropout, 
  "gene1&gene4" =total_cell*(1-gene1_dropout)*(1-gene4_dropout)*gene2_dropout*gene3_dropout, 
  "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout*gene4_dropout, 
  "gene2&gene4" =total_cell*(1-gene2_dropout)*(1-gene4_dropout)*gene1_dropout*gene3_dropout, 
  "gene3&gene4" =total_cell*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene2_dropout, 
  "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*gene4_dropout,
  "gene1&gene2&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene4_dropout)*gene3_dropout,
  "gene1&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout, 
  "gene2&gene3&gene4"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout,
  "gene1&gene2&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout))
data <- UpSetR::fromExpression(input)
colnames(data)<- names(listInput)
obj_upset<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
pdf("./00_Figure/Fig4/Fig4F-promoter1-combination_group_recluster-upsetR_Expecption.pdf", width=8, height=4)
upset(
    data,
    c("Or63-b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[1],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:2],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:3],color="#4BA9D1",fill="#4BA9D1"),
        upset_query(intersect=obj_upset[1:4],color="#4BA9D1",fill="#4BA9D1")
    )
)
dev.off()
# four gene exp violin plot in C1 
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C1"])
data <- FetchData(obj,vars = obj_upset,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter1-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()


# for C2: promoter2 
obj_C2<- subset(obj,idents="C2")
ORN_count<-obj_C2@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./00_Figure/Fig4/Fig4F-promoter2-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[2],  color=color_for_cluster[2],fill=color_for_cluster[2]),
        upset_query(intersect=obj_upset[2:3],color=color_for_cluster[2],fill=color_for_cluster[2]),
        upset_query(intersect=obj_upset[2:4],color=color_for_cluster[2],fill=color_for_cluster[2])
    ))
dev.off()

# simulation 
obj_barcode<-colnames(obj_C2)
ORN_count<-obj_C2@assays$raw_RNA
ORN_count<-ORN_count[obj_features,]
ORN_matrix<-as.matrix(ORN_count)
data<- as.data.frame(ORN_matrix)
data$features<- rownames(data)
data_long<-melt(data, id.vars = c("features"),
               measure.vars = c(colnames(data)[-length(colnames(data))]),
               variable.name = c('barcode'),
               value.name = 'value')
data_long<- data_long[data_long$value>0,]
listInput <- split(data_long$barcode,data_long$features)
total_cell <- length(colnames(ORN_matrix))
gene1<- names(listInput)[1]
gene2<- names(listInput)[2]
gene3<- names(listInput)[3]
gene4<- names(listInput)[4]
gene1_dropout <- 1-length(listInput[[gene1]])/length(colnames(ORN_matrix))
gene2_dropout <- 1-length(listInput[[gene2]])/length(colnames(ORN_matrix))
gene3_dropout <- 1-length(listInput[[gene3]])/length(colnames(ORN_matrix))
gene4_dropout <- 1-length(listInput[[gene4]])/length(colnames(ORN_matrix))
input <- c(
  "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout*gene4_dropout,
  "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout*gene4_dropout,
  "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout*gene4_dropout,
  "gene4" =total_cell*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene3_dropout,
  "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout*gene4_dropout, 
  "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout*gene4_dropout, 
  "gene1&gene4" =total_cell*(1-gene1_dropout)*(1-gene4_dropout)*gene2_dropout*gene3_dropout, 
  "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout*gene4_dropout, 
  "gene2&gene4" =total_cell*(1-gene2_dropout)*(1-gene4_dropout)*gene1_dropout*gene3_dropout, 
  "gene3&gene4" =total_cell*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene2_dropout, 
  "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*gene4_dropout,
  "gene1&gene2&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene4_dropout)*gene3_dropout,
  "gene1&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout, 
  "gene2&gene3&gene4"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout,
  "gene1&gene2&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout))
data <- UpSetR::fromExpression(input)
colnames(data)<- names(listInput)

pdf("./00_Figure/Fig4/Fig4F-promoter2-combination_group_recluster-upsetR_Expecption.pdf", width=8, height=4)
obj_upset<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63-b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[2],  color=color_for_cluster[2],fill=color_for_cluster[2]),
        upset_query(intersect=obj_upset[2:3],color=color_for_cluster[2],fill=color_for_cluster[2]),
        upset_query(intersect=obj_upset[2:4],color=color_for_cluster[2],fill=color_for_cluster[2])
    ))
dev.off()

# four gene exp violin plot in C2
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C2"])
data <- FetchData(obj,vars = obj_upset,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter2-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()



# for C3: promoter3 
obj_C3<- subset(obj,idents="C3")
ORN_count<-obj_C3@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./00_Figure/Fig4/Fig4F-promoter3-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[3],  color=color_for_cluster[3],fill=color_for_cluster[3]),
        #upset_query(intersect=obj_upset[3:3],color=color_for_cluster[3],fill=color_for_cluster[3]),
        upset_query(intersect=obj_upset[3:4],color=color_for_cluster[3],fill=color_for_cluster[3])
    ))
dev.off()

# simulation 
obj_barcode<-colnames(obj_C3)
ORN_count<-obj_C3@assays$raw_RNA
ORN_count<-ORN_count[obj_features,]
ORN_matrix<-as.matrix(ORN_count)
data<- as.data.frame(ORN_matrix)
data$features<- rownames(data)
data_long<-melt(data, id.vars = c("features"),
               measure.vars = c(colnames(data)[-length(colnames(data))]),
               variable.name = c('barcode'),
               value.name = 'value')
data_long<- data_long[data_long$value>0,]
listInput <- split(data_long$barcode,data_long$features)
total_cell <- length(colnames(ORN_matrix))
gene1<- names(listInput)[1]
gene2<- names(listInput)[2]
gene3<- names(listInput)[3]
gene4<- names(listInput)[4]
gene1_dropout <- 1-length(listInput[[gene1]])/length(colnames(ORN_matrix))
gene2_dropout <- 1-length(listInput[[gene2]])/length(colnames(ORN_matrix))
gene3_dropout <- 1-length(listInput[[gene3]])/length(colnames(ORN_matrix))
gene4_dropout <- 1-length(listInput[[gene4]])/length(colnames(ORN_matrix))
input <- c(
  "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout*gene4_dropout,
  "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout*gene4_dropout,
  "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout*gene4_dropout,
  "gene4" =total_cell*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene3_dropout,
  "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout*gene4_dropout, 
  "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout*gene4_dropout, 
  "gene1&gene4" =total_cell*(1-gene1_dropout)*(1-gene4_dropout)*gene2_dropout*gene3_dropout, 
  "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout*gene4_dropout, 
  "gene2&gene4" =total_cell*(1-gene2_dropout)*(1-gene4_dropout)*gene1_dropout*gene3_dropout, 
  "gene3&gene4" =total_cell*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene2_dropout, 
  "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*gene4_dropout,
  "gene1&gene2&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene4_dropout)*gene3_dropout,
  "gene1&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout, 
  "gene2&gene3&gene4"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout,
  "gene1&gene2&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout))
data <- UpSetR::fromExpression(input)
colnames(data)<- names(listInput)

pdf("./00_Figure/Fig4/Fig4F-promoter3-combination_group_recluster-upsetR_Expecption.pdf", width=8, height=4)
obj_upset<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63-b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[3],  color=color_for_cluster[3],fill=color_for_cluster[3]),
        #upset_query(intersect=obj_upset[2:3],color=color_for_cluster[3],fill=color_for_cluster[3]),
        upset_query(intersect=obj_upset[3:4],color=color_for_cluster[3],fill=color_for_cluster[3])
    ))
dev.off()

# four gene exp violin plot in C3
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C3"])
data <- FetchData(obj,vars = obj_upset,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter3-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()


# for C4: promoter4
obj_C4<- subset(obj,idents="C4")
ORN_count<-obj_C4@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or63_b = names(which(ORN_matrix[4,]>0)), 
        LOC410603 = names(which(ORN_matrix[3,]>0)), 
        LOC107963999 = names(which(ORN_matrix[2,]>0)), 
        LOC100578045 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./00_Figure/Fig4/Fig4F-promoter4-combination_group_recluster-upsetR_Observation.pdf", width=8, height=4)
obj_upset<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63_b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[4],  color=color_for_cluster[4],fill=color_for_cluster[4])#,
        #upset_query(intersect=obj_upset[3:3],color=color_for_cluster[3],fill=color_for_cluster[3]),
        #upset_query(intersect=obj_upset[3:4],color=color_for_cluster[3],fill=color_for_cluster[3])
    ))
dev.off()

# simulation 
obj_barcode<-colnames(obj_C4)
ORN_count<-obj_C4@assays$raw_RNA
ORN_count<-ORN_count[obj_features,]
ORN_matrix<-as.matrix(ORN_count)
data<- as.data.frame(ORN_matrix)
data$features<- rownames(data)
data_long<-melt(data, id.vars = c("features"),
               measure.vars = c(colnames(data)[-length(colnames(data))]),
               variable.name = c('barcode'),
               value.name = 'value')
data_long<- data_long[data_long$value>0,]
listInput <- split(data_long$barcode,data_long$features)
total_cell <- length(colnames(ORN_matrix))
gene1<- names(listInput)[1]
gene2<- names(listInput)[2]
gene3<- names(listInput)[3]
gene4<- names(listInput)[4]
gene1_dropout <- 1-length(listInput[[gene1]])/length(colnames(ORN_matrix))
gene2_dropout <- 1-length(listInput[[gene2]])/length(colnames(ORN_matrix))
gene3_dropout <- 1-length(listInput[[gene3]])/length(colnames(ORN_matrix))
gene4_dropout <- 1-length(listInput[[gene4]])/length(colnames(ORN_matrix))
input <- c(
  "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout*gene4_dropout,
  "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout*gene4_dropout,
  "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout*gene4_dropout,
  "gene4" =total_cell*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene3_dropout,
  "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout*gene4_dropout, 
  "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout*gene4_dropout, 
  "gene1&gene4" =total_cell*(1-gene1_dropout)*(1-gene4_dropout)*gene2_dropout*gene3_dropout, 
  "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout*gene4_dropout, 
  "gene2&gene4" =total_cell*(1-gene2_dropout)*(1-gene4_dropout)*gene1_dropout*gene3_dropout, 
  "gene3&gene4" =total_cell*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene2_dropout, 
  "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*gene4_dropout,
  "gene1&gene2&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene4_dropout)*gene3_dropout,
  "gene1&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout, 
  "gene2&gene3&gene4"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout,
  "gene1&gene2&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout))
data <- UpSetR::fromExpression(input)
colnames(data)<- c( "LOC100578045" ,"LOC410603" ,   "Or63-b","LOC107963999")

pdf("./00_Figure/Fig4/Fig4F-promoter4-combination_group_recluster-upsetR_Expecption.pdf", width=8, height=4)
obj_upset<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
upset(
    data,
    c("Or63-b","LOC410603","LOC107963999","LOC100578045"),
    sort_sets=FALSE,
    mode = "exclusive_intersection",
    sort_intersections = FALSE,
    intersections = list(obj_upset[1:4],obj_upset[1:3],obj_upset[1:2],obj_upset[1],obj_upset[c(1,3,4)],
      obj_upset[c(1,2,4)],obj_upset[c(1,3)],obj_upset[c(1,4)],
      obj_upset[2:4],obj_upset[2:3],obj_upset[2],obj_upset[c(2,4)],obj_upset[3:4],obj_upset[3],obj_upset[4]),
    queries=list(
        upset_query(intersect=obj_upset[4],  color=color_for_cluster[4],fill=color_for_cluster[4])#,
        #upset_query(intersect=obj_upset[2:3],color=color_for_cluster[3],fill=color_for_cluster[3]),
        #upset_query(intersect=obj_upset[3:4],color=color_for_cluster[3],fill=color_for_cluster[3])
    ))
dev.off()

# four gene exp violin plot in C3
DefaultAssay(obj) <- "raw_RNA"
selected_cells <- names(Idents(obj)[Idents(obj) == "C4"])
data <- FetchData(obj,vars = obj_upset,cells = selected_cells ,slot = "counts")
long_data <- melt(data)
pdf("./00_Figure/Fig4/Fig4F-promoter4-combination_group_recluster-geneexp.pdf", width=4, height=4)
ggplot(long_data,aes(x = variable, y = value)) + geom_violin() + geom_boxplot(width=0.1,cex=1.2)+theme_classic()
dev.off()
