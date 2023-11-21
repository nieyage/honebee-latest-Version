
# Fig3F two promoter  (first)
obj<- subset(ORN,idents=c("p3:3"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC107965761","LOC102655285")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3G-a-coexp-heatmap-SCT.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,
     border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_features<- gene
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!="p3:3"){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  genes.use = obj_features
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3G-b-multi_OR_without_nopower_trackplot.pdf",width=12,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream =500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
set<- colors_for_exp_pattern[c(2,1)]
p1
dev.off()

# Fig3F two promoter (second)
cluster="p3:4_1" #p2:21
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot

DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-MP2-coexp-heatmap-SCT.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
# Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
##Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
## link peaks to genes
#obj <- LinkPeaks(
#  object = obj,
#  peak.assay = "peaks_ORN_subcluster",
#  expression.assay = "raw_RNA",
#  genes.use = obj_features
#)
#######Visulize track and RNA exp######
idents.plot <- Idents(obj)
## plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-MP2-multi_OR_without_nopower_trackplot2.pdf",width=12,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream =500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()
# Fig3F two promoter (third)
cluster="p2:21" #p2:21
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot

DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-MP3-coexp-heatmap-SCT.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
# Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
##Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
## link peaks to genes
#obj <- LinkPeaks(
#  object = obj,
#  peak.assay = "peaks_ORN_subcluster",
#  expression.assay = "raw_RNA",
#  genes.use = obj_features
#)
#######Visulize track and RNA exp######
idents.plot <- Idents(obj)
## plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-MP3-multi_OR_without_nopower_trackplot2.pdf",width=12,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream =2000,
  annotation = TRUE,
  extend.downstream = 2000,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()
#Fig3H single promoter (first)
cluster="p1:14"
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-a-coexp-heatmap-SCT-singlepromoter.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  genes.use = obj_features
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-b-multi_OR_without_nopower_trackplot-singlepromoter.pdf",width=10,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()


#Fig3H single promoter (second)
cluster="p4:5"
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-SP2-coexp-heatmap-SCT-singlepromoter.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  genes.use = obj_features
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-SP2-multi_OR_without_nopower_trackplot-singlepromoter.pdf",width=10,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()

# # version 2023.10.7
# Fig3H:
pdf('./00_Figure/Fig3/Fig3H-multiple-OR-FeaturePlot.pdf', width=16, height=4)
# Visualize co-expression of two features simultaneously
FeaturePlot(ORN, features = c("LOC102653695", "LOC102653615"),cols=c("lightgrey", "#E31A1C", "#4DAE49"), max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("p3_3:Or85b_LOC102655285")
dev.off()

# Fig3I:
DefaultAssay(ORN) <- "SCT"
### Plotting scatter plot: 1. single-cell level, 2. cluster level
ORN$LOC102653695_UMIs <- ORN@assays$SCT@counts['LOC102653695',]
ORN$LOC102653615_UMIs <- ORN@assays$SCT@counts['LOC102653615',]
#.....................................................................................
#  LOC107965761 vs. LOC102655285 UMI
# ....................................................................................
#  single-cell level
library(cowplot)
pmain <- ORN@meta.data %>%
  ggplot( aes(LOC102653695_UMIs, LOC102653615_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = ORN@meta.data, aes(x = LOC102653695_UMIs), fill="#B31416") +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ORN@meta.data, aes(x = LOC102653615_UMIs), fill="#4DAE49") + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.3, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.3, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/Fig3/Fig3I-3-LOC102653695vsLOC102653615_UMI.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 1000, height =1000,p3)








