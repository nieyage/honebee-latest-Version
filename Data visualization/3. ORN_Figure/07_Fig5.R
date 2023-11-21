# Fig5A plot the tree and label 39-42:

colors_for_exp_pattern<- c("#476D87","#E95C59")
obj<- ORN
DefaultAssay(obj)<-"integratedRNA_onecluster"
embeddings <- Embeddings(object = obj, reduction = "obj_features_pca")[,1:50]
data.dims <- lapply(X = levels(x = obj), FUN = function(x) {
    cells <- WhichCells(object = obj, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = obj)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))

library(ggtree)
pdf("./00_Figure/Fig5/Fig5A-ORN-tree-cosine.pdf",width=5,height=14)
ggtree(data.tree,layout="circular") + 
geom_tiplab()+ 
geom_hilight(node=c(68,98,99),fill = "blue",alpha = 0.6)
dev.off()
# Fig5E: 
obj<- readRDS("./05_ORN_cluster2/05_combination_group_recluster/obj_recluster.rds")
Idents(obj)<-obj$Annotation
identityMapping <- c('C1' = 'C39', 'C2' = 'C40','C3' = 'C41', 'C4' = 'C42')

obj <- RenameIdents(obj, identityMapping )
obj_features<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
barcode_label<-data.frame(barcode=colnames(obj),label=Idents(obj))
DefaultAssay(obj)<-"raw_RNA"
obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))
barcode_label<- barcode_label[order(barcode_label$label),]
obj_data<- obj_data[rownames(barcode_label),]


C39_data<-obj_data[rownames(barcode_label[barcode_label$label=="C39",]),]
C40_data<-obj_data[rownames(barcode_label[barcode_label$label=="C40",]),]
C41_data<-obj_data[rownames(barcode_label[barcode_label$label=="C41",]),]
C42_data<-obj_data[rownames(barcode_label[barcode_label$label=="C42",]),]

C39_data<- C39_data[which(rowSums(C39_data)>0),]
C40_data<- C40_data[which(rowSums(C40_data)>0),]
C41_data<- C41_data[which(rowSums(C41_data)>0),]
C42_data<- C42_data[which(rowSums(C42_data)>0),]

C39_data<- C39_data[order(C39_data$LOC100578045,C39_data$LOC107963999,C39_data$LOC410603,C39_data$`Or63-b`),]
C40_data<- C40_data[order(C40_data$LOC100578045,C40_data$LOC107963999,C40_data$LOC410603,C40_data$`Or63-b`),]
C41_data<- C41_data[order(C41_data$LOC100578045,C41_data$LOC107963999,C41_data$LOC410603,C41_data$`Or63-b`),]
C42_data<- C42_data[order(C42_data$LOC100578045,C42_data$LOC107963999,C42_data$LOC410603,C42_data$`Or63-b`),]


library(pheatmap)
# 
# clusterMatrix <- function(input_matrix) {
#   # Define the clustering method and other parameters
#   clustering_method <- "complete"  # You can change this to other methods like "ward.D", "single", etc.
#   # Perform clustering
#   p <- pheatmap(
#     input_matrix,
#     clustering_method = clustering_method,
#     cluster_cols = F,
#     cluster_rows = T,
#   )
#   clustered_matrix <- input_matrix[p$tree_row$order,]
#   # Return the clustered matrix
#   return(clustered_matrix)
# }
# smooth_column <- function(col) {
#   smoothed <- numeric(length(col))
#   for (i in 2:(length(col) - 1)) {
#     smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
#   }
#   smoothed[1] <- (col[1] + col[2]) / 2
#   smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
#   return(smoothed)
# }
# 
# C39_data_clustered<- clusterMatrix(C39_data)
# C39_data_clustered_smoothed_data <- as.data.frame(lapply(C39_data_clustered, smooth_column))
# rownames(C39_data_clustered_smoothed_data)<- rownames(C39_data_clustered)
# C40_data_clustered<- clusterMatrix(C40_data)
# C40_data_clustered_smoothed_data <- as.data.frame(lapply(C40_data_clustered, smooth_column))
# rownames(C40_data_clustered_smoothed_data)<- rownames(C40_data_clustered)
# C41_data_clustered<- clusterMatrix(C41_data)
# C41_data_clustered_smoothed_data <- as.data.frame(lapply(C41_data_clustered, smooth_column))
# rownames(C41_data_clustered_smoothed_data)<- rownames(C41_data_clustered)
# C42_data_clustered<- clusterMatrix(C42_data)
# C42_data_clustered_smoothed_data <- as.data.frame(lapply(C42_data_clustered, smooth_column))
# rownames(C42_data_clustered_smoothed_data)<- rownames(C42_data_clustered)

# clustered_smoothed_data<- rbind(C39_data_clustered_smoothed_data,
# 	C40_data_clustered_smoothed_data,
# 	C41_data_clustered_smoothed_data,
# 	C42_data_clustered_smoothed_data)
# 
# 对每一列进行平滑
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[6:8])
barcode_label_pheatmap<-data.frame(label=c(rep("C1",nrow(C39_data)),rep("C2",nrow(C40_data)),rep("C3",nrow(C41_data)),rep("C4",nrow(C42_data))))
last_data_heatmap<- rbind(C39_data,C40_data,C41_data,C42_data)
rownames(barcode_label_pheatmap)<-rownames(last_data_heatmap)
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)

pdf("./00_Figure/Fig5/Fig5E_4gene_heatmap.pdf",height=3,width=8)
pheatmap(t(last_data_heatmap),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "#CC0000"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()

colors_list<- colorRampPalette(c("white", "#CC0000"))(100)

pdf("./00_Figure/Fig4/Fig4E_4gene_heatmap_pink.pdf",height=3,width=8)
pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colors_list[10:100],
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colors_list[5:100],
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colors_list[1:100],
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()




# Fig5 UCSC track:
# LG2 gene exp dotplot 
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/Unsupervised_ORN_cluster_WNN_add_subcluster.rds")
Idents(ORN)<- ORN$subcluster
DefaultAssay(ORN) <- "integratedRNA_onecluster"
object<- ORN
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
#cluster order by tree
m<-ggtree(data.tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

ORN$subcluster<-factor(ORN$subcluster,levels=cluster_order)
Idents(ORN)<-factor(ORN$subcluster,levels=cluster_order)


G2<- chemoreceptor_info_data[chemoreceptor_info_data$seqnames=="Group2",]$gene_name
receptor.dot <- DotPlot(ORN, features = G2) + #scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
dotplot_data<-receptor.dot$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>5&&dotplot_data[i,]$avg.exp.scaled >= 2.5){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
#dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% dotplot_data$features.plot])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]

dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))


DefaultAssay(ORN)<- "SCT"
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/G2_OR_UnSupervised_WNN_dotplot-signif-feature_orderbytree.pdf",width=15, height=20)
p<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p

dev.off()

data<- ORN@assays$RNA[c("Or9","Or10","Or11"),]
Or9 <- colnames(data[,which(data[1,]>0)])
Or10 <- colnames(data[,which(data[2,]>0)])
Or11 <- colnames(data[,which(data[3,]>0)])

Or9_10<- intersect(Or9,Or10)


library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN)<-"peaks"
clear_multiple_classes_order_id<- unique(dotplot_data$id)
pdf("./05_ORN_cluster2/02_second_cluster/05_add_subcluster_info/OR_without_nopower_trackplot.pdf",width=10,height=10)
for (cluster in clear_multiple_classes_order_id){
print(cluster)
obj<-subset(ORN,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
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
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
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
}
dev.off()


chrLG2:1056100-2513400


chrLG2:9976700-10160000

region1:chrLG2:1056100-1186600
library(ggplot2)
# 读取GTF文件，假设GTF文件名为"yourfile.gtf"
gtf_data <- read.csv("/Users/fraya/Documents/project/honeybee/chemoreceptor/OR_gene_df_gene_exp_pattern_Group2.csv", header = FALSE)

# 指定要绘制的区域范围，这里假设绘制chr1的1-10000范围内的基因
chr <- "Group2"
start_pos <- 1056100
end_pos <- 1186600
head(gtf_data)
## 筛选符合区域范围的基因
genes <- subset(gtf_data, V3 == chr & V4 >= start_pos & V5 <= end_pos)

ggplot(genes, aes(xmin = V4, xmax = V5, 
                  y = V3, fill = V7, label = V1, forward = 1)) +
  geom_gene_arrow() +
  facet_wrap(~ V3, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  geom_gene_label(align = "centre",min.size = 2,reflow=T)

genes


chr <- "Group2"
start_pos <- 9976700
end_pos <- 10160000
head(gtf_data)
## 筛选符合区域范围的基因
genes <- subset(gtf_data, V3 == chr & V4 >= start_pos & V5 <= end_pos)

ggplot(genes, aes(xmin = V4, xmax = V5, 
                  y = V3, fill = V7, label = V1, forward = 1)) +
  geom_gene_arrow() +
  facet_wrap(~ V3, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  geom_gene_label(align = "centre",min.size = 2,reflow=T)

genes



















