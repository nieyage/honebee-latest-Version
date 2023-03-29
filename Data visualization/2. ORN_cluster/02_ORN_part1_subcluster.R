#Step1: ORN cluster by OR gene 
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

onecluster <- readRDS("./05_ORN_cluster/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")
DefaultAssay(onecluster) <- "integratedRNA_onecluster"
Idents(onecluster) <- "seurat_clusters";
all_cluster<-levels(onecluster$seurat_clusters)
# make the trans dist tree 
object <- onecluster
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

pdf("./05_ORN_cluster/02_second_cluster/second_cluster-ORN-tree-cosine.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()
pdf("./05_ORN_cluster/02_second_cluster/second_cluster-ORN-tree-cosine-nocircular.pdf",width=8,height=12)
ggtree(data.tree) + geom_tiplab()+ geom_treescale()
dev.off()
m<- ggtree(data.tree) + geom_tiplab()+ geom_treescale()

cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

 [1] "11" "24" "23" "26" "37" "1"  "4"  "38" "39" "20" "7"  "30" "12" "34" "6" 
[16] "10" "28" "29" "21" "36" "31" "25" "2"  "22" "16" "18" "0"  "3"  "8"  "42"
[31] "13" "19" "27" "41" "15" "14" "5"  "33" "32" "35" "43" "40" "9"  "17"

one_classes <- c("29 ","28","31","33","36","39","41")
multiple_stop_cluster <- as.character(c(25,34,38))
need2subcluster <- setdiff(all_cluster,c(one_classes,multiple_stop_cluster))


# 4 parts to subcluster 
part1 <- cluster_order[1:7]
part2 <- cluster_order[8:20]
part3 <- cluster_order[21:32]
part4 <- cluster_order[33:44]
need2subcluster 
part1 <- part1[which(part1%in%need2subcluster)]
part2 <- part2[which(part2%in%need2subcluster)]
part3 <- part3[which(part3%in%need2subcluster)]
part4 <- part4[which(part4%in%need2subcluster)]

part1_subcluster <- subset(onecluster,idents=part1)
part2_subcluster <- subset(onecluster,idents=part2)
part3_subcluster <- subset(onecluster,idents=part3)
part4_subcluster <- subset(onecluster,idents=part4)
part1_subcluster
part2_subcluster
part3_subcluster
part4_subcluster

# For part1_subcluster
# vst
DefaultAssay(part1_subcluster) <- "RNA"
part1_subcluster <- FindVariableFeatures(part1_subcluster, selection.method = "vst",features=100)
top500 <- head(VariableFeatures(part1_subcluster),100)
hvf.info <- HVFInfo(object = part1_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[101:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/Find_var_RNA_top100.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
#RNA analysis

DefaultAssay(part1_subcluster) <- "integratedRNA_onecluster"
part1_subcluster <- RunPCA(part1_subcluster,features= top100 )
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/ElbowPlot_top100.pdf")
ElbowPlot(part1_subcluster,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
part1_subcluster <- RunTSNE(
  object = part1_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:50
)
part1_subcluster <- FindNeighbors(object = part1_subcluster, reduction = 'pca', dims = 1:50)
part1_subcluster <- FindClusters( object = part1_subcluster, verbose = FALSE, resolution =3,algorithm = 3)
table(part1_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_ORN_cluster_WNN_part1_subcluster.pdf",width=6,height=5)
DimPlot(part1_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part1_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
DefaultAssay(part1_subcluster) <- "raw_RNA"
Idents(part1_subcluster)<-part1_subcluster$seurat_clusters
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_raw_OR_dotplot_rawRNA_part1_subcluster.pdf",width=30, height=8)
p<-DotPlot(part1_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
saveRDS(part1_subcluster,"./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_multiple_classes_part1_subcluster.rds");

# Step4: select the cluster to subcluster 
# part1_subcluster distinguish OR pipeline 
p<-DotPlot(part1_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>30&&dotplot_data[i,]$avg.exp.scaled > 2){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- as.character(cluster_info[cluster_info$Freq==1,1])

# distinguish_multi_OR
DefaultAssay(part1_subcluster)<- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/part1_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(part1_subcluster,idents=cluster);
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
  cluster    log2FC        pvalue
1       0 0.6290646  0.000000e+00
2       7 0.7937772  0.000000e+00
3       8 0.3850521 9.724969e-100
4      10 0.1069210  1.350169e-04
5      21 0.1277053  4.407509e-03


####OR log2FC
## Random select two ORs to calculate the withingroup and between group FC;
features<-unique(as.character(dotplot_data$features.plot))
# add max exp OR label
ORN_count<-part1_subcluster@assays$raw_RNA
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
barcode_label<-c()
for (i in 1:length(colnames(ORN_matrix))){
    barcode_label[[i]]=names(which(ORN_matrix[,i]== max(ORN_matrix[,i])))
    }
names(barcode_label)<-colnames(ORN_matrix)
label_all<-unlist(barcode_label)

#transcriptome distance 
embeddings <- Embeddings(object = part1_subcluster, reduction = "pca")[,1:50]
library(lsa)
trans_dist <- 1-cosine(t(embeddings))
data<-data.frame()
# select features pairs 
remaining_gene<-features;
for (gene1 in features){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    barcode<- names(label_all[which(label_all %in% c(gene1,gene2))])
    gene1_barcode<- names(label_all[which(label_all %in% c(gene1))])
    gene2_barcode<- names(label_all[which(label_all %in% c(gene2))])
    gene1_barcode<- gsub("-..","-1",gene1_barcode)
    gene2_barcode<- gsub("-..","-1",gene2_barcode)
    barcode<- gsub("-..","-1",barcode)
    #barcode_label_subset<-barcode_label[which(barcode_label$label%in% c(gene1,gene2)),]
    # extract the distance of within and between groups 
    within_group<-c(as.numeric(trans_dist[gene1_barcode,gene1_barcode]),as.numeric(trans_dist[gene2_barcode,gene2_barcode]))
    between_group<-c(as.numeric(trans_dist[gene1_barcode,gene2_barcode]),as.numeric(trans_dist[gene2_barcode,gene1_barcode]))
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
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/part1_subcluster_select_log2FC_cutoff_gene.pdf",width=10,height=5)
ggplot(data, aes(x=log2FC)) + xlab("log2FC")+
              geom_density(alpha=.25) + theme_classic() 
d <- density(data$log2FC)
d$x[which.min(abs(diff(d$y)))]
hist(data$log2FC,prob=TRUE)
lines(d, col="red", lty=2)
#v <- optimize(approxfun(d$x,d$y),interval=c(0,1))$minimum
#abline(v=v, col="blue")
#Kmeans
df<-data
km <- kmeans(df$log2FC,centers=5)
df$clust <- as.factor(km$cluster)
library(ggplot2)
ggplot(df, aes(x=log2FC)) + 
  geom_histogram(aes(fill=clust,y=..count../sum(..count..)),
                 binwidth=0.5, color="grey50")+
  stat_density(geom="line", color="red")
dev.off()

last_data<-data.frame()
for(cluster in multiple_classes){
multi_data<-dotplot_data[dotplot_data$id %in%cluster,]
cluster_features<-multi_data$features.plot
tmp_data<-data[data$gene1%in% cluster_features,]
tmp_data<-tmp_data[tmp_data$gene2%in% cluster_features,]
for(i in 1:nrow(tmp_data)){
  gene1<-tmp_data$gene1[i];
  gene2<-tmp_data$gene2[i];
  if(dotplot_data[which(dotplot_data$features.plot==gene1),4]==dotplot_data[which(dotplot_data$features.plot==gene2),4]){
    data_tmp<-tmp_data[i,];
    data_tmp$cluster<-cluster;
   last_data<-rbind(last_data,data_tmp);
  }
}
}
write.csv(last_data,"./05_ORN_cluster/02_second_cluster/01_part1_subcluster/part1_subcluster_multiOR_pair_log2FC.csv")

> last_data
             gene1          gene2     log2FC        pvalue cluster
1   LOC100577888-a LOC100577888-b 0.23625811  2.692319e-80       0
2   LOC100577888-a LOC100577888-c 0.43796226  0.000000e+00       0
3   LOC100577888-a   LOC102653865 0.40979579  0.000000e+00       0
25  LOC100577888-b LOC100577888-c 0.45614982 1.134819e-206       0
26  LOC100577888-b   LOC102653865 0.56817685  0.000000e+00       0
48  LOC100577888-c   LOC102653865 0.58538822  0.000000e+00       0
148    LOC724590-a   LOC100576282 0.54152583  0.000000e+00       7
181   LOC102653716   LOC102653810 0.43326843 5.867691e-260       8
182   LOC102653716          Or107 0.54516467 8.462671e-238       8
183   LOC102653716          Or109 0.45226605  0.000000e+00       8
196   LOC102653810          Or107 0.68238138 1.842543e-253       8
197   LOC102653810          Or109 0.51383523  0.000000e+00       8
210          Or107          Or109 0.62011439  0.000000e+00       8
298   LOC100576839 LOC100576940-a 0.06045662  4.834416e-17      21
299   LOC100576839   LOC102654183 0.02261735  1.362765e-05      21
300 LOC100576940-a   LOC102654183 0.01766221  2.037136e-01      21

## perform sub-clustering on cluster to find additional structure
Idents(part1_subcluster)<- part1_subcluster$seurat_clusters
part1_subcluster<-FindSubCluster(part1_subcluster,0,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
DefaultAssay(part1_subcluster) <- "raw_RNA"
Idents(part1_subcluster)<-part1_subcluster$sub.cluster
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/OR_dotplot_rawRNA_part1_subcluster_recluster.pdf",width=30, height=8)
p<-DotPlot(part1_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_ORN_cluster_WNN_part1_subcluster_last.pdf",width=6,height=5)
DimPlot(part1_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part1_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna")
dev.off()

saveRDS(part1_subcluster,"./05_ORN_cluster/02_second_cluster/01_part1_subcluster/second_multiple_classes_part1_subcluster.rds");

# FigS2 B 
Idents(part1_subcluster)<-part1_subcluster$sub.cluster
cluster_cellnumber<-as.data.frame(table(Idents(part1_subcluster)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(part1_subcluster))
cluster_cellnumber$color<-myUmapcolors[1:length(levels(part1_subcluster))]
cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/part1_subcluster_cellnumber.pdf",width=4,height=8)
p<-ggplot(data = cluster_cellnumber, aes_string(x = "cluster", y = "number", 
        fill = "cluster")) +  xlab(" ") + ylab("# of cells") + 
        scale_fill_manual(values = cluster_cellnumber$color) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5))+coord_flip();
p
#add gene number in plot 
p+geom_text(aes(label = number), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

DefaultAssay(part1_subcluster) <- "raw_RNA"
Idents(part1_subcluster)<-part1_subcluster$sub.cluster
p<-DotPlot(part1_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>30&&dotplot_data[i,]$avg.exp.scaled > 2){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]


DefaultAssay(part1_subcluster) <- "raw_RNA"
Idents(part1_subcluster)<-part1_subcluster$sub.cluster
pdf("./05_ORN_cluster/02_second_cluster/01_part1_subcluster/dotplot_part1_subcluster_last.pdf",width=10, height=6)
p<-DotPlot(part1_subcluster, features = unique(c(Orco,rev(as.character(dotplot_data$features.plot))))) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

