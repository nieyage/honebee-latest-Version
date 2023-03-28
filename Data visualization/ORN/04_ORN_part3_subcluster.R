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
part3 <- cluster_order[21:32]
part3 <- part3[which(part3%in%need2subcluster)]
part3_subcluster <- subset(onecluster,idents=part3)

# For part3_subcluster
# vst
DefaultAssay(part3_subcluster) <- "RNA"
part3_subcluster <- FindVariableFeatures(part3_subcluster, selection.method = "vst",features=200)
top200 <- head(VariableFeatures(part3_subcluster),200)
all_receptor_gene[which(all_receptor_gene%in% top200)]
hvf.info <- HVFInfo(object = part3_subcluster,status = TRUE)
hvf.info<-hvf.info[order(hvf.info$variance.standardized,decreasing=T),]
hvf.info$variable$vst.variable<-as.numeric(hvf.info$variable$vst.variable)
hvf.info$variable$vst.variable[201:length(hvf.info$variable$vst.variable)]=0
hvf.info$variable$vst.variable<-as.logical(hvf.info$variable$vst.variable)
var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) +  1]
hvf.info$var.status <- var.status
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/Find_var_RNA_top200.pdf",width=20,height=6)
ggplot(data = hvf.info,aes(x = mean, y =variance.standardized ) ) +
geom_point(aes(color=var.status))+xlim(0,10)
dev.off()
#RNA analysis

DefaultAssay(part3_subcluster) <- "integratedRNA_onecluster"
part3_subcluster <- RunPCA(part3_subcluster,features= top200 )
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/ElbowPlot_top200.pdf")
ElbowPlot(part3_subcluster,ndims = 50, reduction = "pca")
dev.off()
#build a tSNE visualization
part3_subcluster <- RunTSNE(
  object = part3_subcluster,
  assay = "integratedRNA_onecluster",
  min.dist = 0.001,
  verbose = TRUE,
  check_duplicates = FALSE,
  reduction.name = "tsne.rna",
  reduction.key = "rnatSNE_",
  dims = 1:40
)
part3_subcluster <- FindNeighbors(object = part3_subcluster, reduction = 'pca', dims = 1:40)
part3_subcluster <- FindClusters( object = part3_subcluster, verbose = FALSE, resolution =4,algorithm = 3)
table(part3_subcluster$seurat_clusters)

pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster.pdf",width=6,height=5)
DimPlot(part3_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part3_subcluster, label = TRUE, repel = TRUE,pt.size=0.1,reduction = "tsne.rna",group.by = "seurat_clusters")
dev.off()
DefaultAssay(part3_subcluster) <- "raw_RNA"
Idents(part3_subcluster)<-part3_subcluster$seurat_clusters
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_OR_dotplot_rawRNA.pdf",width=30, height=8)
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()
saveRDS(part3_subcluster,"./05_ORN_cluster/02_second_cluster/03_part3_subcluster/second_multiple_classes_part3_subcluster.rds");

# Step4: select the cluster to subcluster 
# part3_subcluster distinguish OR pipeline 
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>30&&dotplot_data[i,]$avg.exp.scaled > 2){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
cluster_info<-as.data.frame(table(dotplot_data$id))
all_cluster<-cluster_info$Var1;
multiple_classes <- as.character(cluster_info[cluster_info$Freq>1,1])
one_classes <- as.character(cluster_info[cluster_info$Freq==1,1])

# distinguish_multi_OR
DefaultAssay(part3_subcluster)<- "integratedRNA_onecluster"
log2FCdata<-data.frame()
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_distinguish_multi_OR.pdf",width=14,height=6)
for (cluster in multiple_classes){
print(cluster)
obj<-subset(part3_subcluster,idents=cluster);
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
1        1  0.352645335  0.000000e+00
2        2  0.562758873  0.000000e+00
3        4  0.228997815  8.548535e-93
4        5  0.281421161 9.752846e-106
5        6  0.076313988  8.247579e-13
6        8  0.325889887 4.406146e-199
7        9  0.221754196  7.919745e-97
8       10  0.180517854  7.694071e-14
9       12 -0.016581967  2.731463e-01
10      13  0.059545369  7.899671e-04
11      14  0.151805997  5.628395e-08
12      16 -0.042736358  3.295301e-02
13      19  0.275185221  1.282224e-37
14      21  0.109642908  1.289270e-06
15      23  0.007730165  1.376274e-01
16      24  0.097816608  2.019865e-02
17      25 -0.099184507  3.265791e-01


####OR log2FC
## Random select two ORs to calculate the withingroup and between group FC;
features<-unique(as.character(dotplot_data$features.plot))
# add max exp OR label
ORN_count<-part3_subcluster@assays$raw_RNA
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
embeddings <- Embeddings(object = part3_subcluster, reduction = "pca")[,1:50]
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
    gene1_barcode<- gsub("-.*","-1",gene1_barcode)
    gene2_barcode<- gsub("-.*","-1",gene2_barcode)
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
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_select_log2FC_cutoff_gene.pdf",width=10,height=5)
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
write.csv(last_data,"./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_multiOR_pair_log2FC.csv")
> last_data
             gene1        gene2        log2FC        pvalue cluster
52       LOC725861 LOC100577590  0.1781797323  0.000000e+00       1
53       LOC725861 LOC100577671  0.3409243964  0.000000e+00       1
54       LOC725861 LOC102655825  0.2673335121  0.000000e+00       1
55       LOC725861 LOC102655857  0.3354486276  0.000000e+00       1
102   LOC100577590 LOC100577671  0.1156557742 8.318566e-221       1
103   LOC100577590 LOC102655825  0.0075068282  8.195222e-03       1
104   LOC100577590 LOC102655857 -0.0007645725  9.116792e-01       1
151   LOC100577671 LOC102655825  0.3145089234  0.000000e+00       1
152   LOC100577671 LOC102655857  0.2992664200  0.000000e+00       1
199   LOC102655825 LOC102655857  0.0160257271  1.408516e-01       1
292   LOC100576783 LOC100576816  0.6713878030  0.000000e+00       2
293   LOC100576783 LOC100576881  0.6805198435  0.000000e+00       2
294   LOC100576783 LOC113218533  0.5556743059 2.794207e-286       2
295   LOC100576783         Or50  0.0654248740  1.371845e-18       2
296   LOC100576783         Or53  0.5377857653  0.000000e+00       2
337   LOC100576816 LOC100576881  0.5040731885 1.803961e-296       2
338   LOC100576816 LOC113218533  0.4034778645 1.970162e-167       2
339   LOC100576816         Or50  0.0750153328  2.219731e-25       2
340   LOC100576816         Or53  0.4740049095  0.000000e+00       2
381   LOC100576881 LOC113218533  0.4019170780 2.341466e-248       2
382   LOC100576881         Or50  0.1126317704  5.224483e-73       2
383   LOC100576881         Or53  0.5320865054  0.000000e+00       2
424   LOC113218533         Or50  0.0537310420  8.344332e-12       2
425   LOC113218533         Or53  0.5036331563  0.000000e+00       2
466           Or50         Or53  0.1387403032 1.478781e-121       2
547   LOC100577787 LOC102655367  0.3155740522 8.289481e-204       4
548   LOC100577787 LOC102655434  0.0161706855  3.870772e-01       4
549   LOC100577787 LOC102655553  0.0331300882  4.198487e-07       4
550   LOC100577787          Or4  0.2158120362 1.374052e-143       4
551   LOC100577787          Or5  0.2776444604 4.544317e-164       4
586   LOC102655367 LOC102655434  0.2124821422  5.982022e-84       4
587   LOC102655367 LOC102655553  0.2734005886 3.995798e-151       4
588   LOC102655367          Or4  0.2208188364 4.244127e-196       4
589   LOC102655367          Or5  0.2256178968 2.973644e-151       4
624   LOC102655434 LOC102655553  0.0192222855  7.687394e-03       4
625   LOC102655434          Or4  0.2283699530  6.583587e-98       4
626   LOC102655434          Or5  0.2942249504 2.250164e-117       4
661   LOC102655553          Or4  0.3373891312 7.468619e-239       4
662   LOC102655553          Or5  0.4311809854 2.143253e-260       4
697            Or4          Or5  0.1875641591 4.671380e-116       4
766   LOC102656907         Or55  0.1715201916 5.346684e-169       5
767   LOC102656907         Or56  0.1622752844 5.484244e-210       5
768   LOC102656907         Or57  0.1780940551  0.000000e+00       5
769   LOC102656907         Or58  0.1131279317  3.412628e-22       5
799           Or55         Or56  0.1030223335 1.420914e-101       5
800           Or55         Or57  0.1791275866  0.000000e+00       5
801           Or55         Or58  0.3417472176  0.000000e+00       5
831           Or56         Or57  0.1041553281 2.621208e-159       5
832           Or56         Or58  0.2801202111 6.316209e-288       5
862           Or57         Or58  0.1739824735 1.521848e-120       5
921   LOC102653637 LOC102653703  0.0397452436  1.062332e-37       6
976           Or25         Or26  0.1598360571  0.000000e+00       8
977           Or25         Or27  0.1576230344  0.000000e+00       8
1002          Or26         Or27  0.2029939013  0.000000e+00       8
1051  LOC100576984 LOC100577068  0.2195352997  0.000000e+00       9
1052  LOC100576984         Or41  0.0048998010  5.861738e-02       9
1074  LOC100577068         Or41  0.2074812501  0.000000e+00       9
1117          Or51         Or52  0.0056898629  7.431072e-03      10
1156  LOC100576914 LOC100576944  0.0088407373  1.800152e-02      12
1191  LOC100577226         Or35  0.1461908601 4.246792e-271      13
1222  LOC100577634         Or10  0.4067867954  0.000000e+00      14
1261  LOC100577715         Or14  0.0104019784  1.698742e-05      16
1306     LOC726097       Or30-b  0.1622338616 4.507210e-159      21
1317          Or11         Or19  0.3107844930 1.767445e-208      23
1324        Or30-a    LOC726535  0.3344703007  0.000000e+00      24
11561 LOC100576914 LOC100576944  0.0088407373  1.800152e-02      25

## perform sub-clustering on cluster to find additional structure
Idents(part3_subcluster)<- part3_subcluster$seurat_clusters
part3_subcluster<-FindSubCluster(part3_subcluster,1,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
part3_subcluster<-FindSubCluster(part3_subcluster,2,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
part3_subcluster<-FindSubCluster(part3_subcluster,4,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.8,algorithm = 1)
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
part3_subcluster<-FindSubCluster(part3_subcluster,5,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.5,algorithm = 1)
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
part3_subcluster<-FindSubCluster(part3_subcluster,8,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.4,algorithm = 1)
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
part3_subcluster<-FindSubCluster(part3_subcluster,9,"integratedRNA_onecluster_nn",subcluster.name = "sub.cluster",resolution = 0.4,algorithm = 1)
table(part3_subcluster$sub.cluster)
DefaultAssay(part3_subcluster) <- "raw_RNA"
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/OR_dotplot_rawRNA_part3_subcluster_recluster.pdf",width=30, height=8)
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

saveRDS(part3_subcluster,"./05_ORN_cluster/02_second_cluster/03_part3_subcluster/second_multiple_classes_part3_subcluster.rds");

# FigS2 B 
Idents(part3_subcluster)<-part3_subcluster$sub.cluster
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_last.pdf",width=6,height=5)
DimPlot(part3_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna",group.by = "orig.ident")
DimPlot(part3_subcluster,cols=myUmapcolors, label = TRUE, repel = TRUE,pt.size=0.8,reduction = "tsne.rna")
dev.off()

cluster_cellnumber<-as.data.frame(table(Idents(part3_subcluster)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(part3_subcluster))
cluster_cellnumber$color<-myUmapcolors[1:length(levels(part3_subcluster))]
cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/part3_subcluster_cellnumber.pdf",width=4,height=8)
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

DefaultAssay(part3_subcluster) <- "raw_RNA"
p<-DotPlot(part3_subcluster, features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
   #if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled > 1.5){dotplot_data[i,]$state="Yes"};
   if(dotplot_data[i,]$pct.exp>30&&dotplot_data[i,]$avg.exp.scaled > 2){dotplot_data[i,]$state="Yes"};
   #if(dotplot_data[i,]$pct.exp>20&&dotplot_data[i,]$avg.exp.scaled > 2.4){dotplot_data[i,]$state="Yes"};
 }
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]

DefaultAssay(part3_subcluster) <- "raw_RNA"
pdf("./05_ORN_cluster/02_second_cluster/03_part3_subcluster/dotplot_part3_subcluster_last.pdf",width=10, height=8)
p<-DotPlot(part3_subcluster, features = unique(c(Orco,rev(as.character(dotplot_data$features.plot))))) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

