 
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
ORN<- readRDS("./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower.rds")
DefaultAssay(ORN)<-"raw_RNA"


# Fig2A 
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
#ORN$subcluster<-factor(ORN$subcluster,levels=cluster_order)
pdf("./00_Figure/Fig2A-Supervised_ORN_cluster_WNN_remove_nopower.pdf",width=10,height=6)
DimPlot(ORN, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = F, label.size = 5, repel = TRUE)
dev.off();

# Fig2B # of cells 
Idents(ORN)<-ORN$subcluster
cluster_cellnumber<-as.data.frame(table(Idents(ORN)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(ORN))
cluster_cellnumber$color<-c(myUmapcolors,myUmapcolors)[1:length(levels(ORN))]
cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)

pdf("./00_Figure/Fig2B-remove_nopower_ORN_cluster_cellnumber.pdf",width=6,height=8)
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

# Fig2C Top4-abundant-OR-FeaturePlot
DefaultAssay(ORN)<-"raw_RNA"
pdf('./00_Figure/Fig2C-Top4-abundant-OR-FeaturePlot-1.pdf', width=16, height=5)
# Visualize co-expression of two features simultaneously
FeaturePlot(ORN, features = c("LOC102653782", "LOC102656904"),max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("p4:0")
FeaturePlot(ORN, features = c("LOC410603", "LOC107963999"),max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("25")
FeaturePlot(ORN, features = c("LOC102653637", "LOC102653703"),max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("p3:6")
dev.off()
pdf('./00_Figure/Fig2C-Top4-abundant-OR-FeaturePlot-2.pdf', width=5.5, height=5)
FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("Or12") ,order=TRUE, ncol = 1)+ggtitle("Or12")
dev.off()

# Fig2D
DefaultAssay(ORN)<-"RNA"
pdf('./00_Figure/Fig2D-Orco-FeaturePlot_1.pdf', width=9, height=8)
FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("Or2") ,order=TRUE, ncol = 1)
dev.off()

library(scCustomize)
pdf('./00_Figure/Fig2D-Orco-FeaturePlot_2.pdf', width=14, height=4)
# Visualize co-expression of two features simultaneously
p1<-FeaturePlot(ORN, features = c("LOC552552"),max.cutoff =3, order=TRUE,)
p2<-FeaturePlot(ORN, features = c("LOC726019"),max.cutoff =3, order=TRUE,)
p3<-FeaturePlot(ORN, features = c("LOC551704"),max.cutoff =3, order=TRUE,)
#Plot_Density_Joint_Only(seurat_object = ORN, features = c("Or2", "LOC552552"))
p1|p2|p2
dev.off()

receptorSCT.data = ORN@assays$RNA[dotplot_feature,]
OrcoL<- c("Or2","LOC552552","LOC726019","LOC551704")
for (Orco in OrcoL) {
  print(Orco)
  positiveL <- names( receptorSCT.data[Orco,][receptorSCT.data[Orco,] >= 1] )
  positiveL2 <- names( receptorSCT.data[Orco,][receptorSCT.data[Orco,] >= 2] )
  ORN@meta.data[paste0(Orco, '_norExp1')] <- 
    as.numeric(
      lapply(rownames(ORN@meta.data), function(x){
        # print(x)
        if (as.character(x) %in% positiveL) {
          return(1)
        } else {return(0)}
      })
    )
  
  if (Orco == 'Or2') {
    ORN@meta.data[paste0(Orco, '_norExp2')] <- 
      as.numeric(
        lapply(rownames(ORN@meta.data), function(x){
          # print(x)
          if (as.character(x) %in% positiveL2) {
            return(1)
          } else {return(0)}
        })
      )
  }
}
library(VennDiagram)
dataset1 <-  row.names(ORN@meta.data[ORN@meta.data$Or2_norExp1 == 1,])
dataset2 <-  row.names(ORN@meta.data[ORN@meta.data$LOC552552_norExp1 == 1,])
dataset3 <-  row.names(ORN@meta.data[ORN@meta.data$LOC726019_norExp1 == 1,])
dataset4 <-  row.names(ORN@meta.data[ORN@meta.data$LOC551704_norExp1 == 1,])
name1 <- 'Or2 , norm.exp > 1'
name2 <- 'LOC552552 , norm.exp > 1'
name3 <- 'LOC726019 , norm.exp > 1'
name4 <- 'LOC551704 , norm.exp > 1'
library(UpSetR)
listInput <- list(
        Or2 = dataset1, 
        LOC552552 = dataset2, 
        LOC726019 = dataset3,
        LOC551704 = dataset4)
pdf("./00_Figure/Fig2D-Orco-upsetR.pdf", width=6, height=4)
upset(fromList(listInput), nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()




# Fig2E
DefaultAssay(ORN) <- "integratedRNA_onecluster"
Idents(ORN)<-ORN$subcluster
object<-ORN;
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

#add cluster info transcript distance tree
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
pdf("./00_Figure/Fig2E-a-remove_nopower-cluster-ORN-tree-cosine.pdf",width=12,height=14)
tree <- groupOTU(data.tree, .node=multiOR_cluster)
ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale()
dev.off()
#cluster order by tree
m<-ggtree(tree) + geom_tiplab()+ geom_treescale()
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)

DefaultAssay(ORN)<-"raw_RNA"
ORN_withpower$subcluster<-factor(ORN_withpower$subcluster,levels=cluster_order)
Idents(ORN)<-factor(ORN_withpower$subcluster,levels=cluster_order)

DefaultAssay(ORN)<-"raw_RNA"
p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>80&&dotplot_data[i,]$avg.exp.scaled> 0){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp.scaled> 1){dotplot_data[i,]$state="Yes"};
  if(dotplot_data[i,]$pct.exp>25&&dotplot_data[i,]$avg.exp.scaled>2.4){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))
write.csv(dotplot_data,"./05_ORN_cluster/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
DefaultAssay(ORN)<-"raw_RNA"
pdf("./00_Figure/Fig2E-b-remove_nopower-dotplot-orderbytree.pdf",width=22, height=14)
p<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
# change color bar for p2
p2<-DotPlot(ORN,features = dotplot_feature) +  xlab('') + ylab('') +  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) +
scale_color_gradientn(colours = c('#008080', '#FF00FF'),  name = 'Average\nexpression', oob = scales::squish)
p2
dev.off()

saveRDS(ORN,"./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_orderbytree.rds")