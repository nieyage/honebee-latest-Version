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
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

DefaultAssay(ORN)<-"raw_RNA"
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C"  , "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
# Fig2A:Unsupervised clustering ORN 
onecluster <- readRDS("./05_ORN_cluster2/02_second_cluster/ORN_integrated_antenna_withOr2_second_top500.rds")
pdf("./00_Figure/Fig2/Fig2A-Unsupervised_ORN_cluster_WNN.pdf",width=6,height=6)
DimPlot(onecluster, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = T, label.size = 5, repel = TRUE) & NoLegend() 
dev.off()

# Fig2B:supervised clustering ORN 
pdf("./00_Figure/Fig2/Fig2B-Supervised_ORN_cluster_WNN_remove_nopower.pdf",width=9,height=6)
DimPlot(ORN, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = F, label.size = 5, repel = TRUE)
dev.off();

# Fig2C # of cells 
Idents(ORN)<-ORN$subcluster
cluster_cellnumber<-as.data.frame(table(Idents(ORN)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(ORN))
cluster_cellnumber$color<-c(myUmapcolors,myUmapcolors)[1:length(levels(ORN))]
cluster_cellnumber<-cluster_cellnumber[order(cluster_cellnumber$number,decreasing=F),]
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=cluster_cellnumber$cluster)

pdf("./00_Figure/Fig2/Fig2C-remove_nopower_ORN_cluster_cellnumber.pdf",width=6,height=9)
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

# Fig2D
colors_for_exp_pattern<- c("#476D87","#E95C59")
DefaultAssay(ORN) <- "integratedRNA_onecluster"
object<- ORN
DefaultAssay(object)<-"RNA"
obj<-FindVariableFeatures(object, selection.method = "vst")
top <- head(VariableFeatures(obj),500)

DefaultAssay(obj)<-"raw_RNA"
obj<-ScaleData(obj,rownames(obj))
DefaultAssay(obj)<-"integratedRNA_onecluster"
object <- RunPCA(obj,features=c(all_receptor_gene,top),reduction.name="obj_features_pca") 
embeddings <- Embeddings(object = object, reduction = "obj_features_pca")[,1:50]
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
ORN<- object
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
library(ggtree)
pdf("./00_Figure/Fig2/Fig2D-a-remove_nopower-cluster-ORN-tree-cosine.pdf",width=12,height=14)
tree <- groupOTU(data.tree, .node=multiOR_cluster)
ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
dev.off()
m<-ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)
Idents(ORN)<-factor(ORN$subcluster,levels=cluster_order)

DefaultAssay(ORN)<-"SCT"
p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp>= 1){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% dotplot_feature])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))
write.csv(dotplot_data,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
DefaultAssay(ORN)<-"SCT"

dotplot_feature<- c(Orco,dotplot_feature)
pdf("./00_Figure/Fig2/Fig2D-b-remove_nopower-dotplot-orderbytree.pdf",width=25, height=14)
p<-DotPlot(ORN,features = dotplot_feature,cols=c("lightgrey","#0000CC")) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
p
dev.off()

# > dotplot_feature[which(dotplot_feature%in% IR_gene)]
"LOC552552" "LOC726019" "LOC551704"  
ORN$subcluster<- factor(ORN$subcluster,levels=cluster_order)
saveRDS(ORN,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")

# Fig2E two single OR and two coexp OR pairs
DefaultAssay(ORN)<-"raw_RNA"
pdf('./00_Figure/Fig2/Fig2E-a-single-OR-FeaturePlot.pdf', width=9, height=4)
p1<- FeaturePlot(ORN,cols=c("lightgrey","#0000CC"), reduction = 'tsne.rna',max.cutoff = 10,features = c("Or12") ,order=TRUE, ncol = 1)+ggtitle("p1_2:Or12")
p2<- FeaturePlot(ORN,cols=c("lightgrey","#0000CC"), reduction = 'tsne.rna',max.cutoff = 10,features = c("LOC724763") ,order=TRUE, ncol = 1)+ggtitle("27:LOC724763(Or9a)")
p1|p2
dev.off()
pdf('./00_Figure/Fig2/Fig2E-b-multiple-OR-FeaturePlot.pdf', width=16, height=4)
# Visualize co-expression of two features simultaneously
FeaturePlot(ORN, features = c("LOC107965761", "LOC102655285"),cols=c("lightgrey", "#E31A1C", "#4DAE49"), max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("p3_3:Or85b_LOC102655285")
dev.off()
DefaultAssay(ORN)<-"RNA"
pdf('./00_Figure/Fig2/Fig2E-c-multiple-OR-FeaturePlot.pdf', width=5, height=8)
p1<- FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("LOC107965761") ,order=TRUE, ncol = 1)+ggtitle("Or85b")
p2<- FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("LOC102655285") ,order=TRUE, ncol = 1)+ggtitle("LOC102655285")
p1/p2
dev.off()
DefaultAssay(ORN) <- "SCT"
### Plotting scatter plot: 1. single-cell level, 2. cluster level
ORN$LOC107965761_UMIs <- ORN@assays$SCT@counts['LOC107965761',]
ORN$LOC102655285_UMIs <- ORN@assays$SCT@counts['LOC102655285',]
#.....................................................................................
#  LOC107965761 vs. LOC102655285 UMI
# ....................................................................................
#  single-cell level
library(cowplot)
pmain <- ORN@meta.data %>%
  ggplot( aes(LOC107965761_UMIs, LOC102655285_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = ORN@meta.data, aes(x = LOC107965761_UMIs), fill="#B31416") +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ORN@meta.data, aes(x = LOC102655285_UMIs), fill="#4DAE49") + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.3, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.3, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/Fig2/Fig2G-3-LOC107965761vsLOC102655285_UMI.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 1000, height =1000,p3)

library(scCustomize)
# Fig2F Orco Violin plot 
Orco<- c("Or2","LOC552552","LOC551704","LOC726019")
DefaultAssay(ORN) <- "SCT"
Idents(ORN)<-ORN$subcluster
colors_list <- c(myUmapcolors,myUmapcolors)
pdf('./00_Figure/Fig2/Fig2F-Orcocoreceptor_VlnPlot_RNA.pdf',width=25, height=6)
#VlnPlot(ORN, features = Orco, ncol = 1, pt.size = 0.1)
Stacked_VlnPlot(seurat_object = ORN, features = Orco, x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = colors_list)
dev.off()

# Fig2G Orco and Other Ir gene RMSD heatmap 

# fly Ir 


***fly and honybee coreceptor 
ls *.pdb > pdb_list.txt
for file1 in $(<pdb_list.txt)
do
 for file2 in $(<pdb_list.txt)
 do 
  python ../../RMSD.py $file1 $file2 >> RMSD_result.txt
 done
done
sed -i '/^Matrix/d' RMSD_result.txt 

RMSD <- read.table("RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("pdb_list.txt")
ID <- gsub(".pdb","",pdb$V1)
ID<-gsub(".*_","",ID)

data<- matrix(ncol=length(ID),nrow=length(ID))
rownames(data)<- ID
colnames(data)<- ID
i=1
for (row in 1:length(ID)){
    for(col in 1:length(ID)){
        data[row,col]=RMSD[i];
        i=i+1
    }
}

# plot heatmap 
library(pheatmap)
label_pheatmap<- data.frame(species=c(rep("D.melanogaster",5),rep("Apis mellifera",4)))
rownames(label_pheatmap) <- colnames(data)
data[data==0]<- 0.5
data<- log2(data)
ann_colors<-list(
species = c("D.melanogaster"="#7985B6", "Apis mellifera"="#C7B6E1"))

pdf("./honeybee_OR_feature_RMSD.pdf",width=12,height=12)
pheatmap(data,
      cluster_cols = T,
      cluster_rows = T,
      color = colorRampPalette(c("#1C214F","#3082BD","#F1ECEC"))(100),
      #annotation_col = label_pheatmap,
      #annotation_colors = ann_colors,
      #annotation_row = label_pheatmap,
      #annotation_legend = TRUE,
      show_rownames=T,
      show_colnames=T
 )
dev.off()

# version 2023.10.7
# Fig2C: cluster trans dist tree + chemoreceptor gene heatmap + cluster number 

# make number 1:60 to replace subcluster id 

ORN$cell_group <- as.numeric(factor(ORN$subcluster))
Idents(ORN)<- ORN$cell_group

cell_group<- levels(Idents(ORN))
subcluster<- levels(ORN$subcluster)
names(cell_group)<- subcluster

# Fig2C-1: cluster trans dist tree

colors_for_exp_pattern<- c("#476D87","#E95C59")
DefaultAssay(ORN) <- "integratedRNA_onecluster"
object<- ORN
DefaultAssay(object)<-"RNA"
obj<-FindVariableFeatures(object, selection.method = "vst")
top <- head(VariableFeatures(obj),500)
DefaultAssay(obj)<-"raw_RNA"
obj<-ScaleData(obj,rownames(obj))
DefaultAssay(obj)<-"integratedRNA_onecluster"
object <- RunPCA(obj,features=c(all_receptor_gene,top),reduction.name="obj_features_pca") 
embeddings <- Embeddings(object = object, reduction = "obj_features_pca")[,1:50]
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
ORN<- object
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])

multiOR_cluster<- cell_group[multiOR_cluster]

library(ggtree)
pdf("./00_Figure/Fig2/Fig2C-1-remove_nopower-cluster-ORN-tree-cosine.pdf",width=5,height=14)
tree <- groupOTU(data.tree, .node=multiOR_cluster)
ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
dev.off()

# Fig2C-2: chemoreceptor gene heatmap
m<-ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)
Idents(ORN)<-factor(ORN$cell_group,levels=cluster_order)

# change the OR gene name 
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")

DefaultAssay(ORN)<- "SCT"
all_gene<- rownames(ORN)

RenameGenesSeurat_SCT <- function(obj = ls.Seurat[[i]], newnames = tmp) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts and @data ")
  SCT <- obj@assays$SCT
  if (nrow(SCT) == length(newnames)) {
    if (length(SCT@counts)) SCT@counts@Dimnames[[1]]            <- newnames
    if (length(SCT@data)) SCT@data@Dimnames[[1]]                <- newnames
    #if (length(SCT@scale.data)) SCT@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(SCT) != nrow(newnames)"}
  obj@assays$SCT <- SCT
  return(obj)
}

for (i in 1:length(all_gene)){
    if(all_gene[i] %in% ORgene_name_trans$OR_gene){
        tmp=ORgene_name_trans[which(ORgene_name_trans$OR_gene==all_gene[i]),]$last_name;
        all_gene[i]=tmp
    }
}

ORN_trans <- RenameGenesSeurat_SCT(ORN,newnames = all_gene)

DefaultAssay(ORN)<-"raw_RNA"
p<-DotPlot(ORN,features = all_receptor_gene) +  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) 
dotplot_data<-p$data;
#filter dotplot ;
dotplot_data$state<-"No";
for(i in 1:nrow(dotplot_data)){
  if(dotplot_data[i,]$pct.exp>35&&dotplot_data[i,]$avg.exp>= 1){dotplot_data[i,]$state="Yes"};
}
dotplot_data<-dotplot_data[which(dotplot_data$state=="Yes"),]
dotplot_data<-dotplot_data[-which(dotplot_data$features.plot%in% Orco),];
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf<- gtf[gtf$type=="transcript",]
gtf_data<- as.data.frame(gtf[gtf$gene_name%in% all_receptor_gene])
gtf_data<-gtf_data[order(gtf_data$seqnames,gtf_data$start),]
gene_order<- unique(gtf_data$gene_name)

dotplot_data$features.plot<- factor(dotplot_data$features.plot,levels=gene_order)
dotplot_data<- dotplot_data[order(dotplot_data$id),]
dotplot_feature<-unique(c(Orco,rev(as.character(dotplot_data$features.plot))))

for (i in 1:length(dotplot_feature)){
    if(dotplot_feature[i] %in% ORgene_name_trans$OR_gene){
        tmp=ORgene_name_trans[which(ORgene_name_trans$OR_gene==dotplot_feature[i]),]$last_name;
        dotplot_feature[i]=tmp
    }
}

DefaultAssay(ORN_trans)<-"SCT"
ORN_trans<- ScaleData(ORN_trans,features=dotplot_feature)
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

pdf("./00_Figure/Fig2/Fig2C-2-remove_nopower-chemoreceptor_gene_heatmap-orderbytree.pdf",width=15, height=14)
DoHeatmap(ORN_trans,size = 4,angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = blueYellow)
DoHeatmap(ORN_trans,disp.max = 1,slot = "data",size = 4,angle = 315, features = dotplot_feature,group.colors =c(myUmapcolors,myUmapcolors))+ scale_fill_gradientn(colors = blueYellow)
dev.off()

# Fig2B:supervised clustering ORN 
pdf("./00_Figure/Fig2/Fig2B-changeID-Supervised_ORN_cluster_WNN_remove_nopower.pdf",width=8,height=6)
DimPlot(ORN, cols=c(myUmapcolors,myUmapcolors), reduction = "tsne.rna",  label = F, label.size = 5, repel = TRUE)
dev.off();

# Fig2C-3 # of cells in cluster
Idents(ORN)<-ORN$cell_group
cluster_cellnumber<-as.data.frame(table(Idents(ORN)))
colnames(cluster_cellnumber)<-c("cluster","number")
cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=levels(ORN))
cluster_cellnumber$color<-c(myUmapcolors,myUmapcolors)[1:length(levels(ORN))]

cluster_cellnumber$cluster<- factor(cluster_cellnumber$cluster,levels=as.character(1:60))

pdf("./00_Figure/Fig2/Fig2C-3-remove_nopower_ORN_cluster_cellnumber.pdf",width=4,height=9)
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

# Fig2F: cross species expression pattern statics 
# single OR  vs multiple OR cluster porportion 
single_OR_cluster<- length(levels(ORN))-length(multiOR_cluster)
multiple_OR_cluster<- length(multiOR_cluster)
Apis_mellifera<-c(single_OR_cluster,multiple_OR_cluster);

#fly need to use the public datasets
Drosophila_melanogaster<-c(40,5)

# total:42
Aedes_aegypti<-c(13,26);

cross_species_cluster_number<-data.frame(species=c(rep("Apis mellifera",2),rep("D.melanogaster",2),rep("Ae.aegypti",2)),
  exp_pattern=rep(c("single_OR","multiple_OR"),3),
  cluster_number=c(Apis_mellifera,Drosophila_melanogaster,Aedes_aegypti))
cross_species_cluster_number$exp_pattern<-factor(cross_species_cluster_number$exp_pattern,levels=c("single_OR","multiple_OR"));
cross_species_cluster_number$species<-factor(cross_species_cluster_number$species,levels=c("Apis mellifera","D.melanogaster","Ae.aegypti"))

pdf("./00_Figure/Fig2/Fig2F-cross_species_expression_pattern_statics.pdf",width=4,height=4)
ggplot(data = cross_species_cluster_number, aes_string(x = "species", y = "cluster_number", 
        fill = "exp_pattern")) +  xlab(" ") + ylab("% Percent of cells") + 
        scale_fill_manual(values = colors_for_exp_pattern) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));

p<-ggplot(data = cross_species_cluster_number, aes_string(x = "species", y = "cluster_number", 
        fill = "exp_pattern")) +  xlab(" ") + ylab("# of cluster") + 
        scale_fill_manual(values = colors_for_exp_pattern) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = cluster_number), size = 3, hjust = 0.5, vjust = 3, position = "stack") 
dev.off();




## Fig2F Upset plot for 4 coreceptor barcode 
## 1: Observation
#OrcoL<- c("Or2","LOC552552","LOC551704","LOC726019")
#ORN_count<-ORN@assays$RNA
#ORN_count<-ORN_count[which(rownames(ORN_count)%in%OrcoL),]
#ORN_matrix<-as.matrix(ORN_count)
#ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
#ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
#
#library(UpSetR)
#listInput <- list(
#        Or2 = names(which(ORN_matrix[4,]>0)), 
#        LOC552552 = names(which(ORN_matrix[2,]>0)), 
#        LOC551704 = names(which(ORN_matrix[1,]>0)), 
#        LOC726019 = names(which(ORN_matrix[3,]>0)))
#pdf("./00_Figure/Fig2F-Orco-upsetR_Observation.pdf", width=6, height=4)
#upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
#upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), sets = OrcoL,keep.order = TRUE, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
#upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), sets = OrcoL,keep.order = TRUE, order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE))
#dev.off()

## 2: Expectation by random dropout
## the simulation data 
#total_cell <- length(colnames(ORN_matrix))
#Or2_dropout <- 1-length(listInput$Or2)/length(colnames(ORN_matrix))
#LOC552552_dropout <- 1-length(listInput$LOC552552)/length(colnames(ORN_matrix))
#LOC551704_dropout <- 1-length(listInput$LOC551704)/length(colnames(ORN_matrix))
#LOC726019_dropout <- 1-length(listInput$LOC726019)/length(colnames(ORN_matrix))
#input <- c(
#  "Or2"=total_cell*(1-Or2_dropout)*LOC552552_dropout*LOC551704_dropout*LOC726019_dropout,
#  "LOC552552"=total_cell*(1-LOC552552_dropout)*Or2_dropout*LOC551704_dropout*LOC726019_dropout,
#  "LOC551704"=total_cell*(1-LOC551704_dropout)*Or2_dropout*LOC552552_dropout*LOC726019_dropout,
#  "LOC726019"=total_cell*(1-LOC726019_dropout)*Or2_dropout*LOC552552_dropout*LOC551704_dropout,
#  "Or2&LOC552552"=total_cell*(1-Or2_dropout)*(1-LOC552552_dropout)*LOC551704_dropout*LOC726019_dropout, 
#  "Or2&LOC551704"=total_cell*(1-Or2_dropout)*(1-LOC551704_dropout)*LOC552552_dropout*LOC726019_dropout, 
#  "Or2&LOC726019"=total_cell*(1-Or2_dropout)*(1-LOC726019_dropout)*LOC551704_dropout*LOC552552_dropout, 
#  "LOC552552&LOC551704"=total_cell*(1-LOC552552_dropout)*(1-LOC551704_dropout)*Or2_dropout*LOC726019_dropout, 
#  "LOC552552&LOC726019"=total_cell*(1-LOC552552_dropout)*(1-LOC726019_dropout)*Or2_dropout*LOC551704_dropout, 
#  "LOC551704&LOC726019"=total_cell*(1-LOC551704_dropout)*(1-LOC726019_dropout)*LOC552552_dropout*Or2_dropout, 
#  "Or2&LOC552552&LOC551704"=total_cell*(1-Or2_dropout)*(1-LOC552552_dropout)*(1-LOC551704_dropout)*LOC726019_dropout,
#  "Or2&LOC552552&LOC726019"=total_cell*(1-Or2_dropout)*(1-LOC552552_dropout)*(1-LOC726019_dropout)*LOC551704_dropout,
#  "Or2&LOC726019&LOC551704"=total_cell*(1-Or2_dropout)*(1-LOC726019_dropout)*(1-LOC551704_dropout)*LOC552552_dropout, 
#  "LOC552552&LOC551704&LOC726019"=total_cell*(1-LOC552552_dropout)*(1-LOC551704_dropout)*(1-LOC726019_dropout)*Or2_dropout,
#  "Or2&LOC552552&LOC551704&LOC726019"=total_cell*(1-Or2_dropout)*(1-LOC552552_dropout)*(1-LOC551704_dropout)*(1-LOC726019_dropout)
#)
#data <- UpSetR::fromExpression(input)
#pdf("./00_Figure/Fig2F-Orco-upsetR_Expectation.pdf", width=6, height=4)
#upset(data, nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
#upset(data, nintersects = 30, mb.ratio = c(0.5, 0.5), sets = OrcoL,keep.order = TRUE, order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
#upset(data, nintersects = 30, mb.ratio = c(0.5, 0.5), sets = OrcoL,keep.order = TRUE, order.by = c("freq", "degree"), decreasing = c(TRUE,TRUE))
#dev.off()
