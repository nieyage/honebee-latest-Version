#Identify DEG by FindMarkers, Kendall tau,correlation and mutual information 
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
ORN<- readRDS("./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_order_by_tree_recall_peak.rds")
DefaultAssay(ORN)<-"raw_RNA"
dotplot_data<-read.csv("./05_ORN_cluster/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")

# 1. Find All Markers:
DefaultAssay(ORN)<-"raw_RNA"
ORN <- ScaleData(ORN,features=rownames(ORN))

markers <- FindAllMarkers(ORN, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
table(markers$cluster)
write.csv(markers,"/md01/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster/06_DEG_and_DEP/FindAllMarkers_gene.csv")

# verification
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
signif_markers <- markers[markers$p_val_adj<0.05,] 
top5 <- signif_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC);
pdf("/md01/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster/06_DEG_and_DEP/FindAllMarkers_heatmap.pdf",width=20,height=10)
DoHeatmap(object = ORN,features=top5$gene, label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = ORN,features=top5$gene, label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();
#disp.min=0,disp.max=0.2

# 2. tau cluster specific:
library(VGAM)
ORN_matrix<-as.matrix(GetAssayData(ORN));
#filter the gene only appear in a few cells 
cell_pct = function(data){
    pct<-length(data[data!=0])/length(data)
    return(pct)
}
gene_pct<-apply(ORN_matrix,1,cell_pct)
gene_pass_pct<-names(gene_pct[gene_pct>0.005])
ORN<-NormalizeData(ORN)
ORN_avg<-AverageExpression(
       ORN,
       assays = "raw_RNA",
       features = gene_pass_pct,
       return.seurat = FALSE,
       group.by = "subcluster",
       #add.ident = NULL,
       slot = "data")
ORN_avg<-ORN_avg$raw_RNA
#https://fmicompbio.github.io/swissknife/reference/specificityScore-methods.html#value-1
source("/md01/nieyg/project/honeybee/add_antenna/swissknife-master/R/tissue_specificity_score.R")
library(matrixStats)
gene_tau<-specificityScore(
  ORN_avg,
  method = c("tau", "TSI", "counts"),
  #group = ORN$subcluster,
  thresh = 0,
  expr_values = "logcounts",
  na.rm = FALSE
)
names(gene_tau)<-rownames(ORN_avg)

#plot the cluster specific pheatmap
gene_specific<-names(which(gene_tau>0.945))
gene_specific_data<-as.data.frame(ORN_avg[gene_specific,])
data<-data.frame()
for (gene in gene_specific){
    gene_cluster_avg<-gene_specific_data[gene,]
    gene_specific_cluster<-names(gene_cluster_avg[which(gene_cluster_avg==max(gene_cluster_avg))])
    data_subset<-data.frame(gene,cluster=gene_specific_cluster);
    data<-rbind(data,data_subset)
}
data$cluster<-factor(data$cluster,levels=levels(ORN))
data<-data[order(data$cluster),]
tau1_gene<-data$gene
write.csv(data,"./ORN/07_DEG_and_DEP/tau-0.945_cluster_specfic_data.csv")
#gene_specific_data<-gene_specific_data[do.call(order,gene_specific_data),]
# plot the tau cluster specfic genes heatmap 
#modify_tau_gene<-names(which(gene_tau>0.945))

ORN<-ScaleData(ORN)
pdf("./ORN/07_DEG_and_DEP/modify_tau_gene_cluster_specfic_heatmap.pdf",width=15,height=15)
DoHeatmap(object = ORN,disp.min=0,disp.max=0.05,features=tau1_gene,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = ORN,features=tau1_gene,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();

# 3. the correlation between ORx and other gene:
library(dplyr)
library(BBmisc)
library(mlr)
library(infotheo)
ORN_count<-ORN@assays$raw_RNA
ORN_count<-ORN_count[,which(colnames(ORN_count)%in%colnames(ORN))]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>10,]
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
t_ORN_matrix<-as.data.frame(t(ORN_matrix))

cal_cv=function(x){  
  y=na.omit(x)
  return(sd(y)/mean(y))
}
gene_CV<-apply(t_ORN_matrix, 2, cal_cv)

# density plot 
data<-as.data.frame(gene_CV);
range(data$gene_CV)
#[1]  0.7724479 45.0444225
library(ggplot2)
pdf("./ORN/remove_nopower/DEGandDEP/allgene_cv.pdf",width=10,height=4)
ggplot(data, aes(x=gene_CV)) + xlab("gene_CV")+
              geom_density(alpha=.25) + theme_classic() 
dev.off()
data$gene<-rownames(data)
#select gene_cv>3 as cutoff (OR gene >3)
gene_hCV<-data[data$gene_CV>3,2]
#Or_gene%in% gene_hCV

dotplot_data<-read.csv("./ORN/remove_nopower/dotplot_data_remove_nopower.csv",row.names=1)
gene2calculate<-intersect(gene_pass_pct,gene_hCV)
gene2calculate<-c(gene2calculate,Or_gene)
t_ORN_matrix<-t_ORN_matrix[,unique(gene2calculate)]

# Step: Feature Selection 
# Correlation Matrix
Or_gene<-dotplot_data$features.plot
cor_data<-cor(t_ORN_matrix)
Or_gene_cor_data<-cor_data[Or_gene,]
# Manage correalation results 
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
cluster_info$Var1<-factor(cluster_info$Var1,levels=levels(ORN))
cluster_info<-cluster_info[order(cluster_info$Var1),]
Or_cor_data<-data.frame()
for (i in 1:nrow(cluster_info)){
  if(cluster_info[i,2]==1){
    cluster<-cluster_info[i,1]
    Or_gene<-dotplot_data[dotplot_data$id==cluster,]$features.plot
    cor_value<-Or_gene_cor_data[Or_gene,]
    cor_value<-cor_value[order(cor_value,decreasing=T)];
    top10_cor<-names(cor_value[1:10]);
    Or_cor_data_subset<-data.frame(OR=Or_gene,gene=top10_cor,cluster);
    Or_cor_data<-rbind(Or_cor_data,Or_cor_data_subset);}
  if(cluster_info[i,2]>1){
    cluster<-cluster_info[i,1]
    Or_gene<-dotplot_data[dotplot_data$id==cluster,]$features.plot
    cor_value_all<-Or_gene_cor_data[Or_gene,]
    cor_value<-colSums(cor_value_all)
    cor_value<-cor_value[order(cor_value,decreasing=T)]
    # top 
    top10_cor<-names(cor_value[1:10])
    OR_char<-paste(Or_gene,collapse="_")
    Or_cor_data_subset<-data.frame(OR=OR_char,gene=top10_cor,cluster);
    Or_cor_data<-rbind(Or_cor_data,Or_cor_data_subset);
  }
}
write.csv(Or_cor_data,"./ORN/remove_nopower/DEGandDEP/correlation_top10_data.csv")

pdf("./ORN/remove_nopower/DEGandDEP/Correalation_heatmap-top10.pdf",width=15,height=10)
DoHeatmap(object = ORN,features=Or_cor_data$gene,label=TRUE,size = 3.5,group.colors =myUmapcolors) + 
scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = ORN,features=Or_cor_data$gene,label=TRUE,size = 3.5,group.colors =myUmapcolors) +
scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();

# 4. the Mutinformation between ORx and other gene:
library(infotheo)
Or_mi_data<-data.frame()
for (i in 1:nrow(cluster_info)){
  if(cluster_info[i,2]==1){
    mi_value<-data.frame() 
    cluster<-cluster_info[i,1]
    Or_gene<-dotplot_data[dotplot_data$id==cluster,]$features.plot
    Or_gene_value<-t_ORN_matrix[,Or_gene];
    for (j in 1:ncol(t_ORN_matrix)){
     mivalue<-mutinformation(Or_gene_value,t_ORN_matrix[,j])
     mi_value_subset<-data.frame(Or_gene=Or_gene,mi_value=mivalue,gene=colnames(t_ORN_matrix)[j],cluster)
     mi_value<-rbind(mi_value,mi_value_subset)
   };
   mi_value<-mi_value[order(mi_value$mi_value,decreasing=T),]
   top10_mi<-mi_value[1:10,];
   Or_mi_data<-rbind(Or_mi_data,top10_mi);}
  if(cluster_info[i,2]>1){
    mi_value<-data.frame() 
    cluster<-cluster_info[i,1]
    Or_gene<-dotplot_data[dotplot_data$id==cluster,]$features.plot
    Or_gene_value<-t_ORN_matrix[,Or_gene];
    for (j in 1:ncol(t_ORN_matrix)){
     mivalue<-mutinformation(Or_gene_value,t_ORN_matrix[,j])
     OR_char<-paste(Or_gene,collapse="_")
     mi_value_subset<-data.frame(Or_gene=OR_char,mi_value=mivalue,gene=colnames(t_ORN_matrix)[j],cluster)
     mi_value<-rbind(mi_value,mi_value_subset)
   };
   mi_value<-mi_value[order(mi_value$mi_value,decreasing=T),]
   top10_mi<-mi_value[1:10,];
   Or_mi_data<-rbind(Or_mi_data,top10_mi);
  }
}
write.csv(Or_mi_data,"./ORN/remove_nopower/DEGandDEP/MI_top10_data.csv")
pdf("./ORN/remove_nopower/DEGandDEP/MI_heatmap-top10.pdf",width=15,height=10)
DoHeatmap(object = ORN,features=Or_mi_data$gene,label=TRUE,size = 3.5,group.colors =myUmapcolors) + 
scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = ORN,features=Or_mi_data$gene,label=TRUE,size = 3.5,group.colors =myUmapcolors) +
scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
dev.off();

library(VennDiagram)
library(RColorBrewer)
markers <-read.csv("./ORN/remove_nopower/DEGandDEP/FindAllMarkers_gene.csv",row.names=1)
tau_data<-read.csv("./ORN/remove_nopower/DEGandDEP/tau-0.8_cluster_specfic_data.csv",row.names=1)
cor_data<-read.csv("./ORN/remove_nopower/DEGandDEP/correlation_top10_data.csv",row.names=1)
mi_data <-read.csv("./ORN/remove_nopower/DEGandDEP/MI_top10_data.csv",row.names=1)
FeatureSelection_data<-data.frame()
pdf("./ORN/remove_nopower/DEGandDEP/DEG_4methods_venn..pdf")
for (i in 1:nrow(cluster_info)){
    cluster<-cluster_info[i,1];
    # 1.FindAllMarkers:
    findmarkers<-markers[markers$cluster==cluster,]$gene;
    # 2.tau:
    tau<-tau_data[tau_data$cluster==cluster,]$gene;
    # 3.correlation:
    cor<-cor_data[cor_data$cluster==cluster,]$gene;
    # 4.MI:
    mi<-mi_data[mi_data$cluster==cluster,]$gene;
    Or_gene<-dotplot_data[dotplot_data$id==cluster,]$features.plot
    OR_char<-paste(Or_gene,collapse="_")
    title<-paste(cluster,OR_char,sep=": ")
    vennplot<-venn.diagram(
        x = list(findmarkers,tau,cor,mi),
        category.names = c("findmarkers","tau","top10_cor","top10_mi"),
        filename =NULL,
        fill = brewer.pal(7, "Set2")[1:4],
        alpha = 0.50,
        main=title,
        output=TRUE
      )
    print(grid.draw(vennplot))
    grid.newpage();
}
dev.off()

#ORN cluster OR heatmap
DefaultAssay(ORN)<-"SCT"
dotplot_feature<-unique(dotplot_data$features.plot)
dotplot_feature_Orco<-c(dotplot_feature,"Or2","brp")
ORN<-ScaleData(ORN,features=dotplot_feature_Orco)
#pdf("./ORN/remove_nopower/OR_gene_heatmap.pdf",width=20,height=13)
#DoHeatmap(object = ORN,features=dotplot_feature_Orco,label=TRUE,size = 3.5,group.colors =myUmapcolors) + 
#scale_fill_gradientn(colors = c( "white", "#E41A1C"))
#DoHeatmap(object = ORN,features=dotplot_feature_Orco,label=TRUE,size = 3.5,group.colors =myUmapcolors) +
#scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
library(dittoSeq)
ORN$subcluster<-factor(ORN$subcluster,levels=levels(ORN))
# heatmap add sample and subcluster 
pdf("./ORN/remove_nopower/OR_gene_heatmap.pdf",width=20,height=13)
dittoHeatmap(ORN, dotplot_feature_Orco,scaled.to.max = TRUE,annot.by = c("subcluster", "orig.ident"),cluster_row=FALSE)
dev.off();

# 4 methods upsetR
library(UpSetR)
listInput <- list(
        MI = mi_data$gene, 
        tau = tau_data$gene, 
        FindAllMarkers = markers$gene,
        correlation = cor_data$gene)
pdf("./ORN/remove_nopower/DEGandDEP/DEG_4methods_upsetR.pdf")
upset(fromList(listInput), nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()

overlap<-Reduce(intersect,listInput)
rownames(tau_data)<-tau_data$gene

tau1_gene_decription<-all_gene[which(all_gene$gene_name%in% tau1_gene),]




all_gene<-read.csv("/md01/nieyg/project/honeybee/bulk-RNAseq/blast/onetermdata_adddecription.out")
TF<-all_gene[grep("transcriptional regulator",all_gene$V10),]





