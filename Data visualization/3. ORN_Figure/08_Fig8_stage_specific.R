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

## stage specific cluster
ORN$orig.ident<-factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
data<-as.data.frame(table(ORN$orig.ident,ORN$subcluster))
ORN_sample_number<-as.data.frame(table(ORN$orig.ident))
NE_all_number<-ORN_sample_number[which(ORN_sample_number$Var1=="NE"),2]
Nurse_all_number<-ORN_sample_number[which(ORN_sample_number$Var1=="Nurse"),2]
Forager_all_number<-ORN_sample_number[which(ORN_sample_number$Var1=="Forager"),2]
total<-NE_all_number+Nurse_all_number+Forager_all_number
#  chisq.test;
#  fisher.test more sensitive 
fisher<-data.frame()
for(i in unique(data$Var2)){
    Clu<-data[data$Var2==i,];
    Clu<-Clu$Freq
    all<-c(NE_all_number,Nurse_all_number,Forager_all_number)
    df <- as.data.frame(rbind(Clu,all))
    colnames(df) <- c("NE","Nurse","Forager")
    print(i)
    #fisher.test(df)
    #chisq.test(df)
    p_value<-fisher.test(df)$p.value;
    #p_value<-chisq.test(df)$p.value;
    #X_squared<-fisher.test(df)$statistic;
    fisher[1,i]<-p_value;
    #chisq[2,i]<-X_squared;
}
rownames(fisher)<-c("p_value")
Idents(ORN)<-ORN$subcluster
#  sample propotion in each cluster;
library(dittoSeq)
p<-dittoBarPlot(ORN, "orig.ident", group.by = "subcluster",data.out = TRUE)
# change levels
p$data$grouping<-factor(p$data$grouping,levels=levels(factor(ORN$subcluster)))
p$data$label<-factor(p$data$label,levels=c("NE","Nurse","Forager"))
por_data<-p$data;
por_data$grouping<-as.character(por_data$grouping)
por_data<-por_data[por_data$grouping%in%colnames(fisher),]
por_data$p_value<-100
for(i in 1:length(por_data$grouping)){
  df<-fisher[,por_data$grouping[i]];
  por_data$p_value[i]=df[1]
  }
por_data$signif<-" "
for(i in 1:length(por_data$grouping)){
  if(por_data$p_value[i]<0.05){
    por_data$signif[i]="*"
  }
  };
por_data$grouping<-factor(por_data$grouping,levels=levels(Idents(ORN)))
one_class_por_data<-por_data

pdf("./05_ORN_cluster2/06_stage_specific/stage_specific_cluster_proportion.pdf",height=2,width=12)
ggplot(data = por_data, aes_string(x = "grouping", y = "percent", 
        fill = "label")) +  xlab("orig.ident") + ylab("Percent of cells") + 
        scale_fill_manual(values = c("#5CC8F2","#009E73","#E69F00")) + 
        geom_col()+ scale_y_continuous(breaks = c(0, Forager_all_number/total,1-NE_all_number/total, 1))+ coord_cartesian(ylim = c(0,1))+
geom_hline(aes(yintercept=Forager_all_number/total),linetype="dashed",col="black")+
geom_hline(aes(yintercept=1-NE_all_number/total),linetype="dashed",col="black")+
geom_text(aes(label=signif,y =0))+theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0.5))
dev.off()

## Stage specific 
signif<-por_data[por_data$p_value<0.05,]
global_ratio<-c(NE_all_number/total,Nurse_all_number/total,Forager_all_number/total)
signif$global_ratio<-c(rep(global_ratio,nrow(signif)/3))
signif$FC<-signif$percent/signif$global_ratio;

#all_class_Data<-rbind(one_class_data,multiple_classes_data)
cluster<-unique(signif$grouping)
#cluster<-cluster[-12]
all_class_Data<-dotplot_data
cluster_specific<-data.frame()
for (i in cluster){
  signif_subset<-signif[which(signif$grouping==i),];
  sample_specific<-signif_subset[which(signif_subset$FC==max(signif_subset$FC)),]$label
  OR_specfic<-all_class_Data[which(all_class_Data$id==i),]$features.plot;
  if(length(OR_specfic)!=1){
    OR_specfic<-paste(OR_specfic,collapse=",")
  }
  cluster_specific_subset<-data.frame(cluster=i,OR_specfic,sample_specific)
  cluster_specific<-rbind(cluster_specific,cluster_specific_subset)
} 
write.csv(cluster_specific,"./05_ORN_cluster2/06_stage_specific/cluster_specific.csv")


library(scCustomize)
# this is batch effect in NE,Nurse and Forager
DefaultAssay(ORN)<-"RNA";
#ORN<-NormalizeData(ORN)
#pdf("./05_ORN_cluster2/06_stage_specific/global_batch_profile_among_stage.pdf",width=10,height=20)
#p1 <- Stacked_VlnPlot(seurat_object = ORN, features = c("nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","nCount_raw_RNA","nFeature_raw_RNA","nCount_SCT","nFeature_SCT"),group.by ="orig.ident", 
#   x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
#p2 <- Stacked_VlnPlot(seurat_object = ORN, features = c("percent.mt","nucleosome_signal","TSS.enrichment","nCount_peaks","nFeature_peaks","nCount_peaks_ORN_subcluster","nFeature_peaks_ORN_subcluster"),group.by ="orig.ident", 
#  x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
#p1|p2
#dev.off()


specfic <- cluster_specific$cluster
# SCT 
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_OR_exp_SCT_data.pdf",width=5,height=10)
for(cluster in specfic){
	obj<-subset(ORN,idents=cluster);
	obj_feature<- strsplit(cluster_specific[cluster_specific$cluster==cluster,2],",")[[1]]
  #p<-VlnPlot(obj,cols =c("#5CC8F2","#009E73","#E69F00"),features=obj_feature,group.by ="orig.ident",stack = F,pt.size = 1)
	p1<-Stacked_VlnPlot(seurat_object = obj,assay = "SCT",slot="data",features = obj_feature,group.by ="orig.ident", x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
	print(p1)
}
dev.off()
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_OR_exp_SCT_counts.pdf",width=5,height=15)
for(cluster in specfic){
  obj<-subset(ORN,idents=cluster);
  obj_feature<- strsplit(cluster_specific[cluster_specific$cluster==cluster,2],",")[[1]]
  #p<-VlnPlot(obj,cols =c("#5CC8F2","#009E73","#E69F00"),features=obj_feature,group.by ="orig.ident",stack = F,pt.size = 1)
  p1<-Stacked_VlnPlot(seurat_object = obj,assay = "SCT",slot="counts",features = obj_feature,group.by ="orig.ident", x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
  print(p1)
}
dev.off()
pdf("./05_ORN_cluster2/06_stage_specific/OR2_exp_SCT_counts.pdf",width=20,height=15)
p1 <- Stacked_VlnPlot(ORN,colors_use =c("#5CC8F2","#009E73","#E69F00"),assay = "SCT",slot="counts",features=Orco,group.by ="orig.ident",stack = F,pt.size = 1)
p2 <- Stacked_VlnPlot(ORN,colors_use =c("#5CC8F2","#009E73","#E69F00"),assay = "SCT",slot="data",features=Orco,group.by ="orig.ident",stack = F,pt.size = 1)
p3 <- Stacked_VlnPlot(ORN,colors_use =c("#5CC8F2","#009E73","#E69F00"),assay = "RNA",slot="counts",features=Orco,group.by ="orig.ident",stack = F,pt.size = 1)
p4 <- Stacked_VlnPlot(ORN,colors_use =c("#5CC8F2","#009E73","#E69F00"),assay = "RNA",slot="data",features=Orco,group.by ="orig.ident",stack = F,pt.size = 1)
p5 <- Stacked_VlnPlot(ORN,colors_use =c("#5CC8F2","#009E73","#E69F00"),assay = "raw_RNA",slot="counts",features=Orco,group.by ="orig.ident",stack = F,pt.size = 1)
p6 <- Stacked_VlnPlot(ORN,colors_use =c("#5CC8F2","#009E73","#E69F00"),assay = "raw_RNA",slot="data",features=Orco,group.by ="orig.ident",stack = F,pt.size = 1)
p1|p2|p3|p4|p5|p6
dev.off()

# OR exp in their cluster among stages
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_OR_exp_raw_RNA_data_modify.pdf",width=10,height=10)
for(cluster in specfic){
  obj<-subset(ORN,idents=cluster);
  obj_feature<- strsplit(cluster_specific[cluster_specific$cluster==cluster,2],",")[[1]]
  #for(i in 1:length(obj_feature))
  #p<-VlnPlot(obj,cols =c("#5CC8F2","#009E73","#E69F00"),features=obj_feature,group.by ="orig.ident",stack = F,pt.size = 1)
  p1<-VlnPlot_scCustom(num_columns = 2,seurat_object = obj,assay = "SCT",slot="counts",features = obj_feature,group.by ="orig.ident", colors_use = c("#5CC8F2","#009E73","#E69F00"))

  print(p1)
}
dev.off()
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_OR_exp_raw_RNA_data_modify.pdf",width=5,height=3)
obj<-subset(ORN,idents="p4:0_0");
VlnPlot_scCustom(num_columns = 1,seurat_object = obj,assay = "SCT",slot="counts",features = "Or56",group.by ="orig.ident", colors_use = c("#5CC8F2","#009E73","#E69F00"))
obj<-subset(ORN,idents="p4:0_2");
VlnPlot_scCustom(num_columns = 1,seurat_object = obj,assay = "SCT",slot="counts",features = "Or55",group.by ="orig.ident", colors_use = c("#5CC8F2","#009E73","#E69F00"))
obj<-subset(ORN,idents="p2:17");
VlnPlot_scCustom(num_columns = 1,seurat_object = obj,assay = "SCT",slot="counts",features = "Or109",group.by ="orig.ident", colors_use = c("#5CC8F2","#009E73","#E69F00"))
dev.off()

#sample cell number in obj 
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_cellnumber.pdf",width=4,height=4)
for(cluster in specfic){
	obj<-subset(ORN,idents=cluster);
	data<-as.data.frame(table(obj$orig.ident))
    colnames(data)<-c("Sample","cellnumber");
    p<- ggplot(data, aes(x=Sample,y=cellnumber,fill=Sample )) +  
        geom_bar(stat = "identity" ) +
        scale_fill_manual(values = c("#5CC8F2","#009E73","#E69F00") ) +
        theme_light();
    print(p)
}

dev.off()



# stage DEG in cluster 

library(AnnotationHub)
library(biomaRt)
library(dplyr)
library(goseq)
library(DOSE)
library(stringr)
library(clusterProfiler)
library(GenomicRanges)
library(AnnotationDbi)
Apis_mellifera.OrgDb <-loadDb("/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/Apis_mellifera_AH102515.OrgDb")
columns(Apis_mellifera.OrgDb)

ORN<-ScaleData(ORN,features=rownames(ORN))
for(cluster in specfic){
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "raw_RNA";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
  markers <- markers[which(markers$p_val<0.05),]
  write.csv(markers,paste0("./05_ORN_cluster2/06_stage_specific/Cluster",cluster,"_specific_markers.csv",sep=""))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
  obj<-ScaleData(obj,features=rownames(obj))
  p<- DoHeatmap(object = obj,features=top10$gene,label=T, group.colors =c("#E69F00","#55B4E9","#009E73"),disp.min = -0.5,disp.max = 0.5,size = 2,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/Cluster",cluster,"_specific_heatmap.pdf",sep=""))
  print(p)
  dev.off()
  # Go and KEGG 
for (i in unique(markers$cluster)){
  print(i);
  gene_in_cluster<-markers[markers$cluster==i,7]
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  #GO
  pdf(paste0("Cluster",i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.1,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("Cluster",i,"_GO_",type,".csv"))
      }
  }
  dev.off()
}
}



all_ego<-rbind(NE_ego,Nurse_ego,Forager_ego)
write.csv(all_ego,"./ORN/06_stage_specific/Global_Stage_specific_markers_GO.csv")

top5_term <- all_ego %>% group_by(celltype) %>% top_n(n = 5, wt = pvalue);
top5_term$celltype<-factor(top5_term$celltype,levels=c("NE","Nurse","Forager"))
library(ggplot2)
pdf("./ORN/06_stage_specific/Global_Stage_specific_markers_GO.pdf",width=10,height=10)
p <- ggplot(top5_term,aes(y=Count,x=Description,fill=pvalue)) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(celltype~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 14),
            legend.position="right",
            legend.title = element_text(size=18),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
p
dev.off()



### true  or false DEG 
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# random select 3 cluster as control 
pdf("./05_ORN_cluster2/06_stage_specific/cluster_DEG_pval_MAplot-raw_RNA.pdf",width=24,height=6)
for (cluster in levels(ORN)){
  print(cluster)
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "raw_RNA";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  # p distribution 
  p_val_data<- data.frame()
  for(nrow in 1:nrow(markers)){
    p_val<- markers[nrow,]$p_val;
    filter_data <- markers[markers$p_val<p_val,]
    number<- as.data.frame(table(filter_data$cluster))
    p_val_data_subset<- data.frame(p_val,NE=number[1,2],Nurse=number[2,2],Forager=number[3,2])
    p_val_data<- rbind(p_val_data,p_val_data_subset)
  }
  p_val_data_long<- melt(p_val_data,"p_val")
  colnames(p_val_data_long)<- c("p_val","Stage","Freq")
  p1<-ggplot(p_val_data_long, aes(x=p_val, y=Freq, color=Stage,shape=Stage,labels=p_val)) + 
  geom_line()+
  geom_point(size=4)+
   theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 20),
        legend.text=element_text(size = 20),
        legend.title=element_text(size = 20),
        legend.position="top")+theme_bw()+
   ggtitle(cluster)+
  scale_color_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]))
   # MA plot
   MA_data<- markers[,c(1,2,6,7)]
   Avg <- AverageExpression(obj,features=MA_data$gene,assays = "raw_RNA")
   MA_data<- MA_data[names(rowMeans(Avg$raw_RNA)),]
   MA_data$baseMean<- rowMeans(Avg$raw_RNA)
   MA_data$mlgPval<- -log10(MA_data$p_val)
   xlim_max<- summary(MA_data$baseMean)[5]+100;
   p2<- ggplot(MA_data, aes(x = baseMean, y = avg_log2FC)) +
        geom_point(aes(color = cluster, size = mlgPval), show.legend = TRUE)+
        scale_radius(range = c(.1, 2)) +
        xlim(c(0,xlim_max))+theme_bw()+
        scale_color_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]))
   data<-as.data.frame(table(obj$orig.ident))
   colnames(data)<-c("Sample","cellnumber");
    p4<- ggplot(data, aes(x=Sample,y=cellnumber,fill=Sample )) +  
        geom_bar(stat = "identity" ) +
        scale_fill_manual(values =  c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]) ) +
        theme_light();
   p3<- p1+p2+p4
   print(p3)
}
dev.off()

# SCT matrix (counts)
# random select 3 cluster as control 
pdf("./05_ORN_cluster2/06_stage_specific/cluster_DEG_pval_MAplot-SCT.pdf",width=18,height=6)
for (cluster in levels(ORN)){
  print(cluster)
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "SCT";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  # p distribution 
  p_val_data<- data.frame()
  for(nrow in 1:nrow(markers)){
    p_val<- markers[nrow,]$p_val;
    filter_data <- markers[markers$p_val<p_val,]
    number<- as.data.frame(table(filter_data$cluster))
    p_val_data_subset<- data.frame(p_val,NE=number[1,2],Nurse=number[2,2],Forager=number[3,2])
    p_val_data<- rbind(p_val_data,p_val_data_subset)
  }
  p_val_data_long<- melt(p_val_data,"p_val")
  colnames(p_val_data_long)<- c("p_val","Stage","Freq")
  p1<-ggplot(p_val_data_long, aes(x=p_val, y=Freq, color=Stage,shape=Stage,labels=p_val)) + 
  geom_line()+
  geom_point(size=4)+
   theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 20),
        legend.text=element_text(size = 20),
        legend.title=element_text(size = 20),
        legend.position="top")+theme_bw()+ggtitle(cluster)+
  scale_color_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]))
   # MA plot
   MA_data<- markers[,c(1,2,6,7)]
   Avg <- AverageExpression(obj,features=MA_data$gene,assays = "SCT")
   MA_data<- MA_data[names(rowMeans(Avg$SCT)),]
   MA_data$baseMean<- rowMeans(Avg$SCT)
   MA_data$mlgPval<- -log10(MA_data$p_val)
   xlim_max<- summary(MA_data$baseMean)[5]+100;
   p2<- ggplot(MA_data, aes(x = baseMean, y = avg_log2FC)) +
        geom_point(aes(color = cluster, size = mlgPval), show.legend = TRUE)+
        scale_radius(range = c(.1, 2)) +
        xlim(c(0,xlim_max))+theme_bw()+
        scale_color_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]))
   data<-as.data.frame(table(obj$orig.ident))
   colnames(data)<-c("Sample","cellnumber");
    p4<- ggplot(data, aes(x=Sample,y=cellnumber,fill=Sample )) +  
        geom_bar(stat = "identity" ) +
        scale_fill_manual(values =  c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]) ) +
        theme_light();
   p3<- p1+p2+p4
   print(p3)
}
dev.off()

# global 
pdf("./05_ORN_cluster2/06_stage_specific/all_cluster_cluster_DEG_pval_MAplot-raw_RNAvsSCT.pdf",width=18,height=6)
  obj<- ORN
  DefaultAssay(obj)<- "SCT";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  # p distribution 
  p_val_data<- data.frame()
  for(nrow in 1:nrow(markers)){
    p_val<- markers[nrow,]$p_val;
    filter_data <- markers[markers$p_val<p_val,]
    number<- as.data.frame(table(filter_data$cluster))
    p_val_data_subset<- data.frame(p_val,NE=number[1,2],Nurse=number[2,2],Forager=number[3,2])
    p_val_data<- rbind(p_val_data,p_val_data_subset)
  }
  p_val_data_long<- melt(p_val_data,"p_val")
  colnames(p_val_data_long)<- c("p_val","Stage","Freq")
  p1<-ggplot(p_val_data_long, aes(x=p_val, y=Freq, color=Stage,shape=Stage,labels=p_val)) + 
  geom_line()+
  geom_point(size=4)+
   theme(axis.title=element_text(size = 20),
        axis.text=element_text(size = 20),
        legend.text=element_text(size = 20),
        legend.title=element_text(size = 20),
        legend.position="top")+theme_bw()+
   ggtitle(cluster)+
  scale_color_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]))
   # MA plot
   MA_data<- markers[,c(1,2,6,7)]
   Avg <- AverageExpression(obj,features=MA_data$gene,assays = "SCT")
   MA_data<- MA_data[names(rowMeans(Avg$SCT)),]
   MA_data$baseMean<- rowMeans(Avg$SCT)
   MA_data$mlgPval<- -log10(MA_data$p_val)
   xlim_max<- summary(MA_data$baseMean)[5]+100;
   p2<- ggplot(MA_data, aes(x = baseMean, y = avg_log2FC)) +
        geom_point(aes(color = cluster, size = mlgPval), show.legend = TRUE)+
        scale_radius(range = c(.1, 2)) +ggtitle("SCT")+
        xlim(c(0,xlim_max))+theme_bw()+
        scale_color_manual(values = c(myUmapcolors[17],myUmapcolors[1],myUmapcolors[2]))
   p3<- p1+p2
   print(p3)
dev.off()

# UpsetR for NE gene in different cluster 
NE_list<- list()
Nurse_list<- list()
Forager_list<- list()

for (cluster in cluster_specific$cluster){
  print(cluster)
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "raw_RNA";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  NE_gene<- markers[markers$cluster=="NE",]$gene;
  Nurse_gene<- markers[markers$cluster=="Nurse",]$gene
  Forager_gene<- markers[markers$cluster=="Forager",]$gene
  NE_list<- c(NE_list,as.data.frame(NE_gene))
  Nurse_list<- c(Nurse_list,as.data.frame(Nurse_gene))
  Forager_list<- c(Forager_list,as.data.frame(Forager_gene))
}
  obj<- ORN
  DefaultAssay(obj)<- "raw_RNA";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  NE_gene<- markers[markers$cluster=="NE",]$gene;
  Nurse_gene<- markers[markers$cluster=="Nurse",]$gene
  Forager_gene<- markers[markers$cluster=="Forager",]$gene
  NE_list<- c(NE_list,as.data.frame(NE_gene))
  Nurse_list<- c(Nurse_list,as.data.frame(Nurse_gene))
  Forager_list<- c(Forager_list,as.data.frame(Forager_gene))
names(NE_list)<- c(cluster_specific$cluster,"global")
names(Nurse_list)<- c(cluster_specific$cluster,"global")
names(Forager_list)<- c(cluster_specific$cluster,"global")

library(UpSetR)
pdf("./05_ORN_cluster2/06_stage_specific/UpsetR_Stage_DEG_raw_RNA.pdf",width=16,height=6)
upset(fromList(NE_list),nsets=8,nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
upset(fromList(Nurse_list),nsets=8, nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
upset(fromList(Forager_list),nsets=8, nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()






