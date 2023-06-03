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
# p4:0_0
cluster<- "p4:0_0"
obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "raw_RNA";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
  markers <- markers[which(markers$p_val<0.05),]
  write.csv(markers,paste0("./05_ORN_cluster2/06_stage_specific/Cluster",cluster,"_specific_markers.csv",sep=""))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
  obj<-ScaleData(obj,features=rownames(obj))

# p distribution 
