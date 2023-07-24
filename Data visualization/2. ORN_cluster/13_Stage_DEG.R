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

# different Stage
# before correct
Idents(ORN)<- ORN$orig.ident
NE<- subset(ORN,idents="NE")
NE<- as.matrix(NE@assays$raw_RNA@counts)
cell_attr<- data.frame(n_umi=colSums(NE),n_gene= colSums(NE>0))
p1<- ggplot(cell_attr,aes(n_umi,n_gene))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("NE")+xlim(0,15000)+ylim(0,4000)
Nurse<- subset(ORN,idents="Nurse")
Nurse<- as.matrix(Nurse@assays$raw_RNA@counts)
cell_attr<- data.frame(n_umi=colSums(Nurse),n_gene= colSums(Nurse>0))
p2<- ggplot(cell_attr,aes(n_umi,n_gene))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("Nurse")+xlim(0,15000)+ylim(0,4000)
Forager<- subset(ORN,idents="Forager")
Forager<- as.matrix(Forager@assays$raw_RNA@counts)
cell_attr<- data.frame(n_umi=colSums(Forager),n_gene= colSums(Forager>0))
p3<- ggplot(cell_attr,aes(n_umi,n_gene))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("Forager")+xlim(0,15000)+ylim(0,4000)
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/global_stage_numi_vs_ngene_point.pdf",width=12,height=4)
p1|p2|p3
dev.off()

#v2 regularization
ORN_correct <- Seurat::SCTransform(ORN,
	assay="RNA",
	new.assay.name = "SCT_v2",
	ncells = 4000,
    residual.features = NULL,# all gene
    variable.features.n = 3000,
    do.scale = FALSE,
    do.center = TRUE,
    return.only.var.genes = FALSE,
    seed.use = 1448145,
    verbose = TRUE,
	vst.flavor="v2")
saveRDS(ORN_correct,"./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_SCT_v2.rds")

ORN_correct<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_SCT_v2.rds")

# after correct
Idents(ORN_correct)<- ORN_correct$orig.ident
NE<- subset(ORN_correct,idents="NE")
NE<- as.matrix(NE@assays$SCT_v2@counts)
cell_attr<- data.frame(n_umi=colSums(NE),n_gene= colSums(NE>0))
p1<- ggplot(cell_attr,aes(n_umi,n_gene))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("NE")+xlim(2000,3600)+ylim(0,2500)
Nurse<- subset(ORN_correct,idents="Nurse")
Nurse<- as.matrix(Nurse@assays$SCT_v2@counts)
cell_attr<- data.frame(n_umi=colSums(Nurse),n_gene= colSums(Nurse>0))
p2<- ggplot(cell_attr,aes(n_umi,n_gene))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("Nurse")+xlim(1500,3600)+ylim(0,2500)
Forager<- subset(ORN_correct,idents="Forager")
Forager<- as.matrix(Forager@assays$SCT_v2@counts)
cell_attr<- data.frame(n_umi=colSums(Forager),n_gene= colSums(Forager>0))
p3<- ggplot(cell_attr,aes(n_umi,n_gene))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("Forager")+xlim(1500,3600)+ylim(0,2500)
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/global_stage_numi_vs_ngene_point_SCTv2.pdf",width=12,height=4)
p1|p2|p3
dev.off()

# OR GENE IN THE CLUSTER
library(scCustomize)
cluster_specific<- read.csv("./05_ORN_cluster2/06_stage_specific/cluster_specific.csv")

specfic <- cluster_specific$cluster
# SCT 
Idents(ORN_correct)<- ORN_correct$subcluster
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_specific_OR_exp_SCTv2_counts.pdf",width=5,height=10)
for(cluster in specfic){
	obj<-subset(ORN_correct,idents=cluster);
	obj_feature<- strsplit(cluster_specific[cluster_specific$cluster==cluster,]$OR_specfic,",")[[1]]
	p1<-Stacked_VlnPlot(seurat_object = obj,assay = "SCT_v2",slot="counts",features = obj_feature,group.by ="orig.ident", x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
	regular_df<- data.frame()
	for(i in obj_feature){
		tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
			scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
			gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
			)
		tmp$stage<- obj$orig.ident
		regular_df<- rbind(regular_df,tmp)
	}
	regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
	p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
	geom_point(alpha=0.7,shape=16,size=1)+
	geom_density_2d(size=0.7,color="cadetblue1")+
	facet_wrap(.~factor(gene))+
    scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )+ggtitle(cluster)
	p3<- p1/p2
	print(p3)
}
dev.off()

# OR exp in their cluster among stages
ORN_correct$orig.ident<- factor(ORN_correct$orig.ident,levels=c("NE","Nurse","Forager"))
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_OR_exp_raw_RNA_data_modify.pdf",width=5,height=8)
cluster<- "p4:0_0";
obj_feature<-"Or57"
obj<-subset(ORN_correct,idents=cluster);
	p1<-Stacked_VlnPlot(seurat_object = obj,assay = "SCT_v2",slot="counts",features = obj_feature,group.by ="orig.ident", x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
	regular_df<- data.frame()
	for(i in obj_feature){
		tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
			scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
			gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
			)
		tmp$stage<- obj$orig.ident
		regular_df<- rbind(regular_df,tmp)
	}
	regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
	p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
	geom_point(alpha=0.7,shape=16,size=1)+
	geom_density_2d(size=0.7,color="cadetblue1")+
	facet_wrap(.~factor(gene))+
    scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )+ggtitle(cluster)
	p3<- p1/p2
	print(p3)
cluster<- "p4:0_2";
obj_feature<-"Or55"
obj<-subset(ORN_correct,idents=cluster);
	p1<-Stacked_VlnPlot(seurat_object = obj,assay = "SCT_v2",slot="counts",features = obj_feature,group.by ="orig.ident", x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
	regular_df<- data.frame()
	for(i in obj_feature){
		tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
			scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
			gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
			)
		tmp$stage<- obj$orig.ident
		regular_df<- rbind(regular_df,tmp)
	}
	regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
	p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
	geom_point(alpha=0.7,shape=16,size=1)+
	geom_density_2d(size=0.7,color="cadetblue1")+
	facet_wrap(.~factor(gene))+
    scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )+ggtitle(cluster)
	p3<- p1/p2
	print(p3)
cluster<- "p2:17";
obj_feature<-"Or109"
obj<-subset(ORN_correct,idents=cluster);
	p1<-Stacked_VlnPlot(seurat_object = obj,assay = "SCT_v2",slot="counts",features = obj_feature,group.by ="orig.ident", x_lab_rotate = TRUE,plot_spacing = 0.3, plot_legend = T,colors_use = c("#5CC8F2","#009E73","#E69F00"))
	regular_df<- data.frame()
	for(i in obj_feature){
		tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
			scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
			gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
			)
		tmp$stage<- obj$orig.ident
		regular_df<- rbind(regular_df,tmp)
	}
	regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
	p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
	geom_point(alpha=0.7,shape=16,size=1)+
	geom_density_2d(size=0.7,color="cadetblue1")+
	facet_wrap(.~factor(gene))+
    scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )+ggtitle(cluster)
	p3<- p1/p2
	print(p3)
dev.off()

#sample cell number in obj 
Idents(ORN)<- ORN$subcluster
ORN$orig.ident<- factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
pdf("./05_ORN_cluster2/06_stage_specific/cluster_specific_cellnumber.pdf",width=4,height=3)
for(cluster in specfic){
	obj<-subset(ORN,idents=cluster);
	data<-as.data.frame(table(obj$orig.ident))
    colnames(data)<-c("Sample","cellnumber");
    p<- ggplot(data, aes(x=Sample,y=cellnumber,fill=Sample )) +  
        geom_bar(stat = "identity" ) +
        ggtitle(cluster)+
        scale_fill_manual(values = c("#5CC8F2","#009E73","#E69F00") ) +
        theme_light();
    print(p)
}
dev.off()

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

# global Stage DEG 
obj<- ORN_correct
  DefaultAssay(obj)<- "SCT_v2";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  markers <- markers[which(markers$p_val<0.05),]
  write.csv(markers,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/global_specific_markers_FC0.5.csv",sep=""))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
  #obj<-ScaleData(obj,features=rownames(obj))
  p<- DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =c("#E69F00","#55B4E9","#009E73"),disp.min = -0.5,disp.max = 0.5,size = 2,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster_specific_heatmap_FC0.5.pdf",sep=""),width=8,height=16)
  print(p)
  dev.off()

table(markers$cluster)

  gene_in_cluster<-markers[markers$cluster=="NE",]$gene
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  NE_ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
  gene_in_cluster<-markers[markers$cluster=="Nurse",]$gene
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  Nurse_ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
  gene_in_cluster<-markers[markers$cluster=="Forager",]$gene
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  Forager_ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
NE_ego<- as.data.frame(NE_ego)
NE_ego$stage<-"NE"
Nurse_ego<- as.data.frame(Nurse_ego)
Nurse_ego$stage<-"Nurse"
Forager_ego<- as.data.frame(Forager_ego)
Forager_ego$stage<-"Forager"

all_ego<-rbind(NE_ego,Nurse_ego,Forager_ego)
write.csv(all_ego,"./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Global_Stage_specific_markers_GO.csv")

top5_term <- all_ego %>% group_by(stage) %>% top_n(n = 8, wt = pvalue);
top5_term$stage<-factor(top5_term$stage,levels=c("NE","Nurse","Forager"))
library(ggplot2)
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Global_Stage_specific_markers_GO.pdf",width=16,height=10)
p <- ggplot(top5_term,aes(y=Count,x=Description,fill=pvalue)) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(stage~.,scales = "free",space = "free") + 
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
# true DEG or false DEG 
top5<- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC);
obj_feature<-top5$gene
regular_df<- data.frame()
for(i in obj_feature){
	tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
		scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
		gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
		)
	tmp$stage<- obj$orig.ident
	regular_df<- rbind(regular_df,tmp)
}
regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
geom_point(alpha=0.3,shape=16,size=0.7)+
geom_density_2d(size=0.3,color="cadetblue1")+
facet_wrap(.~factor(gene,levels=obj_feature),ncol =5)+
   scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/global_Stage_DEG_top5_validation.pdf",width=15,height=9)
p2
dev.off()

# stage DEG in cluster 
ORN<- ORN_correct
for(cluster in specfic){
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "SCT_v2";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
  markers <- markers[which(markers$p_val<0.05),]
  write.csv(markers,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,"_specific_markers.csv",sep=""))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
  #obj<-ScaleData(obj,features=rownames(obj))
  p<- DoHeatmap(object = obj,features=top10$gene,label=T, group.colors =c("#55B4E9","#009E73","#E69F00"),disp.min = -0.5,disp.max = 0.5,size = 2,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,"_specific_heatmap.pdf",sep=""))
  print(p)
  dev.off()
  # Go and KEGG 
for (i in unique(markers$cluster)){
  print(i);
  gene_in_cluster<-markers[markers$cluster==i,7]
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  #GO
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  # validation
  # true DEG or false DEG 
  top3<- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC);
  obj_feature<-top3$gene
  regular_df<- data.frame()
  for(i in obj_feature){
  	tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
  		scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
  		gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
  		)
  	tmp$stage<- obj$orig.ident
  	regular_df<- rbind(regular_df,tmp)
  }
  regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
  p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
  geom_point(alpha=0.8,shape=16,size=1)+
  geom_density_2d(size=0.3,color="cadetblue1")+
  facet_wrap(.~factor(gene,levels=obj_feature),ncol =3)+
  scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,"_DEG_validation.pdf",sep=""),width=6,height=6)
  print(p2)
  dev.off()

}
}

# UpsetR for NE gene in different cluster 
NE_list<- list()
Nurse_list<- list()
Forager_list<- list()

for (cluster in cluster_specific$cluster){
  print(cluster)
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "SCT_v2";
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
  DefaultAssay(obj)<- "SCT_v2";
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
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/UpsetR_Stage_DEG_raw_RNA.pdf",width=16,height=6)
upset(fromList(NE_list),nsets=8,nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
upset(fromList(Nurse_list),nsets=8, nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
upset(fromList(Forager_list),nsets=8, nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()


# upSetR plot with bulk RNA seq
library(DESeq2)
raw_count<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/05_ORN_cluster2/06_stage_specific/UMI_corrected/rawcount_Worker_filter.csv",row.names=1)
rownames(raw_count)<- raw_count[,1]
raw_count<- raw_count[,2:10]

condition <- factor(c(rep("other",3),rep("other",3),rep("Forager",3)))
colData <- data.frame(row.names=colnames(raw_count),condition=condition)
colData
countData <- raw_count[apply(raw_count, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <-results(dds,contrast = c("condition","Forager","other"))
Forager_DEG<- rownames(subset(res, padj < 0.05 & log2FoldChange >1))

condition <- factor(c(rep("other",3),rep("Nurse",3),rep("other",3)))
colData <- data.frame(row.names=colnames(raw_count),condition=condition)
colData
countData <- raw_count[apply(raw_count, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <-results(dds,contrast = c("condition","Nurse","other"))
Nurse_DEG<- rownames(subset(res, padj < 0.05 & log2FoldChange >1))

condition <- factor(c(rep("NE",3),rep("other",3),rep("other",3)))
colData <- data.frame(row.names=colnames(raw_count),condition=condition)
colData
countData <- raw_count[apply(raw_count, 1, sum) > 1 , ] 
countData <- na.omit(countData)
dds <- DESeqDataSetFromMatrix(countData,colData,formula(~condition)) 
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <-results(dds,contrast = c("condition","NE","other"))
NE_DEG<- rownames(subset(res, padj < 0.05 & log2FoldChange >1))

  obj<- ORN
  DefaultAssay(obj)<- "SCT_v2";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
  NE_gene<- markers[markers$cluster=="NE",]$gene;
  Nurse_gene<- markers[markers$cluster=="Nurse",]$gene
  Forager_gene<- markers[markers$cluster=="Forager",]$gene

  stage_DEG_list<- c(
    bulk_NE<-as.data.frame(NE_DEG),
    scRNA_NE<- as.data.frame(NE_gene),
    bulk_Nurse<-as.data.frame(Nurse_DEG),
    scRNA_Nurse<- as.data.frame(Nurse_gene),
    bulk_Forager<-as.data.frame(Forager_DEG),
    scRNA_Forager<- as.data.frame(Forager_gene)
    )
  names(stage_DEG_list)<-c("bulk_NE","scRNA_NE","bulk_Nurse","scRNA_Nurse","bulk_Forager","scRNA_Forager")
 library(UpSetR)
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/bulk_scRNA_UpsetR_Stage_DEG.pdf",width=10,height=6)
upset(fromList(stage_DEG_list),nsets=8,nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()

DefaultAssay(ORN)<-"raw_RNA"
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Featureplot.pdf",width=8,height=8)
  p1<-FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("Mrjp1") ,order=TRUE, ncol = 1)
  p3<-FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("Obp17") ,order=TRUE, ncol = 1)
  p4<-FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("Mrjp7") ,order=TRUE, ncol = 1)  
  p2<-FeaturePlot(ORN, reduction = 'tsne.rna',max.cutoff = 10,features = c("Mrjp4") ,order=TRUE, ncol = 1)
  f1<-p1|p2
  f2<-p3|p4
  f1/f2
dev.off()



# stage DEG in cluster P4:9 and p4:12(Or151,and Or152)
ORN<- ORN_correct
for(cluster in c("p4:9","p4:12")){
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "SCT_v2";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
  markers <- markers[which(markers$p_val<0.05),]
  write.csv(markers,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,"_specific_markers.csv",sep=""))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
  #obj<-ScaleData(obj,features=rownames(obj))
  p<- DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =c("#55B4E9","#009E73","#E69F00"),disp.min = -0.5,disp.max = 0.5,size = 2,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,"_specific_heatmap.pdf",sep=""))
  print(p)
  dev.off()
  # Go and KEGG 
for (i in unique(markers$cluster)){
  print(i);
  gene_in_cluster<-markers[markers$cluster==i,7]
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  #GO
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,i,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,i,"_GO_",type,".csv"))
      }
  }
  dev.off()
  # validation
  # true DEG or false DEG 
  top3<- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC);
  obj_feature<-top3$gene
  regular_df<- data.frame()
  for(i in obj_feature){
    tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_SCT_v2+1),
      scale.residual=obj@assays$SCT_v2@scale.data[rownames(obj@assays$SCT_v2@scale.data)==i,drop=TRUE],
      gene=rep(i,length(log10(obj@meta.data$nCount_SCT_v2+1)))
      )
    tmp$stage<- obj$orig.ident
    regular_df<- rbind(regular_df,tmp)
  }
  regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
  p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
  geom_point(alpha=0.8,shape=16,size=1)+
  geom_density_2d(size=0.3,color="cadetblue1")+
  facet_wrap(.~factor(gene,levels=obj_feature),ncol =3)+
  scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Cluster",cluster,"_DEG_validation.pdf",sep=""),width=6,height=6)
  print(p2)
  dev.off()

}
}

# merge all ORN cluster DEG for GO and KEGG 
# compare with bulk DEG 
NE_list<- list()
Nurse_list<- list()
Forager_list<- list()
for (cluster in levels(ORN)){
  print(cluster)
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "SCT_v2";
  Idents(obj)<- obj$orig.ident
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  NE_gene<- markers[markers$cluster=="NE",]$gene;
  Nurse_gene<- markers[markers$cluster=="Nurse",]$gene
  Forager_gene<- markers[markers$cluster=="Forager",]$gene
  NE_list<- c(NE_list,as.data.frame(NE_gene))
  Nurse_list<- c(Nurse_list,as.data.frame(Nurse_gene))
  Forager_list<- c(Forager_list,as.data.frame(Forager_gene))
}
NE_all<- unique(unlist(NE_list))
Nurse_all<- unique(unlist(Nurse_list))
Forager_all<- unique(unlist(Forager_list))

  gene_in_cluster<-Forager_all
  stage<-"Forager"
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  #GO
  pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/global_DEG/",stage,"_GO.pdf",sep=""))
  for (type in c("BP","CC","MF")){
    print(type);
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
    data<-as.data.frame(ego)
    if(!nrow(data)==0){
      p<-barplot(ego, showCategory=20);
      print(p);
      write.csv(ego,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/global_DEG/Cluster",stage,"_GO_",type,".csv"))
      }
  }
  dev.off()
# compare with bulk DEG 
  stage_DEG_list<- c(
    bulk_NE<-as.data.frame(NE_DEG),
    scRNA_NE<- as.data.frame(NE_all),
    bulk_Nurse<-as.data.frame(Nurse_DEG),
    scRNA_Nurse<- as.data.frame(Nurse_all),
    bulk_Forager<-as.data.frame(Forager_DEG),
    scRNA_Forager<- as.data.frame(Forager_all)
    )
  names(stage_DEG_list)<-c("bulk_NE","scRNA_NE","bulk_Nurse","scRNA_Nurse","bulk_Forager","scRNA_Forager")
 library(UpSetR)
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/global_DEG/bulk_scRNA-merge_DEG-UpsetR_Stage_DEG.pdf",width=10,height=6)
upset(fromList(stage_DEG_list),nsets=8,nintersects = 100, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
dev.off()


# Stage DEG present number in ORN cluster 
NE_all<- unlist(NE_list)
Nurse_all<- unlist(Nurse_list)
Forager_all<- unlist(Forager_list)

pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/stage_DEG_freq.pdf",width=6,height=4)
cluster_number<-as.data.frame(table(table(NE_all)))
#cluster_number$Var1<-as.character(cluster_number$Var1)
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of DEG") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of NE-specific gene occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Nurse_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of DEG") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Nurse-specific gene occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Forager_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of DEG") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Forager-specific gene occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
dev.off();

# SNMP: 感觉神经元膜蛋白：SNMP是一类与果蝇嗅觉感受器神经元中的Orco和Or接收受体的蛋白质，它对动力化合物的作用起到重要的调节作用。
# OBP and CSP;
# 离子通道相关；
# dCREB2 (cAMP response element-binding protein 2)：编码转录因子，参与嗅觉学习和记忆的形成过程。
# Trp (Transient receptor potential) 基因家族：这是一类与感知环境刺激（如温度、机械和化学刺激）相关的离子通道基因。其中的一些基因在果蝇嗅觉神经元中表达，参与感知气味的过程。

Trp<- c("LOC410823","LOC409631","LOC412391","LOC724608","LOC107964339","LOC107964338","LOC726193","LOC726724","LOC726423","HsTRPA","LOC413265","LOC100578683","LOC551894","LOC552810","LOC408777","LOC552792","LOC726119")
OBP <- grep("Obp",rownames(ORN),value=T)
dCREB2<- c("LOC409401","LOC100576912","LOC726280")

# violin plot split by group 
DefaultAssay(ORN_correct)<- "SCT_v2"
library(scCustomize)
sample_colors <- c("#5CC8F2","#009E73","#E69F00")
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/neuron_activated_gene_violin_dCREB.pdf",width=25,height=4)
Stacked_VlnPlot(seurat_object = ORN_correct, features = dCREB2, x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
dev.off()
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/neuron_activated_gene_violin_Obp.pdf",width=25,height=10)
Stacked_VlnPlot(seurat_object = ORN_correct, features = OBP[1:10], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
Stacked_VlnPlot(seurat_object = ORN_correct, features = OBP[11:19], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
dev.off()

pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/neuron_activated_gene_violin_Trp.pdf",width=25,height=10)
Stacked_VlnPlot(seurat_object = ORN_correct, features = Trp[1:8], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
Stacked_VlnPlot(seurat_object = ORN_correct, features = Trp[9:17], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
dev.off()

# unpg,fkh和Gr21a fkh和Gr63a
TF<- c("")
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/neuron_activated_gene_violin_Trp.pdf",width=25,height=10)
Stacked_VlnPlot(seurat_object = ORN_correct, features = Trp[1:8], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
Stacked_VlnPlot(seurat_object = ORN_correct, features = Trp[9:17], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "orig.ident")
dev.off()





















LOC410823          GB55618, TRPgama, trpgamma  transient receptor potential-gamma protein  transient receptor potential-gamma protein|TRPgamma cation channel    LG2 NC_037639.1 13897382  13976504  minus 21    
LOC409631          GB47942, TRPL transient-receptor-potential-like protein transient-receptor-potential-like protein   LG8 NC_037645.1 4151916 4161254 plus  18    
LOC412391       GB41177, Trpm3  transient receptor potential cation channel, subfamily M, member 3  transient receptor potential cation channel trpm    LG8 NC_037645.1 1766842 1778598 minus 21    
LOC724608      GB41297, TRP  transient receptor potential  transient receptor potential protein    LG5 NC_037642.1 12749226  12758549  plus  18    
LOC107964339    short transient receptor potential channel 6-like short transient receptor potential channel 6-like   LG4 NC_037641.1 12360139  12361936  plus  5   
LOC107964338    transient receptor potential-gamma protein  transient receptor potential-gamma protein    LG4 NC_037641.1 12358704  12360112  plus  2   
LOC726193       GB55067 short transient receptor potential channel 5-like short transient receptor potential channel 5-like   LG12  NC_037649.1 2559279 2589858 minus 12    
LOC726724     GB53163, Pyx2, TRPA5  transient receptor potential channel pyrexia  transient receptor potential channel pyrexia|transient receptor potential cation channel, subfamily A, member 5|transient receptor potential channel pyrexia 2    LG4 NC_037641.1 3132130 3136392 minus 9   
LOC726423    GB46068, Pain transient receptor potential cation channel protein painless  transient receptor potential cation channel protein painless|painless   LG6 NC_037643.1 13287832  13292359  plus  3   
HsTRPA          GB16385, GB50806  hymenoptera-specific transient receptor potential cation channel, subfamily A     LG2             
LOC413265          AmHsTRPA, GB50808 uncharacterized LOC413265 uncharacterized protein LOC413265|Hymenoptera-specific transient receptor potential channel pyrexia-like|transient receptor potential channel A   LG2 NC_037639.1 6244132 6260685 plus  9   
LOC100578683  GB41937 transient receptor potential cation channel subfamily V member 5-like                   
LOC551894    GB49284, TRPML  transient receptor potential mucolipin  mucolipin-3   LG1 NC_037638.1 23263599  23266624  minus 4   
LOC552810    GB45664 uncharacterized LOC552810 LOW QUALITY PROTEIN: uncharacterized protein LOC552810|protein FAM161A|transient receptor potential channel pyrexia   LG5 NC_037642.1 13465261  13473424  minus 13    
LOC408777    GB44092, NompC  no mechanoreceptor potential C  serine/threonine-protein phosphatase 6 regulatory ankyrin repeat subunit A    LG4 NC_037641.1 12345474  12359114  plus  17    
LOC552792    GB44065, Nan  nanchung  transient receptor potential cation channel subfamily V member 5    LG4 NC_037641.1 11948517  11953566  plus  14    
LOC726119    GB40818, Iav  inactive  transient receptor potential cation channel subfamily V member 5|uncharacterized protein LOC726119    LG1 NC_037638.1 6178108 6227750 minus 17    












