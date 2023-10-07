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

# random add the fake fake_stage label
all_barcode<- colnames(ORN)
NE_barcode<- sample(all_barcode,1323)
Nurse_barcode<- sample(setdiff(all_barcode,NE_barcode),1323)
Forager_barcode<- setdiff(all_barcode,c(NE_barcode,Nurse_barcode))

fake_info<- data.frame(barcode=c(NE_barcode,Nurse_barcode,Forager_barcode),
  label=c(rep("NE",1323),rep("Nurse",1323),rep("Forager",1324)))
ORN$fake_stage<- fake_info[match(rownames(ORN@meta.data),fake_info$barcode),]$label
ORN$fake_stage<- factor(ORN$fake_stage,levels=c("NE","Nurse","Forager"))


# different Stage
# before correct
Idents(ORN)<- ORN$fake_stage
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
pdf("./11_stage_permutation/01_DEG/global_stage_numi_vs_ngene_point.pdf",width=12,height=4)
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
saveRDS(ORN_correct,"./11_stage_permutation/01_DEG/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_SCT_v2.rds")

ORN_correct<- readRDS("./11_stage_permutation/01_DEG/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_SCT_v2.rds")

# after correct
Idents(ORN_correct)<- ORN_correct$fake_stage
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
pdf("./11_stage_permutation/01_DEG/global_stage_numi_vs_ngene_point_SCTv2.pdf",width=12,height=4)
p1|p2|p3
dev.off()

ORN<- ORN_correct
Idents(ORN)<- ORN$subcluster
NE_list<- list()
Nurse_list<- list()
Forager_list<- list()
for (cluster in levels(ORN)){
  print(cluster)
  obj<- subset(ORN,idents=cluster);
  DefaultAssay(obj)<- "SCT_v2";
  Idents(obj)<- obj$fake_stage
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  NE_gene<- markers[markers$cluster=="NE",]$gene;
  Nurse_gene<- markers[markers$cluster=="Nurse",]$gene
  Forager_gene<- markers[markers$cluster=="Forager",]$gene
  NE_list<- c(NE_list,as.data.frame(NE_gene))
  Nurse_list<- c(Nurse_list,as.data.frame(Nurse_gene))
  Forager_list<- c(Forager_list,as.data.frame(Forager_gene))
}

# Stage DEG present number in ORN cluster 
NE_all<- unlist(NE_list)
Nurse_all<- unlist(Nurse_list)
Forager_all<- unlist(Forager_list)

pdf("./11_stage_permutation/01_DEG/stage_DEG_freq_permutation.pdf",width=6,height=4)
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
pdf("./11_stage_permutation/01_DEG/neuron_activated_gene_violin_dCREB_permutation.pdf",width=25,height=4)
Stacked_VlnPlot(seurat_object = ORN_correct, features = dCREB2, x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "fake_stage")
dev.off()
pdf("./11_stage_permutation/01_DEG/neuron_activated_gene_violin_Obp_permutation.pdf",width=25,height=10)
Stacked_VlnPlot(seurat_object = ORN_correct, features = OBP[1:10], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "fake_stage")
Stacked_VlnPlot(seurat_object = ORN_correct, features = OBP[11:19], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "fake_stage")
dev.off()

pdf("./11_stage_permutation/01_DEG/neuron_activated_gene_violin_Trp_permutation.pdf",width=25,height=10)
Stacked_VlnPlot(seurat_object = ORN_correct, features = Trp[1:8], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "fake_stage")
Stacked_VlnPlot(seurat_object = ORN_correct, features = Trp[9:17], x_lab_rotate = TRUE, colors_use = sample_colors, split.by = "fake_stage")
dev.off()

# DEP 

# different Stage
# before correct
Idents(ORN)<- ORN$fake_stage
NE<- subset(ORN,idents="NE")
NE<- as.matrix(NE@assays$peaks_ORN_subcluster@counts)
cell_attr<- data.frame(n_umi=colSums(NE),n_peak= colSums(NE>0))
p1<- ggplot(cell_attr,aes(n_umi,n_peak))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("NE")+xlim(0,15000)+ylim(0,4000)
Nurse<- subset(ORN,idents="Nurse")
Nurse<- as.matrix(Nurse@assays$peaks_ORN_subcluster@counts)
cell_attr<- data.frame(n_umi=colSums(Nurse),n_peak= colSums(Nurse>0))
p2<- ggplot(cell_attr,aes(n_umi,n_peak))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("Nurse")+xlim(0,15000)+ylim(0,4000)
Forager<- subset(ORN,idents="Forager")
Forager<- as.matrix(Forager@assays$peaks_ORN_subcluster@counts)
cell_attr<- data.frame(n_umi=colSums(Forager),n_peak= colSums(Forager>0))
p3<- ggplot(cell_attr,aes(n_umi,n_peak))+geom_point(alpha=0.3,shape=16)+geom_density_2d(size=0.3)+ggtitle("Forager")+xlim(0,15000)+ylim(0,4000)
pdf("./11_stage_permutation/02_DEP/global_stage_numi_vs_npeak_point.pdf",width=12,height=4)
p1|p2|p3
dev.off()

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'insects', all_versions = FALSE)
)
pfm_honeybee <- readJASPARMatrix("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/CisBP-honeybee.jaspar", matrixClass=c("PFM", "PWM", "PWMProb"))
# Merge the two PFMatrixList objects
mergedPfmList <- c(pfm, pfm_honeybee)
DefaultAssay(ORN) <- 'peaks_ORN_subcluster'
# add motif information
ORN <- AddMotifs(
  object = ORN,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = mergedPfmList
)
ORN <- RegionStats(ORN, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)

NE_list<- list()
Nurse_list<- list()
Forager_list<- list()

NE_motiflist<- list()
Nurse_motiflist<- list()
Forager_motiflist<- list()

NE_motiflist_da<- list()
Nurse_motiflist_da<- list()
Forager_motiflist_da<- list()
enriched.motifs_last<- data.frame()
differential.activity_last<- data.frame()
Idents(ORN) <- ORN$subcluster

for(cluster in levels(ORN)){
  obj<- subset(ORN,idents=cluster);
  print(cluster)
  DefaultAssay(obj)<- "peaks_ORN_subcluster";
  obj$fake_stage<- factor(obj$fake_stage,levels=c("NE","Nurse","Forager"))
  Idents(obj)<- obj$fake_stage
  obj<-ScaleData(obj,features=rownames(obj))
  # stage DEP in cluster 
  markers <- FindAllMarkers(obj,test.use = 'LR', only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.2,latent.vars = 'nCount_peaks')
  markers <- markers[which(markers$p_val<0.05),]
  write.csv(markers,paste0("./11_stage_permutation/02_DEP/Cluster",cluster,"_specific_markers.csv",sep=""))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
  p<- DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =c("#5CC8F2","#009E73","#E69F00"),disp.min = -0.5,disp.max = 0.5,size = 2,group.by = "fake_stage") + scale_fill_gradientn(colors = c("white", "firebrick3"))+NoLegend()
  pdf(paste0("./11_stage_permutation/02_DEP/Cluster",cluster,"_DEP_specific_heatmap.pdf",sep=""))
  print(p)
  dev.off()
  # validation
  # true DEP or false DEP
  top3<- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC);
  obj_feature<-top3$gene
  regular_df<- data.frame()
  for(i in obj_feature){
    tmp<- data.frame(cell_log_umi=log10(obj@meta.data$nCount_peaks_ORN_subcluster+1),
      scale.residual=obj@assays$peaks_ORN_subcluster@scale.data[rownames(obj@assays$peaks_ORN_subcluster@scale.data)==i,drop=TRUE],
      peak=rep(i,length(log10(obj@meta.data$nCount_peaks_ORN_subcluster+1)))
      )
    tmp$stage<- obj$fake_stage
    regular_df<- rbind(regular_df,tmp)
  }
  regular_df$stage<- factor(regular_df$stage,levels=c("NE","Nurse","Forager"))
  p2<- ggplot(regular_df,aes(cell_log_umi,scale.residual,color=stage))+
  geom_point(alpha=0.8,shape=16,size=1)+
  #geom_density_2d(size=0.3,color="cadetblue1")+
  facet_wrap(.~factor(peak,levels=obj_feature),ncol =3)+
  scale_color_manual(values = c("#5CC8F2","#009E73","#E69F00") )
  pdf(paste0("./11_stage_permutation/02_DEP/Cluster",cluster,"_DEP_validation.pdf",sep=""),width=6,height=6)
  print(p2)
  dev.off()
  # get the list
  NE_peak<- markers[markers$cluster=="NE",]$gene;
  Nurse_peak<- markers[markers$cluster=="Nurse",]$gene
  Forager_peak<- markers[markers$cluster=="Forager",]$gene
  NE_list<- c(NE_list,as.data.frame(NE_peak))
  Nurse_list<- c(Nurse_list,as.data.frame(Nurse_peak))
  Forager_list<- c(Forager_list,as.data.frame(Forager_peak))
  # motif enrichment 
  DefaultAssay(obj)<- "peaks_ORN_subcluster"
  NE_enriched.motifs<- c()
  Nurse_enriched.motifs<- c()
  Forager_enriched.motifs<- c()
  if(length(NE_peak)>0){
  NE_enriched.motifs <- FindMotifs(object = obj,features = NE_peak)
  NE_enriched.motifs <- NE_enriched.motifs[NE_enriched.motifs$pvalue<0.05,]
  NE_motiflist<- c(NE_motiflist,NE_enriched.motifs$motif.name)
    NE_enriched.motifs$cluster<- "NE"

  }
if(length(Nurse_peak)>0){
  Nurse_enriched.motifs <- FindMotifs(object = obj,features = Nurse_peak)
  Nurse_enriched.motifs <- Nurse_enriched.motifs[Nurse_enriched.motifs$pvalue<0.05,]
  Nurse_motiflist<- c(Nurse_motiflist,Nurse_enriched.motifs$motif.name)
    Nurse_enriched.motifs$cluster<- "Nurse"

  }
  if(length(Forager_peak)>0){
  Forager_enriched.motifs <- FindMotifs(object = obj,features = Forager_peak)
  Forager_enriched.motifs <- Forager_enriched.motifs[Forager_enriched.motifs$pvalue<0.05,]
  Forager_motiflist<- c(Forager_motiflist,Forager_enriched.motifs$motif.name)
    Forager_enriched.motifs$cluster<- "Forager"
  }
 obj <- RunChromVAR(object = obj,genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  DefaultAssay(obj) <- 'chromvar'
  differential.activity <- FindAllMarkers(object = obj,only.pos = TRUE,mean.fxn = rowMeans,fc.name = "avg_diff")
  # save the results
  # TF (motif enrichment)
  enriched.motifs<- data.frame()
  enriched.motifs<- rbind(NE_enriched.motifs,Nurse_enriched.motifs,Forager_enriched.motifs)
  if(nrow(enriched.motifs)>0){
  enriched.motifs$subcluster<- cluster;
  enriched.motifs_last<- rbind(enriched.motifs_last,enriched.motifs)
  }
  if(nrow(differential.activity)!=0){
    differential.activity$subcluster<- cluster;
    differential.activity_last<- rbind(differential.activity_last,differential.activity)
  }
}
}
write.csv(enriched.motifs_last,"./11_stage_permutation/03_TF/enriched.motifs_results.csv")
write.csv(differential.activity_last,"./11_stage_permutation/03_TF/differential.activity_results.csv")

# Stage DEP present number in ORN cluster 
NE_all<- unlist(NE_list)
Nurse_all<- unlist(Nurse_list)
Forager_all<- unlist(Forager_list)
pdf("./11_stage_permutation/02_DEP/stage_DEP_freq.pdf",width=6,height=4)
cluster_number<-as.data.frame(table(table(NE_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of DEP") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of NE-specific peak occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Nurse_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of DEP") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Nurse-specific peak occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Forager_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of DEP") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Forager-specific peak occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
dev.off();

# Stage TF(motif enrichment) present number in ORN cluster 
NE_all<- unlist(NE_motiflist)
Nurse_all<- unlist(Nurse_motiflist)
Forager_all<- unlist(Forager_motiflist)

pdf("./11_stage_permutation/03_TF/stage_TF_freq.pdf",width=6,height=4)
cluster_number<-as.data.frame(table(table(NE_all)))
#cluster_number$Var1<-as.character(cluster_number$Var1)
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of TF") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of NE-specific TF occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Nurse_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of TF") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Nurse-specific TF occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Forager_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of TF") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Forager-specific TF occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
dev.off();

# Stage DEG present number in ORN cluster 
NE_all<- unlist(NE_motiflist_da)
Nurse_all<- unlist(Nurse_motiflist_da)
Forager_all<- unlist(Forager_motiflist_da)

pdf("./11_stage_permutation/03_TF/stage_TF_chromVar_freq.pdf",width=6,height=4)
cluster_number<-as.data.frame(table(table(NE_all)))
#cluster_number$Var1<-as.character(cluster_number$Var1)
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of TF") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of NE-specific TF occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Nurse_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of TF") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Nurse-specific TF occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
cluster_number<-as.data.frame(table(table(Forager_all)))
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "Freq")) +  
        xlab("number of cluster") +
        ylab("number of TF") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_point() +
        ggtitle("Frequency of Forager-specific TF occurrence")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
dev.off();

# nothing at all 

# the TF frep high in Forager
M02961-2.00 LOC551896 M06430-2.00 LOC409246     MA0213.1 MA0531.1
                   25                    23           28       28


> names(table(Forager_all)[which(table(Forager_all)>6)])
TF<- c("LOC551896","LOC409246","brk","CTCF")
fly2honeybee<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/05_fly2honeybee.csv")
TF[3:4]<- fly2honeybee[match(TF[3:4],fly2honeybee$fly_gene),]$honeybee_gene_name
ORN_correct<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_SCT_v2.rds")
library(scCustomize)
pdf("./05_ORN_cluster2/06_stage_specific/DEP/Forager_da_shared_TF_exp_violin.pdf",width=25,height=8)
Stacked_VlnPlot(seurat_object = ORN_correct, features = TF, x_lab_rotate = TRUE, split.by = "fake_stage")
dev.off()











