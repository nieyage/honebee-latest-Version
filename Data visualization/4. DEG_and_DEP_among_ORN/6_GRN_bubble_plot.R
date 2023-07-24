library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
# remove nopower cluster,then plot UMAP
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )
ORN<- subset(ORN,idents=setdiff(levels(ORN),"33"))
# part1: bubble for TF enrich in DEP 
tau_data<-read.csv("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_tau-0.85_cluster_specfic_data_peak_ORN.csv",row.names=1)

library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'insects', all_versions = FALSE)
)
pfm_honeybee <-  readJASPARMatrix("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/CisBP-honeybee.jaspar", matrixClass=c("PFM", "PWM", "PWMProb"))

# Merge the two PFMatrixList objects
mergedPfmList <- c(pfm, pfm_honeybee)

library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN) <- 'peaks_ORN_subcluster'
# add motif information
ORN <- AddMotifs(
  object = ORN,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = mergedPfmList
)

# motif enrich in peak  

All_motif_info <- data.frame()
for (cluster in levels(ORN)){
  obj<- subset(ORN,idents=cluster);
  print(cluster)
  cluster_peak <- tau_data[tau_data$cluster==cluster,]$peak;
  if(length(cluster_peak!=0)){
  enriched.motifs <- FindMotifs(obj,features = cluster_peak);
  enriched.motifs$cluster <- cluster;
  enriched.motifs <- enriched.motifs[enriched.motifs$pvalue<0.05,]
  All_motif_info <- rbind(All_motif_info,enriched.motifs)}
}
# motif barplot statistics
table(All_motif_info$cluster)
DefaultAssay(ORN)<- "raw_RNA"
for(nrow in 1:nrow(All_motif_info)){
	if(length(grep("^M",All_motif_info[nrow,]$motif.name))){
		All_motif_info$TF_name[nrow]<- strsplit(All_motif_info$motif.name[nrow]," ")[[1]][2];
	}else{
		if(All_motif_info[nrow,]$motif.name %in% rownames(ORN)){
			All_motif_info$TF_name[nrow]<- All_motif_info[nrow,]$motif.name
		}else{All_motif_info$TF_name[nrow]<- All_motif_info$motif.name[nrow]}
	}
	
}

fly2honeybee<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/05_fly2honeybee.csv")

for (i in which(All_motif_info$TF_name%in%fly2honeybee$fly_gene)){
	All_motif_info$TF_name[i]<- fly2honeybee[match(All_motif_info$TF_name[i],fly2honeybee$fly_gene),]$honeybee_gene_name
}

Avg <- AverageExpression(ORN,features=rownames(ORN),assays = "raw_RNA")
# add exp info 
Avg_data<- Avg$raw_RNA;
colnames(Avg_data)<- levels(ORN)
All_motif_info$exp<- NA
for(nrow in 1:nrow(All_motif_info)){
	cluster<- All_motif_info[nrow,]$cluster;
	if(length(strsplit(All_motif_info$TF_name[nrow],"::")[[1]])!=1){
		TF<- strsplit(All_motif_info$TF_name[nrow],"::")[[1]];
		if(length(which(TF%in%rownames(Avg_data)))){
		exp<- mean(Avg_data[TF,cluster])
        All_motif_info$exp[nrow]<- exp;
		}
		}else{
			TF<- All_motif_info[nrow,]$TF_name;
			if(TF%in% rownames(Avg_data)){
			exp<- Avg_data[TF,cluster]
	        All_motif_info$exp[nrow]<- exp;
			}
	        
		}
}

write.csv(All_motif_info,"./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_enrich.csv")
# get the TF name form motif(long data frame)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
motif_data<- na.omit(All_motif_info)
motif_data$cluster<- factor(motif_data$cluster,levels=levels(ORN))
motif_data$TF_name<- factor(motif_data$TF_name,levels=rev(unique(motif_data$TF_name)))
motif_data<- motif_data[order(motif_data$cluster),]

pdf("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_enrich_bubbleplot.pdf",width=10,height=12)
ggplot(motif_data,aes(x=cluster,y=TF_name))+
  geom_point(aes(size=fold.enrichment,color=log10(exp+1)),alpha=0.8)+
  scale_size(range=c(1,12))+
  scale_color_gradientn(colours=c("#080C26","#A61C68","#D92958","#F26E50","#F2E0D5"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dev.off()


 # chromvar 
 DefaultAssay(ORN)<- "peaks_ORN_subcluster"
  ORN <- RunChromVAR(object = ORN,genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  DefaultAssay(ORN) <- 'chromvar'
  differential.activity <- FindAllMarkers(object = ORN,only.pos = TRUE,mean.fxn = rowMeans,fc.name = "avg_diff")

for (nrow in 1:nrow(differential.activity)){
	gene<- differential.activity[nrow,]$gene;
	gene<- sub("-","_",gene)
	differential.activity$motif.name[nrow]<- mergedPfmList[[gene]]@name 
}
table(differential.activity$cluster)
DefaultAssay(ORN)<- "raw_RNA"
for(nrow in 1:nrow(differential.activity)){
	if(length(grep("^M",differential.activity[nrow,]$motif.name))){
		differential.activity$TF_name[nrow]<- strsplit(differential.activity$motif.name[nrow]," ")[[1]][2];
	}else{
		if(differential.activity[nrow,]$motif.name %in% rownames(ORN)){
			differential.activity$TF_name[nrow]<- differential.activity[nrow,]$motif.name
		}else{differential.activity$TF_name[nrow]<- differential.activity$motif.name[nrow]}
	}
	
}

fly2honeybee<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/05_fly2honeybee.csv")

for (i in which(differential.activity$TF_name%in%fly2honeybee$fly_gene)){
	differential.activity$TF_name[i]<- fly2honeybee[match(differential.activity$TF_name[i],fly2honeybee$fly_gene),]$honeybee_gene_name
}

Avg <- AverageExpression(ORN,features=rownames(ORN),assays = "raw_RNA")
# add exp info 
Avg_data<- Avg$raw_RNA;
colnames(Avg_data)<- levels(ORN)
differential.activity$exp<- NA
for(nrow in 1:nrow(differential.activity)){
	cluster<- differential.activity[nrow,]$cluster;
	if(length(strsplit(differential.activity$TF_name[nrow],"::")[[1]])!=1){
		TF<- strsplit(differential.activity$TF_name[nrow],"::")[[1]];
		if(length(which(TF%in%rownames(Avg_data)))){
		exp<- mean(Avg_data[TF%in%rownames(Avg_data),cluster])
        differential.activity$exp[nrow]<- exp;
		}
		}else{
			TF<- differential.activity[nrow,]$TF_name;
			if(TF%in% rownames(Avg_data)){
			exp<- Avg_data[TF,cluster]
	        differential.activity$exp[nrow]<- exp;
			}
	        
		}
}
write.csv(differential.activity,"./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_differential.activity.csv")
da_data<- na.omit(differential.activity)
da_data$cluster<- factor(da_data$cluster,levels=levels(ORN))
da_data$TF_name<- factor(da_data$TF_name,levels=rev(unique(da_data$TF_name)))
da_data<- da_data[order(da_data$cluster),]
pdf("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_differential.activity_bubbleplot.pdf",width=10,height=16)
ggplot(da_data,aes(x=cluster,y=TF_name))+
  geom_point(aes(size=avg_diff,color=log10(exp+1)),alpha=0.8)+
  scale_size(range=c(1,12))+
  scale_color_gradientn(colours=c("#080C26","#A61C68","#D92958","#F26E50","#F2E0D5"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dev.off()

# part2: bubble for TF enrich in OR promoter 





# part3: 








