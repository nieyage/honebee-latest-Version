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
			exp<- mean(Avg_data[TF[which(TF%in%rownames(Avg_data))],cluster])
		  exp_in_cluster<- exp
		  exp_other <- (sum(Avg_data[TF[which(TF%in%rownames(Avg_data))],])-exp_in_cluster*length(TF))/58
      All_motif_info$FC[nrow]<- exp_in_cluster/exp_other;
      All_motif_info$exp[nrow]<- exp;
		}
		}else{
			TF<- All_motif_info[nrow,]$TF_name;
			if(TF%in% rownames(Avg_data)){
			exp_in_cluster<- Avg_data[TF,cluster]
			exp_other <- (sum(Avg_data[TF,])-exp_in_cluster)/58
	    All_motif_info$exp[nrow]<- exp;
	    All_motif_info$FC[nrow]<- exp_in_cluster/exp_other;
			}
	        
		}
}


write.csv(All_motif_info,"./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_enrich.csv")

All_motif_info<- read.csv("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_enrich.csv")
# get the TF name form motif(long data frame)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
motif_data<- na.omit(All_motif_info)
motif_data$cluster<- factor(motif_data$cluster,levels=levels(ORN))
motif_data$TF_name<- factor(motif_data$TF_name,levels=rev(unique(motif_data$TF_name)))
motif_data<- motif_data[order(motif_data$cluster),]

# read the motif TF family info 
data<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/motifs-v10nr_clust-nr.flybase-m0.001-o0.0-trans2honeybee.tbl",sep="\t")

TF_list<- unique(motif_data$TF_name)

data<- data[which(data$gene_name%in%TF_list), ]

pdf("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_enrich_bubbleplot-expFC.pdf",width=10,height=12)
ggplot(motif_data,aes(x=cluster,y=TF_name))+
  geom_point(aes(size=fold.enrichment,color=log2(FC+1)),alpha=0.8)+
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
differential.activity$FC<- NA
for(nrow in 1:nrow(differential.activity)){
	cluster<- differential.activity[nrow,]$cluster;
	if(length(strsplit(differential.activity$TF_name[nrow],"::")[[1]])!=1){
		TF<- strsplit(differential.activity$TF_name[nrow],"::")[[1]];
		if(length(which(TF%in%rownames(Avg_data)))){
		  exp<- mean(Avg_data[TF[which(TF%in%rownames(Avg_data))],cluster])
		  exp_in_cluster<- exp
		  exp_other <- (sum(Avg_data[TF[which(TF%in%rownames(Avg_data))],])-exp_in_cluster*length(TF))/58
      differential.activity$FC[nrow]<- exp_in_cluster/exp_other;
      differential.activity$exp[nrow]<- exp;
		}
		}else{
			TF<- differential.activity[nrow,]$TF_name;
			if(TF%in% rownames(Avg_data)){
			exp<- Avg_data[TF,cluster]
			exp_in_cluster<- Avg_data[TF,cluster]
			exp_other <- (sum(Avg_data[TF,])-exp_in_cluster)/58
	    differential.activity$exp[nrow]<- exp;
	    differential.activity$FC[nrow]<- exp_in_cluster/exp_other;
			}
	        
		}
}
write.csv(differential.activity,"./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_differential.activity.csv")

differential.activity<- read.csv("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_differential.activity.csv")
da_data<- na.omit(differential.activity)
da_data$cluster<- factor(da_data$cluster,levels=levels(ORN))
da_data$TF_name<- factor(da_data$TF_name,levels=rev(unique(da_data$TF_name)))
da_data<- da_data[order(da_data$cluster),]
pdf("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_motif_differential.activity_bubbleplot-FC.pdf",width=10,height=16)
ggplot(da_data,aes(x=cluster,y=TF_name))+
  geom_point(aes(size=avg_diff,color=log10(FC+0.001)),alpha=0.8)+
  scale_size(range=c(1,12))+
  scale_color_gradientn(colours=c("#080C26","#A61C68","#D92958","#F26E50","#F2E0D5"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dev.off()

# part2: bubble for TF enrich in OR promoter 
cp /md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/*fa .

cat CisBP-honeybee.jaspar JASPAR2022_CORE_insects_non-redundant_pfms_jaspar.txt > merged_honeybee_pfam.jaspar
ls *fa > OR_promoter.list
for file1 in $(<OR_promoter.list)
do
  TOBIAS TFBScan --motifs merged_honeybee_pfam.jaspar --naming id --fasta ./promoter_fa/$file1 --cores 16 --outfile ./TOBIAS_out/$file1.txt
done

rename _promoter.fa.txt .txt *txt

for i in *.txt
do
echo
awk -F " " '{print FILENAME,"\t",$0}' $i >> TOBIAS_merge.txt
done


# filter TF by the TFBS 
TOBIAS_merge<- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_promoter_TOBIAS/TOBIAS_out/TOBIAS_merge.txt")
TOBIAS_merge$OR<- unlist(strsplit(TOBIAS_merge$V1,split="\\."))[seq(1,length(unlist(strsplit(TOBIAS_merge$V1,split="\\."))),2)]
pfm_honeybee <- readJASPARMatrix("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_promoter_TOBIAS/merged_honeybee_pfam.jaspar", matrixClass=c("PFM", "PWM", "PWMProb"))
motif2TF<- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_promoter_TOBIAS/motif2TF.txt")
# get the TF info 
TOBIAS_merge$TF<- motif2TF[match(TOBIAS_merge$V5,motif2TF$V1),]$V2
fly2honeybee<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/05_fly2honeybee.csv")
TOBIAS_merge$TF_name<- TOBIAS_merge$TF

for (i in which(TOBIAS_merge$TF%in%fly2honeybee$fly_gene)){
	TOBIAS_merge$TF_name[i]<- fly2honeybee[match(TOBIAS_merge$TF[i],fly2honeybee$fly_gene),]$honeybee_gene_name
}

Avg <- AverageExpression(ORN,features=rownames(ORN),assays = "raw_RNA")
# add exp info 
Avg_data<- Avg$raw_RNA;
colnames(Avg_data)<- levels(ORN)

dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_data<- dotplot_data[dotplot_data$id!="33",]
TOBIAS_merge<- TOBIAS_merge[TOBIAS_merge$OR %in% dotplot_data$features.plot,]


TOBIAS_merge$exp<- NA
TOBIAS_merge$FC<- NA
for(nrow in 1:nrow(TOBIAS_merge)){
	OR<- TOBIAS_merge[nrow,]$OR
	cluster<- unique(dotplot_data[dotplot_data$features.plot==OR,]$id)
	if(length(cluster)==1){
	if(length(strsplit(TOBIAS_merge$TF_name[nrow],"::")[[1]])!=1){
		TF<- strsplit(TOBIAS_merge$TF_name[nrow],"::")[[1]];
		if(length(which(TF%in%rownames(Avg_data)))){
		  exp<- mean(Avg_data[TF%in%rownames(Avg_data),cluster])
		  exp_in_cluster<- exp
		  exp_other <- (sum(Avg_data[TF%in%rownames(Avg_data),])-exp_in_cluster)/68
      TOBIAS_merge$FC[nrow]<- exp_in_cluster/exp_other;
      TOBIAS_merge$exp[nrow]<- exp;
		}
		}else{
			TF<- TOBIAS_merge[nrow,]$TF_name;
			if(TF%in% rownames(Avg_data)){
			exp<- Avg_data[TF,cluster]
			exp_in_cluster<- Avg_data[TF,cluster]
			exp_other <- (sum(Avg_data[TF,])-exp_in_cluster)/58
	    TOBIAS_merge$exp[nrow]<- exp;
	    TOBIAS_merge$FC[nrow]<- exp_in_cluster/exp_other;
			}
	}
	# 
	if(length(cluster)>1){
		if(length(strsplit(TOBIAS_merge$TF_name[nrow],"::")[[1]])!=1){
		TF<- strsplit(TOBIAS_merge$TF_name[nrow],"::")[[1]];
		if(length(which(TF%in%rownames(Avg_data)))){
		  exp<- mean(Avg_data[TF[which(TF%in%rownames(Avg_data))],cluster])
		  exp_in_cluster<- exp
		  exp_other <- (sum(Avg_data[TF[which(TF%in%rownames(Avg_data))],])-exp_in_cluster*length(cluster))/58
      TOBIAS_merge$FC[nrow]<- exp_in_cluster/exp_other;
      TOBIAS_merge$exp[nrow]<- exp;
		}
		}else{
			TF<- TOBIAS_merge[nrow,]$TF_name;
			if(TF%in% rownames(Avg_data)){
			exp<- mean(Avg_data[TF,cluster])
			exp_in_cluster<- exp
			exp_other <- (sum(Avg_data[TF,])-exp_in_cluster*length(cluster))/58
	    TOBIAS_merge$exp[nrow]<- exp;
	    TOBIAS_merge$FC[nrow]<- exp_in_cluster/exp_other;
			}
	}
		}
}

ORN_avg<-AverageExpression(
       ORN,
       assays = "raw_RNA",
       features = TOBIAS_merge$TF_name,
       return.seurat = FALSE,
       group.by = "subcluster",
       #add.ident = NULL,
       slot = "data")
ORN_avg<-ORN_avg$raw_RNA
colnames(ORN_avg)<- levels(ORN)
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
TOBIAS_merge$tau<- gene_tau[TOBIAS_merge$TF_name]


write.csv(TOBIAS_merge,"./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/TOBIAS_merge_OR_bubbleplot.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
da_data<- na.omit(TOBIAS_merge)
da_data$OR<- factor(da_data$OR,levels=dotplot_feature)
summary(da_data$V6)
top5 <- da_data %>% group_by(OR) %>% top_n(n = 5, wt = V6);
da_data<- da_data[da_data$V6>8,]
da_data<- da_data[order(da_data$OR),]
da_data$TF_name<- factor(da_data$TF_name,levels=unique(da_data$TF_name))

pdf("./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/TOBIAS_merge_OR_bubbleplot_bubbleplot-FC.pdf",width=10,height=16)
ggplot(da_data,aes(x=OR,y=TF_name))+
  geom_point(aes(size=tau,color=log2(FC+1)),alpha=0.8)+
  scale_size(range=c(1,12))+
  scale_color_gradientn(colours=c("#080C26","#A61C68","#D92958","#F26E50","#F2E0D5"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
dev.off()

# part3: 
# TF add the family info 
TF_list_info<- read.table("./05_ORN_cluster2/07_DEG_and_DEP/TF_Information.txt",sep="\t",header=T)
TOBIAS_merge$Family_Name<- TF_list_info[match(TOBIAS_merge$V5,TF_list_info$Motif_ID),]$Family_Name




differential.activity$Motif_ID<- 
differential.activity$Family_Name<- TF_list_info[match(TOBIAS_merge$V5,TF_list_info$Motif_ID),]$Family_Name





