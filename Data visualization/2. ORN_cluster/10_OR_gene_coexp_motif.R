library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)
set.seed(1234);
library(pheatmap)

coexp_motif<- read.csv("./07_coexp_motif/Or_pair_wehave_heatmapdata.csv")
coexp_motif<- coexp_motif[,-4]
coexp_motif<- coexp_motif[!duplicated(coexp_motif),]
coexp_motif2<- coexp_motif[,c(2,1,3)]
colnames(coexp_motif2)<- c("file1","file2","count")
coexp_motif3<- rbind(coexp_motif,coexp_motif2)
gene<- unique(c(coexp_motif$file1,coexp_motif$file2))

data<- matrix(ncol=length(gene),nrow=length(gene))
rownames(data)<- gene
colnames(data)<- gene
for (i in 1:length(gene)){
	for (j in 1:length(gene)){
		gene1<- gene[i]
		gene2<- gene[j]
		if(gene1==gene2){data[i,j]==NA}else{
		data[i,j]<- coexp_motif3[which(coexp_motif3$file1==gene1&coexp_motif3$file2==gene2),3] }
	}
}
library(pheatmap)
OR_count<- read.csv("./07_coexp_motif/Or_pair_wehave_singlecount.csv")
for(g in gene){
	data[g,g]<- OR_count[which(OR_count$file1==g),2]
}

count=t(scale(t(data),scale = T,center = F))

pdf("./07_coexp_motif/coexp_motif.pdf",width=11,height=10)
pheatmap(log2(data+1),
      cluster_cols = T,
      cluster_rows = T,
      color = colorRampPalette(c("white", "firebrick3"))(100),
      #annotation_col = barcode_label_pheatmap,
      #annotation_colors = ann_colors,
      #annotation_row = barcode_label_pheatmap,
      annotation_legend = TRUE,
      show_rownames=T,
      show_colnames=T
 )
dev.off()


# cluster
TF<- c("Twi","LOC726594","LOC552079","LOC551928","LOC552079","LOC726999")
Avg <- AverageExpression(ORN,features=TF,assays = "raw_RNA")
library(pheatmap)
count=t(scale(t(Avg$raw_RNA),scale = T,center = F))
pdf("./07_coexp_motif/P3_3_TF_exp.pdf",width=15,height=3)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
dev.off();
