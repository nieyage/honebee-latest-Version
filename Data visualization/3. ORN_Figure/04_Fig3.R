library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)

set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_orderbytree.rds")
DefaultAssay(ORN)<-"raw_RNA"

#Fig3A OR correlation  heatmap and sequence tree 
# sequence tree 
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
dotplot_data<-read.csv("./05_ORN_cluster/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
#OR2 is placed in the last column;
dotplot_feature_OR2<-c(Orco,dotplot_feature)
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
# corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature_OR2,]
# color by the group info 
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
rownames(Group_info)<-Group_info$gene_name
groupInfo <- split(row.names(Group_info), Group_info$seqnames)
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
tree<-as.phylo(clust)
data.tree <- groupOTU(tree, groupInfo)
pdf("./00_Figure/Fig3A-OR_sequence_protein_similarity-tree_add_groupinfo_heatmap.pdf",width=25,height=16)
ggtree(data.tree,ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab(size=3) + theme(legend.position = "right")
dev.off()
m<-ggtree(data.tree,ladderize = FALSE, branch.length = "none",aes(color=group))
gene_order<-na.omit(m$data[order(m$data$y),]$label)
gene_order<-as.character(gene_order)
library(pheatmap) 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_data <- dist_percent[rev(gene_order),rev(gene_order)]
pdf("./00_Figure/Fig3A-OR_in_dotplot_sequence_correlation.pdf",width=17,height=16)
pheatmap(dist_data,
         cluster_cols = F,
         cluster_rows = F,
         #color = colorRampPalette(c("#F1ECEC", "#3082BD","#1C214F"))(100),
         #annotation_col = barcode_label_pheatmap,
         #annotation_colors = ann_colors,
         #annotation_row = barcode_label_pheatmap,
         annotation_legend = TRUE,
         show_rownames=T,
         show_colnames=T
    )
dev.off()







##add a OR symbol label 
## add max_exp OR label for each cell
#ORN_count<-ORN@assays$RNA
#ORN_count<-ORN_count[which(rownames(ORN_count)%in%dotplot_feature),]
#ORN_matrix<-as.matrix(ORN_count)
##ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
##ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
#barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("Not_sure",length(colnames(ORN_matrix))))
#  for (i in 1:length(colnames(ORN_matrix))){
#    if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
#    barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
#}
#  }
#ORN$OR_label<-barcode_label$label
#object<-ORN
#embeddings <- Embeddings(object = object, reduction = "pca")[,1:50]
#Idents(object)<-ORN$OR_label
#data.dims <- lapply(X = levels(x = object), FUN = function(x) {
#    cells <- WhichCells(object = object, idents = x)
#    if (length(x = cells) == 1) {
#        cells <- c(cells, cells)
#    }
#    temp <- colMeans(x = embeddings[cells, ])
#})
#data.dims <- do.call(what = "cbind", args = data.dims)
#colnames(x = data.dims) <- levels(x = object)
#last_data<-data.dims[,-which(colnames(data.dims)=="Not_sure")]
#cor.data.dims <- as.data.frame(abs(cor(last_data)))
#cor.data.dims <- cor.data.dims[rev(gene_order),rev(gene_order)]
#pdf("./00_Figure/Fig3A-OR_in_dotplot_cosine_correlation.pdf",width=17,height=16)
#pheatmap(cor.data.dims,
#         cluster_cols = F,
#         cluster_rows = F,
#         #color = colorRampPalette(c("#F1ECEC", "#3082BD","#1C214F"))(100),
#         #annotation_col = barcode_label_pheatmap,
#         #annotation_colors = ann_colors,
#         #annotation_row = barcode_label_pheatmap,
#         annotation_legend = TRUE,
#         show_rownames=T,
#         show_colnames=T
#    )
#dev.off()

# Fig3B:
cluster_number<-as.data.frame(table(table(dotplot_data$id)))
cluster_number$percent<- cluster_number$Freq/sum(cluster_number$Freq)*100;
pdf("./00_Figure/Fig3B-multi_OR_percent.pdf",width=4,height=4)
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "percent")) +  
        xlab("# chemosensory receptor expressed") +
        ylab("% percent of cluster") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity", width = 0.6, fill="#69b3a2") +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
round(cluster_number$percent, 2)
p+geom_text(aes(label = round(percent, 2)), size = 3, hjust = 0.5, vjust = 3) 
dev.off();
# Fig3C:
OR_number<-as.data.frame(table(table(dotplot_data$features.plot)));
#OR_number<-OR_number[-1,]
OR_number$percent<- OR_number$Freq/sum(OR_number$Freq)*100;
pdf("./00_Figure/Fig3C-OR_exp_in_multicluster_percent.pdf",width=4,height=4)
p<-ggplot(data = OR_number, aes_string(x = "Var1", y = "percent")) +  
        xlab("# of cluster OR expressed") +
        ylab("% percent of chemosensory receptor") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity", width = 0.6, fill="#BC80BD") +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
p+geom_text(aes(label = round(percent, 2)), size = 3, hjust = 0.5, vjust = 3) 
dev.off();

#Fig 3D,3E
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature,]
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
tree<-as.phylo(clust)
# the similarity of OR pair in the same cluster and random OR pair 
# All OR pair similarity matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
same_cluster <- c()
not_same_cluster <- c()
remaining_gene<-dotplot_feature
for(gene1 in dotplot_feature){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id
    dist<-dist_percent[gene1,gene2]
    if(length(which(duplicated(cluster_id)))){
        same_cluster<-c(same_cluster,dist)
    }else{not_same_cluster<-c(not_same_cluster,dist)}
  }
}
type<-c(rep("same_cluster",length(same_cluster)),rep("not_same_cluster",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("same_cluster","not_same_cluster"))
library(ggpubr)
library(cowplot)
pdf("./00_Figure/Fig3D-OR_sequence_similarity_distribution.pdf",width=8,height=5)
p1<-ggplot(data, aes(x=var, fill=type)) + xlab("% OR sequence similarity")+
          geom_density(alpha=.25) + theme_classic()+theme(legend.position="top")+
scale_color_manual(values =c("#F55050","#86A3B8"))
p2<-ggboxplot(data, x="type", y="var", color = "type",width=0.6, notch = T)+
stat_compare_means()+theme(legend.position="none")+ylab("% OR sequence similarity")+
scale_color_manual(values =c("#F55050","#86A3B8"))
plot_grid(p1,p2,labels = c(" "," "),rel_widths = c(1.5, 1))
dev.off()
#the distance of OR pair in the same cluster and random OR pair 
DefaultAssay(ORN)<-"ATAC"
gene_transcript<-Annotation(ORN)[which(Annotation(ORN)$type=="transcript"),];
OR_gene_transcript<-gene_transcript[gene_transcript$gene_name%in%dotplot_feature,]
OR_gene_transcript<-OR_gene_transcript[!duplicated(OR_gene_transcript$gene_name),]
AmelchrNameLength <- read.table("/md01/nieyg/ref/10X/honeybee/de_novo_antenna/Amel_antenan/star/chrNameLength.txt", sep="\t", stringsAsFactors=F) 
AmelchrNameLength<-AmelchrNameLength[which(AmelchrNameLength$V1%in%as.character(unique(OR_gene_transcript@seqnames@values))),]
data<-as.data.frame(sort(OR_gene_transcript));
data<-data[,c(1:3,12,4,5)]
write.table(data,"OR_gene_transcript.bed",sep="\t",row.names=F,col.names=F)

# #Shell
# sed -i 's/"//g' OR_gene_transcript.bed;
# bedtools sort -i OR_gene_transcript.bed > OR_gene_transcript_sorted.bed
# bedtools closest -s -d -io -N -a OR_gene_transcript_sorted.bed -b OR_gene_transcript_sorted.bed > output.bed
# awk '{print $NF,"\t",$1,"\t",$4,"\t",$10}' output.bed > closestOlfrGenes.txt
# sort -n closestOlfrGenes.txt | awk '$1 > 0 {print $0}' > sortedClosestOlfrGenes.txt

# back R 
# make OR pair 
OR_pair<-data.frame()
remaining_gene<-dotplot_feature
for(gene1 in dotplot_feature){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    data_subset<-data.frame(gene1=gene1,gene2=gene2);
    OR_pair<-rbind(OR_pair,data_subset)
  }
}
OR_pair$gene1_seqname <- data[match(OR_pair$gene1,data$gene_name),1]
OR_pair$gene2_seqname <- data[match(OR_pair$gene2,data$gene_name),1]

#Same seqname?
for (i in 1:nrow(OR_pair)){
    if(OR_pair[i,]$gene1_seqname==OR_pair[i,]$gene2_seqname){
        OR_pair$same_seqname[i]="Yes";}
        else{OR_pair$same_seqname[i]="No"}
};
#Same cluster?
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$gene1;
    gene2<-OR_pair[i,]$gene2;
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id;
    if(length(which(duplicated(cluster_id)))){
        OR_pair$same_cluster[i]="Yes";
    }else{OR_pair$same_cluster[i]="No";}
    
};

# The probability of OR pair appearing on the same chromosome
#> table(OR_pair$same_seqname,OR_pair$same_cluster)
#     
#        No  Yes
#  No  1709    0
#  Yes 1099   42
#        No  Yes
#  No  4046    0
#  Yes 1607  125

probability<-data.frame(OR_pair_type=c("same cluster","different cluster"),
    probability=c(125/125,1607/5653));

pdf("./00_Figure/Fig3E-The_probability_OR_pair_appearing_on_same_chromosome.pdf",width=4,height=4)
p<-ggplot(data = probability, aes_string(x = "OR_pair_type", y = "probability")) +  
        xlab("OR pair type") +
        ylab("The probability on the same chromosome") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity", width = 0.6, fill="#80B1D3") +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
p+geom_text(aes(label = round(probability, 2)), size = 3, hjust = 0.5, vjust = 3) +
scale_fill_manual(values =c("#F55050","#86A3B8"))
dev.off();
# the distance 
same_chrom<-OR_pair[OR_pair$same_seqname=="Yes",]
for (i in 1:nrow(same_chrom)){
    gene1<-same_chrom[i,]$gene1;
    gene2<-same_chrom[i,]$gene2;
    OR_pair_transcript<-OR_gene_transcript[OR_gene_transcript$gene_name%in%c(gene1,gene2),]
    OR_pair_transcript<-sort(OR_pair_transcript);
    same_chrom$dist[i]<-width(gaps(OR_pair_transcript))[2];
}

same_chrom$same_cluster<-factor(same_chrom$same_cluster,levels=c("Yes","No"))
library(ggpubr)
library(cowplot)
pdf("./00_Figure/Fig3E-OR_transcript_distance_distribution.pdf",width=5,height=5)
#p1<-ggplot(same_chrom, aes(x=dist, fill=same_cluster)) + xlab("OR transcript distance")+
#          geom_density(alpha=.25,outlier.shape = NA) + theme_classic()+theme(legend.position="top")
#p2<-ggboxplot(same_chrom, x="same_cluster", y="dist", color = "same_cluster",width=0.6, notch = T)+
#stat_compare_means()+theme(legend.position="none")+ylab("OR transcript distance")
#plot_grid(p1,p2,labels = c(" "," "),rel_widths = c(1.5, 1))
ggboxplot(same_chrom, x="same_cluster", y="dist",color = "same_cluster",width=0.6, notch = T)+
stat_compare_means()+theme(legend.position="none")+ylab("OR transcript distance")+
scale_color_manual(values =c("#F55050","#86A3B8"))
dev.off()

# Fig3F



