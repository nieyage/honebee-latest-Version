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
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C"  , "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
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
#OR2 is placed in the last column;
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
# corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature,]
# color by the group info 
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
# Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
# rownames(Group_info)<-Group_info$gene_name
# groupInfo <- split(row.names(Group_info), Group_info$seqnames)
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)
#data.tree <- groupOTU(tree, groupInfo)
data.tree <- tree
pdf("./00_Figure/Fig3/Fig3A-a-OR_sequence_protein_similarity-tree_add_groupinfo_heatmap.pdf",width=25,height=16)
ggtree(data.tree,ladderize = FALSE, branch.length = "none") + geom_tiplab(size=3) + 
theme(legend.position = "right")+ scale_color_manual (values =myUmapcolors[5:30]) +  geom_treescale()
dev.off()
m<-ggtree(data.tree,ladderize = FALSE, branch.length = "none")+ geom_tiplab(size=3) + theme(legend.position = "right")+ scale_color_manual (values =myUmapcolors[5:30]) 
gene_order<-na.omit(m$data[order(m$data$y),]$label)
gene_order<-as.character(gene_order);

# OR color bar in heatmap 

library(pheatmap)
DefaultAssay(ORN)<-"SCT"
matrix<- ORN@assays$SCT[dotplot_feature,]
cor_data<-cor(as.data.frame(t(matrix)))
dist_data <- cor_data[rev(gene_order),rev(gene_order)]
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
label_pheatmap<- data.frame(Group=Group_info$seqnames)
rownames(label_pheatmap) <- Group_info$gene_name

ann_colors<-list(Group =    c("Group1"=myUmapcolors[1],    "Group2"=myUmapcolors[2],    "Group3"=myUmapcolors[3],    "Group4"=myUmapcolors[4],    "Group5"=myUmapcolors[5],    "Group6"=myUmapcolors[6],    "Group7"=myUmapcolors[7],    "Group8"=myUmapcolors[8],    "Group9"=myUmapcolors[9],    "Group10"=myUmapcolors[10],
    "Group11"=myUmapcolors[11],    "Group12"=myUmapcolors[12],    "Group13"=myUmapcolors[13],    "Group14"=myUmapcolors[14],    "Group15"=myUmapcolors[15],    "Group16"=myUmapcolors[16],    "GroupUN3"=myUmapcolors[17],
    "GroupUN226"=myUmapcolors[18],    "GroupUN243"=myUmapcolors[19],    "GroupUN248"=myUmapcolors[20]))

pdf("./00_Figure/Fig3/Fig3A-b-OR_in_dotplot_sequence_correlation.pdf",width=17,height=16)
pheatmap(dist_data,
        border = F,
         cluster_cols = F,
         cluster_rows = F,
         cellwidth = 8, cellheight = 8,
         color = colorRampPalette(c("white", "#3082BD","#1C214F"))(100),
         annotation_legend = TRUE,
         annotation_colors = ann_colors,
         annotation_row = label_pheatmap,
         show_rownames=T,
         show_colnames=T
    )
dev.off()

# Fig3B:
cluster_number<-as.data.frame(table(table(dotplot_data$id)))
cluster_number$Var1<-as.character(cluster_number$Var1)
cluster_number[4,1]<-">3"
cluster_number[4,2]<-2
cluster_number<-cluster_number[1:4,]
cluster_number$percent<- cluster_number$Freq/sum(cluster_number$Freq)*100;
cluster_number$Var1<- factor(cluster_number$Var1,levels=c("1","2","3",">3"))


pdf("./00_Figure/Fig3/Fig3B-multi_OR_percent.pdf",width=3,height=4)
p<-ggplot(data = cluster_number, aes_string(x = "Var1", y = "percent")) +  
        xlab("# chemosensory receptor expressed") +
        ylab("% percent of cluster") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity",width=0.6,color = 'black', fill='grey') +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
round(cluster_number$percent, 2)
p+geom_text(aes(label = Freq), size = 3, hjust = 0.5, vjust = 3) 
dev.off();
# Fig3C:
OR_number<-as.data.frame(table(table(dotplot_data$features.plot)));
#OR_number<-OR_number[-1,]
OR_number$percent<- OR_number$Freq/sum(OR_number$Freq)*100;
pdf("./00_Figure/Fig3/Fig3C-OR_exp_in_multicluster_percent.pdf",width=3,height=4)
p<-ggplot(data = OR_number, aes_string(x = "Var1", y = "percent")) +  
        xlab("# of cluster OR expressed") +
        ylab("% percent of chemosensory receptor") + 
        #scale_fill_manual(values = "#FFED6F") + 
        geom_bar( stat = "identity",width=0.6,color = 'black', fill='grey') +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5));
p
#add number in plot 
p+geom_text(aes(label = Freq), size = 3, hjust = 0.5, vjust = 3) 
dev.off();


colors_for_exp_pattern<- c("#476D87","#E95C59")
#Fig 3D: Sequence similarity

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
type<-c(rep("co-exp",length(same_cluster)),rep("non co-exp",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data<-data.frame(type,var)
data$type<-factor(data$type,levels=c("non co-exp","co-exp"))
library(ggpubr)
library(cowplot)

pdf("./00_Figure/Fig3/Fig3D-OR_sequence_similarity_distribution.pdf",width=3,height=3)
p1<-ggplot(data, aes(x=var, fill=type)) + xlab("% OR sequence similarity")+
          geom_density(alpha=.25) + theme_classic()+theme(legend.position="top")+
scale_color_manual(values =colors_for_exp_pattern)
p2<-ggboxplot(data, x="type", y="var", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("% OR sequence similarity")+
scale_color_manual(values =colors_for_exp_pattern)
p2
dev.off()


# the RMSD (structure similarity)
RMSD <- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/rename_pdb/RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/rename_pdb/pdb_list.txt")
ID <- gsub(".pdb","",pdb$V1)
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
#data[data==0]<- 0.5
#data<- log2(data)
same_cluster <- c()
not_same_cluster <- c()
remaining_gene<-rownames(data)
for(gene1 in remaining_gene){
  remaining_gene<-remaining_gene[-which(remaining_gene%in%gene1)]
  for (gene2 in remaining_gene){
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id
    dist<-data[gene1,gene2]
    if(length(which(duplicated(cluster_id)))){
        same_cluster<-c(same_cluster,dist)
    }else{not_same_cluster<-c(not_same_cluster,dist)}
  }
}
type<-c(rep("co-exp",length(same_cluster)),rep("non co-exp",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data2<-data.frame(type,var)
data2$type<-factor(data2$type,levels=c("non co-exp","co-exp"))
pdf("./00_Figure/Fig3/Fig3E-OR_structure_similarity_distribution.pdf",width=3,height=3)
ggboxplot(data2, x="type", y="var", color = "type",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("OR RMSD")+
scale_color_manual(values = colors_for_exp_pattern)
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
        OR_pair$same_cluster[i]="co-exp";
    }else{OR_pair$same_cluster[i]="non co-exp";}
    
};

# The probability of OR pair appearing on the same chromosome
#> table(OR_pair$same_seqname,OR_pair$same_cluster)
      co-exp non co-exp
  No       0       1749
  Yes     46        416


probability<-data.frame(OR_pair_type=c("non co-exp","co-exp"),
    probability=c(416/2165,46/46),count=c(46,416));
probability$OR_pair_type<-factor(probability$OR_pair_type,levels=c("non co-exp","co-exp"))
pdf("./00_Figure/Fig3/Fig3F-The_probability_OR_pair_appearing_on_same_chromosome.pdf",width=6,height=3)
p<-ggplot(data = probability, aes_string(x = "OR_pair_type", y = "probability",color="OR_pair_type")) +  
        xlab("OR pair type") +
        ylab("The probability on the same chromosome") + 
        geom_bar(stat = "identity", width = 0.6,fill="white") +
        theme_classic()+guides(color='none')+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5)) + 
        scale_color_manual(values = colors_for_exp_pattern);
#add number in plot 
a<- p+geom_text(aes(label = count), size = 3, hjust = 0.5, vjust = 3)
# the distance 
same_chrom<-OR_pair[OR_pair$same_seqname=="Yes",]
for (i in 1:nrow(same_chrom)){
    gene1<-same_chrom[i,]$gene1;
    gene2<-same_chrom[i,]$gene2;
    OR_pair_transcript<-OR_gene_transcript[OR_gene_transcript$gene_name%in%c(gene1,gene2),]
    OR_pair_transcript<-sort(OR_pair_transcript);
    same_chrom$dist[i]<-width(gaps(OR_pair_transcript))[2];
}

same_chrom$same_cluster<-factor(same_chrom$same_cluster,levels=c("non co-exp","co-exp"))
b<- ggboxplot(same_chrom, x="same_cluster", y="dist",color = "same_cluster",width=0.6, notch = F)+
stat_compare_means()+theme(legend.position="none")+ylab("OR transcript distance")+
scale_color_manual(values = colors_for_exp_pattern)
a|b 
dev.off()

# Fig3F two promoter  (first)
obj<- subset(ORN,idents=c("p3:3"))
DefaultAssay(obj)<-"raw_RNA"
gene<- c("LOC107965761","LOC102655285")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3G-a-coexp-heatmap-SCT.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,
     border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_features<- gene
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!="p3:3"){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  genes.use = obj_features
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3G-b-multi_OR_without_nopower_trackplot.pdf",width=12,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream =500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
set<- colors_for_exp_pattern[c(2,1)]
p1
dev.off()

# Fig3F two promoter (second)
cluster="p3:4_1" #p2:21
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot

DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-MP2-coexp-heatmap-SCT.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
# Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
##Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
## link peaks to genes
#obj <- LinkPeaks(
#  object = obj,
#  peak.assay = "peaks_ORN_subcluster",
#  expression.assay = "raw_RNA",
#  genes.use = obj_features
#)
#######Visulize track and RNA exp######
idents.plot <- Idents(obj)
## plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-MP2-multi_OR_without_nopower_trackplot2.pdf",width=12,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream =500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()
# Fig3F two promoter (third)
cluster="p2:21" #p2:21
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot

DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-MP3-coexp-heatmap-SCT.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
# Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
##Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
## link peaks to genes
#obj <- LinkPeaks(
#  object = obj,
#  peak.assay = "peaks_ORN_subcluster",
#  expression.assay = "raw_RNA",
#  genes.use = obj_features
#)
#######Visulize track and RNA exp######
idents.plot <- Idents(obj)
## plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-MP3-multi_OR_without_nopower_trackplot2.pdf",width=12,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream =2000,
  annotation = TRUE,
  extend.downstream = 2000,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()
#Fig3H single promoter (first)
cluster="p1:14"
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-a-coexp-heatmap-SCT-singlepromoter.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  genes.use = obj_features
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-b-multi_OR_without_nopower_trackplot-singlepromoter.pdf",width=10,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()


#Fig3H single promoter (second)
cluster="p4:5"
obj<- subset(ORN,idents=cluster)
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
DefaultAssay(obj)<-"raw_RNA"
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
pdf("./00_Figure/Fig3/Fig3H-SP2-coexp-heatmap-SCT-singlepromoter.pdf",height=3,width=12)
pheatmap(ORN_matrix,
     cluster_cols = T,
     cluster_rows = T,border=F,
     color = colorRampPalette(c("white", "#CC0000"))(100),
     annotation_legend = TRUE,
     show_rownames=T,
     show_colnames=F
  )
dev.off()

# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
obj_barcode<-colnames(obj)
all_barcode<-colnames(ORN)
random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
obj$subcluster<-as.character(obj$subcluster)
for (i in 1:length(obj$subcluster)){
  if(obj$subcluster[i]!=cluster){obj$subcluster[i]="other"}
    }
Idents(obj)<-obj$subcluster
# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
#Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  genes.use = obj_features
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
# plot region 
start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
ranges.show <- paste(seq,start,end,sep="-")
pdf("./00_Figure/Fig3/Fig3H-SP2-multi_OR_without_nopower_trackplot-singlepromoter.pdf",width=10,height=6)
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 500,
  annotation = TRUE,
  extend.downstream = 500,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
dev.off()

# # version 2023.10.7
# Fig3H:
pdf('./00_Figure/Fig3/Fig3H-multiple-OR-FeaturePlot.pdf', width=16, height=4)
# Visualize co-expression of two features simultaneously
FeaturePlot(ORN, features = c("LOC102653695", "LOC102653615"),cols=c("lightgrey", "#E31A1C", "#4DAE49"), max.cutoff =3, blend = TRUE,order=TRUE,)+ggtitle("p3_3:Or85b_LOC102655285")
dev.off()

# Fig3I:
DefaultAssay(ORN) <- "SCT"
### Plotting scatter plot: 1. single-cell level, 2. cluster level
ORN$LOC102653695_UMIs <- ORN@assays$SCT@counts['LOC102653695',]
ORN$LOC102653615_UMIs <- ORN@assays$SCT@counts['LOC102653615',]
#.....................................................................................
#  LOC107965761 vs. LOC102655285 UMI
# ....................................................................................
#  single-cell level
library(cowplot)
pmain <- ORN@meta.data %>%
  ggplot( aes(LOC102653695_UMIs, LOC102653615_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = ORN@meta.data, aes(x = LOC102653695_UMIs), fill="#B31416") +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ORN@meta.data, aes(x = LOC102653615_UMIs), fill="#4DAE49") + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.3, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.3, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/Fig3/Fig3I-3-LOC102653695vsLOC102653615_UMI.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 1000, height =1000,p3)











