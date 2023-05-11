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
DefaultAssay(ORN)<-"raw_RNA"
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

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
dotplot_feature_OR2<-c(Orco,dotplot_feature)
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
## #add a OR symbol label 
## # add max_exp OR label for each cell
# ORN_count<-ORN@assays$RNA
# ORN_count<-ORN_count[which(rownames(ORN_count)%in%dotplot_feature),]
# ORN_matrix<-as.matrix(ORN_count)
# #ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
# #ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
# barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("Not_sure",length(colnames(ORN_matrix))))
#   for (i in 1:length(colnames(ORN_matrix))){
#     if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
#     barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))
# }
#   }
# ORN$OR_label<-barcode_label$label

barcode_label<-c()
for (i in 1:ncol(ORN_matrix)){
  if (max(ORN_matrix[,i])>0){barcode_label[[i]]<-names(which(ORN_matrix[,i]== max(ORN_matrix[,i])))}
  }
names(barcode_label)<-colnames(ORN_matrix)
label_all<-unlist(barcode_label)
embeddings <- Embeddings(object = ORN, reduction = "pca")[,1:50]
data.dims <- lapply(X = dotplot_feature, FUN = function(x) {
    barcode_gene<- names(which(label_all==x))
    cells <- gsub("-1.*","-1",barcode_gene)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- dotplot_feature
#last_data<-data.dims[,-which(colnames(data.dims)=="Not_sure")]
cor.data.dims <- as.data.frame(abs(cor(data.dims)))

# OR sequence corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature,]
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
tree<-as.phylo(clust)
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
groupInfo <- split(Group_info$gene_name, Group_info$seqnames)
data.tree <- groupOTU(tree, groupInfo)

pdf("./00_Figure/Fig3A-OR_sequence_protein_similarity-tree_add_groupinfo_heatmap.pdf",width=25,height=16)
ggtree(data.tree,ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab(size=3) + theme(legend.position = "right")
dev.off()
m<-ggtree(data.tree,ladderize = FALSE, branch.length = "none",aes(color=group))
gene_order<-na.omit(m$data[order(m$data$y),]$label)
gene_order<-as.character(gene_order)

# transcript_dist heatmap 
library(pheatmap)
cor.data.dims <- cor.data.dims[rev(gene_order),rev(gene_order)]
pdf("./00_Figure/Fig3A-OR_in_dotplot_cosine_correlation.pdf",width=17,height=16)
pheatmap(cor.data.dims,
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
# within_group/between_group value heatmap 
embeddings <- Embeddings(object = ORN, reduction = "pca")[,1:50]
trans_dist <- 1-cosine(t(embeddings))
# calculate the cosine simility within group and between groups;
log2FCdata<-matrix(nrow=length(dotplot_feature),ncol=length(dotplot_feature))
colnames(log2FCdata)<- dotplot_feature
rownames(log2FCdata)<- dotplot_feature
for (gene1 in dotplot_feature){
    for(gene2 in dotplot_feature){
        gene1_barcode <- names(which(label_all==gene1))
        gene1_cells <- gsub("-1.*","-1",gene1_barcode)
        gene2_barcode <- names(which(label_all==gene2))
        gene2_cells <- gsub("-1.*","-1",gene2_barcode)
        within_group<-c(trans_dist[gene1_cells,gene1_cells],trans_dist[gene2_cells,gene2_cells])
        between_group<-c(trans_dist[gene1_cells,gene2_cells],trans_dist[gene2_cells,gene1_cells])
        # calculate the FC 
        if(!is.null(mean(between_group))){
        log2FC<-log2(mean(within_group))-log2(mean(between_group));
        log2FCdata[gene1,gene2]<-log2FC
        }

    }
}
log2FCdata <- log2FCdata[rev(gene_order),rev(gene_order)]
pdf("./00_Figure/Fig3A-OR_in_dotplot_cosine-log2FCdata.pdf",width=17,height=16)
pheatmap(log2FCdata,
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
data[data==0]<- 0.5
data<- log2(data)
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

type<-c(rep("same_cluster",length(same_cluster)),rep("not_same_cluster",length(not_same_cluster)))
var<-c(same_cluster,not_same_cluster)
data2<-data.frame(type,var)
data2$type<-factor(data2$type,levels=c("same_cluster","not_same_cluster"))
library(ggpubr)

library(cowplot)
pdf("./00_Figure/Fig3F-OR_structure_similarity_distribution.pdf",width=3,height=3)
ggboxplot(data2, x="type", y="var", color = "type",width=0.6, notch = T)+
stat_compare_means()+theme(legend.position="none")+ylab("% OR structure similarity")+
scale_color_manual(values =c("#F55050","#86A3B8"))
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
#        No  Yes
#  No  4678    1
#  Yes 1880  344

probability<-data.frame(OR_pair_type=c("same cluster","different cluster"),
    probability=c(344/(344+1),1880/(1880+4678)));

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


# the clear co-exp and the Combinations
#filter dotplot ;
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
Freq<- as.data.frame(table(dotplot_data$id));
clear_multiple_classes_order_id <- as.character(Freq[Freq$Freq>1,1])

# Fig3F
# coexp in multiple OR cluster
#Fig 3FGH
# show a typical OR pair (co-expression);
# select a beautiful track to show :
library(ggplot2)
library(gggenes)
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
library(RColorBrewer)
log2FCdata<-data.frame();
DefaultAssay(ORN)<- "integratedRNA_onecluster"
pdf("./05_ORN_cluster2/03_coexp_cluster/distinguish_multi_OR_with_powerful.pdf",width=14,height=6)
for (cluster in clear_multiple_classes_order_id){
print(cluster)
obj<-subset(ORN,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
# add max_exp OR label for each cell
  ORN_count<-obj@assays$SCT
  ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
  ORN_matrix<-as.matrix(ORN_count)
  ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
  ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
  barcode_label<-data.frame(barcode=colnames(ORN_matrix),label=rep("NA",length(colnames(ORN_matrix))))
    for (i in 1:length(colnames(ORN_matrix))){
      if(length(names(which(ORN_matrix[,i]==max(ORN_matrix[,i]))))==1){
      barcode_label[i,2]=names(which(ORN_matrix[,i]==max(ORN_matrix[,i])))}
    }
    barcode_label<-barcode_label[barcode_label$label!="NA",]
    barcode_label<-barcode_label[order(barcode_label$label),];
# p1 cell cosine simility heatmap 
    embeddings <- Embeddings(object = ORN, reduction = "pca")[,1:50]
    embeddings <- embeddings[barcode_label$barcode,]
    trans_dist <- 1-cosine(t(embeddings))
    barcode_label_pheatmap<-data.frame(OR=barcode_label$label)
    rownames(barcode_label_pheatmap)<-barcode_label$barcode
    col<-brewer.pal(12,"Set3")[1:length(unique(barcode_label_pheatmap$OR))]
    names(col)<-unique(barcode_label_pheatmap$OR)
    ann_colors= list(OR = col)
    p1<-pheatmap(trans_dist,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
             color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=F,
             show_colnames=F)
# calculate the cosine simility within group and between groups;
    rownames(barcode_label)<-barcode_label$barcode
    within_group<-c()
    between_group<-c()
    for (i in 1:nrow(trans_dist)){
    for (j in 1:ncol(trans_dist)){
        if(i!=j){
            if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
                within_group<-c(within_group,trans_dist[i,j]);
            }
            else{between_group<-c(between_group,trans_dist[i,j])}
        }
    }
    }
# calculate the FC 
    if(!is.null(median(between_group))){
    log2FC<-log2(median(between_group))-log2(median(within_group));
    test<-wilcox.test(within_group,between_group);
    pvalue<-test$p.value;
    data_subset<-data.frame(cluster,log2FC,pvalue)
    log2FCdata<-rbind(log2FCdata,data_subset)
    }
# plot density line 
# manage data
    type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
    var<-c(within_group,between_group)
    data<-data.frame(type,var)
    data$type<-factor(data$type,levels=c("within-OR","between-OR"))
    p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
                  geom_density(alpha=.25) + theme_classic() 
    # t-test
    p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#raw counts heatmap 
# Heat map of expression  value for Or25-27 for cluster 5_1,5_2,14-1 and others( random select a few as control) 
    obj_barcode<-colnames(obj)
    #all_barcode<-colnames(ORN)
    #random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
    #obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
    #for (i in 1:length(obj$subcluster)){
    #  if(obj$subcluster[i]!=cluster){
    #        {obj$subcluster[i]="other"}
    #    }
    #}
    #DefaultAssay(obj)<-"SCT"
    obj_data<-as.data.frame(t(obj@assays$SCT[obj_features,]))
    obj_data<-obj_data[order(obj_data[,1],obj_data[,2],decreasing=T),]
    barcode_info<-data.frame(obj$subcluster)
    rownames(barcode_info)<-colnames(obj)
    p4<-pheatmap(t(obj_data),
             cluster_cols = T,
             cluster_rows = T,
             color = colorRampPalette(c("white", "red"))(100),
             annotation_col = barcode_info,
             #annotation_colors = ann_colors,
             #annotation_row = barcode_info,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
    # plot OR heatmap 
    #obj<-ScaleData(obj);
    #p4<-dittoHeatmap(obj,obj_features, annot.by = "subcluster",slot ="count",scaled.to.max = F)
    # merge 4 plots 
    top_right<-plot_grid(p2,p3,labels = c(" "," "),rel_widths = c(2, 1))
    right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" "," "))
    last<-plot_grid(p1$gtable, right, labels = c(' ', ''), label_size = 12, ncol = 2)
    title <- ggdraw() + 
      draw_label(
        paste("Cluster",cluster,":","log2FC=",log2FC),
        fontface = 'bold',
        x = 0,
        hjust = 0
      )
    add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
    print(add_title)
    }
dev.off()

# OR exp correlation in different cells 




# Extract result pages

library(pdftools)
# inputs
infile <- "./05_ORN_cluster2/03_coexp_cluster/distinguish_multi_OR_with_powerful.pdf"  # input pdf
outfile <-  "./05_ORN_cluster2/03_coexp_cluster/distinguish_multi_OR_with_powerful_extract.pdf"
num <- pdf_length(infile)/3
pdf_subset(infile, pages = (1:num)*3, output = outfile)


# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
pdf("./05_ORN_cluster2/03_coexp_cluster/multi_OR_without_nopower_trackplot.pdf",width=10,height=10)
for (cluster in clear_multiple_classes_order_id){
print(cluster)
obj<-subset(ORN,idents=cluster);
obj_features<-dotplot_data[dotplot_data$id==cluster,]$features.plot
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
p1<-CoveragePlot(
  object = obj,
  region = ranges.show,
  window = 150,
  extend.upstream = 100,
  annotation = TRUE,
  extend.downstream = 100,
  tile = TRUE,
  tile.size = 100,
  tile.cells = 30,
  links=F
)
print(p1)
}
dev.off()



cd /md01/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/Orco_pdb

ls *.pdb > pdb_list.txt
for file1 in $(<pdb_list.txt)
do
 for file2 in $(<pdb_list.txt)
 do 
  python ../RMSD.py $file1 $file2 >> RMSD_result.txt
 done
done

conda activate r4-base
R
RMSD <- read.table("RMSD_result.txt")
RMSD<- RMSD$V1
pdb <- read.table("pdb_list.txt")
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
data[data==0]<- 0.5
data<- log2(data)
pdf("./Orco_feature_RMSD.pdf",width=11,height=10)
pheatmap(data,
      cluster_cols = F,
      cluster_rows = F,
      color = colorRampPalette(c("#1C214F","#3082BD","#F1ECEC"))(100),
      #annotation_col = barcode_label_pheatmap,
      #annotation_colors = ann_colors,
      #annotation_row = barcode_label_pheatmap,
      annotation_legend = TRUE,
      show_rownames=T,
      show_colnames=T
 )
dev.off()