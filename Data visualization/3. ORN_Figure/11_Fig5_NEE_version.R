#Fig5_NEE_version:
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

ORN <- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_latest.rds")

#FigCDEF: single,RT,MP

#the distance of OR pair in the same cluster and random OR pair 
DefaultAssay(ORN)<-"ATAC"
gene_transcript<-Annotation(ORN)[which(Annotation(ORN)$type=="transcript"),];
OR_gene_transcript<-gene_transcript[gene_transcript$gene_name%in%dotplot_feature,]
OR_gene_transcript<-OR_gene_transcript[!duplicated(OR_gene_transcript$gene_name),]
AmelchrNameLength <- read.table("/md01/nieyg/ref/10X/honeybee/de_novo_antenna/Amel_antenan/star/chrNameLength.txt", sep="\t", stringsAsFactors=F) 
AmelchrNameLength<-AmelchrNameLength[which(AmelchrNameLength$V1%in%as.character(unique(OR_gene_transcript@seqnames@values))),]
data<-as.data.frame(sort(OR_gene_transcript));
data<-data[,c(1:3,12,4,5)]

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
        OR_pair$same_seqname[i]="same_seqname";}
        else{OR_pair$same_seqname[i]="diff_seqname"}
};

#Same cluster?
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$gene1;
    gene2<-OR_pair[i,]$gene2;
    cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id;
    if(length(which(duplicated(cluster_id)))){
        OR_pair$same_cluster[i]="co-exp_RT";
    }else{OR_pair$same_cluster[i]="single_OR";}
    
};

MP_gene1<- c("LOC102656904","LOC100577101","LOC107965761","LOC725052")
MP_gene2<- c("LOC102656221","Or41","LOC102655285","LOC100576839")
MP_OR_pair<- data.frame(gene1=MP_gene1,gene2=MP_gene2)

MP_gene<- c(MP_gene1,MP_gene2)
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$gene1;
    gene2<-OR_pair[i,]$gene2;
    if(gene1 %in%MP_gene){
    	if(gene2 %in% MP_gene){
    		if(OR_pair[i,6]=="co-exp_RT"){
    			OR_pair[i,6]="co-exp_MP"
    		}
    	}
    } 
};



# The probability of OR pair appearing on the same chromosome
# > table(OR_pair$same_seqname,OR_pair$same_cluster)
               co-exp_MP co-exp_RT single_OR
  diff_seqname         0         0      2341
  same_seqname         4        56       759

library(ggpubr)
probability<-data.frame(OR_pair_type=c("single_OR","co-exp_MP","co-exp_RT"),
    probability=c(759/3100,56/56,4/4),count=c(759,56,4));

colors_for_exp_pattern<- c('#2673A0','#E28118','#C23A39')
probability$OR_pair_type<-factor(probability$OR_pair_type,levels=c("single_OR","co-exp_MP","co-exp_RT"))
pdf("./00_Figure/Fig5/Fig5C-The_probability_OR_pair_appearing_on_same_chromosome_3groups.pdf",width=7,height=4)
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
same_chrom<-OR_pair[OR_pair$same_seqname=="same_seqname",]
for (i in 1:nrow(same_chrom)){
    gene1<-same_chrom[i,]$gene1;
    gene2<-same_chrom[i,]$gene2;
    OR_pair_transcript<-OR_gene_transcript[OR_gene_transcript$gene_name%in%c(gene1,gene2),]
    OR_pair_transcript<-sort(OR_pair_transcript);
    same_chrom$dist[i]<-width(gaps(OR_pair_transcript))[2];
}

same_chrom$same_cluster<-factor(same_chrom$same_cluster,levels=c("single_OR","co-exp_MP","co-exp_RT"))
my_comparisons <- list( c("single_OR", "co-exp_MP"),c("single_OR", "co-exp_RT") ,c("co-exp_MP", "co-exp_RT") )

b<- ggboxplot(same_chrom, x="same_cluster", y="dist",color = "same_cluster",width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("Genomic distance")+
scale_color_manual(values = colors_for_exp_pattern)
a|b 
dev.off()

#Fig 5E: Sequence similarity
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
OR_pair$Sequence_similarity<- NA
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$gene1;
    gene2<-OR_pair[i,]$gene2;
    OR_pair$Sequence_similarity[i]<-dist_percent[gene1,gene2];
}

OR_pair$same_cluster<-factor(OR_pair$same_cluster,levels=c("single_OR","co-exp_MP","co-exp_RT"))

pdf("./00_Figure/Fig5/Fig5E-OR_sequence_similarity_distribution_3groups.pdf",width=4,height=4)
p1<-ggplot(OR_pair, aes(x=Sequence_similarity, fill=same_cluster)) + xlab("% OR sequence similarity")+
          geom_density(alpha=.25) + theme_classic()+theme(legend.position="top")+
scale_color_manual(values =colors_for_exp_pattern)
p2<-ggboxplot(OR_pair, x="same_cluster", y="Sequence_similarity", color = "same_cluster",width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("% OR sequence similarity")+
scale_color_manual(values =colors_for_exp_pattern)
p2
dev.off()

# Fig2M: OR RMSD
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
OR_pair$RMSD<- NA
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$gene1;
    gene2<-OR_pair[i,]$gene2;
    if(gene1%in% ID){
    	if(gene2%in%ID){
    		OR_pair$RMSD[i]<-data[gene1,gene2];
    	}
    }
    
}
tmp_data<- na.omit(OR_pair)
pdf("./00_Figure/Fig5/Fig5F-OR_structure_similarity_distribution_3group.pdf",width=4,height=4)
ggboxplot(tmp_data, x="same_cluster", y="RMSD", color = "same_cluster",width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("OR RMSD")+
scale_color_manual(values = colors_for_exp_pattern)
dev.off()


# Fig5J Intergenic region length 
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']
gene.coords <- gene.coords[gene.coords$gene_name %in% all_receptor_gene]
gene.coords <- gene.coords[gene.coords$type == "transcript"]
gene.coords <- as.data.frame(gene.coords)
strand_p <- gene.coords[gene.coords$strand == "+",]
strand_n <- gene.coords[gene.coords$strand == "-",]
strand_p_gene <- strand_p$gene_name
strand_n_gene <- strand_n$gene_name

# Step1: Intergenic region length
MP_OR<- c(OR_pair[which(OR_pair$same_cluster=="co-exp_MP"),]$gene1,OR_pair[which(OR_pair$same_cluster=="co-exp_MP"),]$gene2)
RT_OR<- unique(c(OR_pair[which(OR_pair$same_cluster=="co-exp_RT"),]$gene1,OR_pair[which(OR_pair$same_cluster=="co-exp_RT"),]$gene2))
single_OR<- setdiff(dotplot_feature,c(MP_OR,RT_OR))

OR_type<- data.frame(OR=c(single_OR,MP_OR,RT_OR),type=c(rep("single_OR",length(single_OR)),rep("co-exp_MP",length(MP_OR)),rep("co-exp_RT",length(RT_OR))))
# for strand "+" OR gene 
strand_p <- strand_p[!duplicated(strand_p$gene_name),]
strand_p <- strand_p[order(strand_p$seqnames,strand_p$start,decreasing=F),]

tmp_data<- rbind(strand_p,strand_n)
tmp_data$exp_type<- OR_type[match(tmp_data$gene_name,OR_type$OR),]$type
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")

tmp_data$OR<- ORgene_name_trans[match(tmp_data$gene_name,ORgene_name_trans$gene_name),]$last_name
tmp_data<- tmp_data[which(!is.na(tmp_data$exp_type)),]
write.csv(tmp_data,"./OR_info.csv")


data<- read.csv("./OR_info.csv")
data<- na.omit(data)
data$dist[data$dist>10000]<- 10000;
library(ggpubr)
library(cowplot)
my_comparisons <- list( c("single_OR", "co-exp_MP"), c("single_OR", "co-exp_RT"), c("co-exp_MP", "co-exp_RT") )
data$exp_type<- factor(data$exp_type,levels=c("single_OR","co-exp_RT","co-exp_MP"))
pdf("./00_Figure/Fig5/Fig5J_singlevsNRTvsRT_Intergenic_region_length_distribution.pdf",width=4,height=4)

ggboxplot(data, x="exp_type", y="dist", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons, method='t.test')+theme(legend.position="none")+ylab("Intergenic region length")
dev.off()
pdf("./Intergenic_region_length_distribution.pdf",width=8,height=4)
ggplot(data, aes(dist)) +
  geom_density()+xlim(0,10000)
ggplot(data, aes(dist)) +
  geom_density()+xlim(0,20000)

dev.off()
# Fig5K: Number of "AAUAAA"
# trans to mRNA

# for strand "+" OR gene 
strand_p <- strand_p[!duplicated(strand_p$gene_name),]
strand_p <- strand_p[order(strand_p$seqnames,strand_p$start,decreasing=F),]
OR_TES_location<- data.frame()
for (i in 1:nrow(strand_p)){
	OR1 <- strand_p[i,]
	if(i!=nrow(strand_p)){
	OR2 <- strand_p[i+1,]
	if(OR1$seqnames==OR2$seqnames){
		OR1_TES_start <- OR1$end+1;
		intersections <- abs(OR1$end-OR2$start);
		OR1_TES_end <- OR2$start-1;
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_start,end=OR1_TES_end,strand=OR1$strand,intersections=intersections)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}
}
}
OR_TES_location_strand_p<- OR_TES_location;

# for strand "-" OR gene 
strand_n <- strand_n[!duplicated(strand_n$gene_name),]
strand_n <- strand_n[order(strand_n$seqnames,strand_n$start,decreasing=T),]

OR_TES_location<- data.frame()
for (i in 1:nrow(strand_n)){
	OR1 <- strand_n[i,]
	if(i!=nrow(strand_n)){
	OR2 <- strand_n[i+1,]
	if(OR1$seqnames==OR2$seqnames){
		OR1_TES_start <- OR1$start-1;
		intersections <- abs(OR1$start-OR2$end);
		OR1_TES_end <- OR2$end+1;
		OR1_TES_location <- data.frame(OR=OR1$gene_name,seqnames=OR1$seqnames,start=OR1_TES_end,end=OR1_TES_start,strand=OR1$strand,intersections=intersections)
		OR_TES_location <- rbind(OR_TES_location,OR1_TES_location)
	}
}
}
OR_TES_location$strand<- "-"
OR_TES_location<- rbind(OR_TES_location,OR_TES_location_strand_p)
OR_TES_location$seqnames<- as.character(OR_TES_location$seqnames)

# only keep the OR in Fig5J
OR_TES_location$OR_symbol<- ORgene_name_trans[match(OR_TES_location$OR,ORgene_name_trans$gene_name),]$last_name
OR_TES_location<- OR_TES_location[which(OR_TES_location$OR_symbol%in% data$OR),]

for(i in 1:nrow(OR_TES_location)){
	OR_TES_location_subset <- OR_TES_location[i,];
	#OR_TES_location_subset$V1<- 1;
	OR<- OR_TES_location_subset$OR_symbol
	write.table(OR_TES_location_subset[,c(2:4,7,5)],paste0("./13_TES_signal/06_Fig5_AAUAAA/",OR,".bed"),row.names=F,col.names=F,sep="\t")
}
write.table(OR_TES_location[,c(2:4,7,5)],"./13_TES_signal/06_Fig5_AAUAAA/OR_TES_location.bed",row.names=F,col.names=F,sep="\t")


# back to shell 
sed -i 's/"//g' *.bed

ls *.bed > OR_promoter.list
for file1 in $(<OR_promoter.list)
do
  bedtools getfasta -fi /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/fasta/genome.fa -bed $file1 -fo ./$file1.fa -s -name
done

cat Or*bed.fa > OR_intergenic_region.fa

library(Biostrings)
library(seqinr)
fasta_file <- "/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/06_Fig5_AAUAAA/OR_intergenic_region.fa"
OR_fasta<-readDNAStringSet(fasta_file, format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

count_polya <- function(dna_sequence) {
  rna_sequence <- RNAString(dna_sequence)
  count <- countPattern("AAUAAA", rna_sequence)
  return(count)
}
polya_counts <- sapply(OR_fasta, function(seq) count_polya(seq))

# 打印结果
result <- data.frame(ID = names(OR_fasta),
                     Polya_Count = polya_counts)

result$ID<- gsub("::.*","",result$ID)
data<- read.csv("./OR_info.csv")
data<- na.omit(data)

data$Polya_Count<- result[match(data$OR,result$ID),]$Polya_Count


library(ggpubr)
library(cowplot)
my_comparisons <- list( c("single_OR", "co-exp_MP"), c("single_OR", "co-exp_RT"), c("co-exp_MP", "co-exp_RT") )
data$exp_type<- factor(data$exp_type,levels=c("single_OR","co-exp_RT","co-exp_MP"))
pdf("./00_Figure/Fig5/Fig5K_singlevsNRTvsRT_Intergenic_region_AAUAAA.pdf",width=4,height=4)
data$Polya_Count[data$Polya_Count>20]<- 20;

ggboxplot(data, x="exp_type", y="Polya_Count", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons, method='t.test')+
theme(legend.position="none")+ylab("The number of AAUAAA")
dev.off()

# The ration of Polya_Count and Intergenic_region_length
data<- read.csv("./OR_info.csv")
data<- na.omit(data)
data$Polya_Count<- result[match(data$OR,result$ID),]$Polya_Count
data$ratio<- data$Polya_Count/data$dist

pdf("./00_Figure/Fig5/Fig5L_singlevsNRTvsRT_ratio_of_AAUAAA_and_Intergenic_region.pdf",width=4,height=4)
#data$Polya_Count[data$Polya_Count>20]<- 20;
data$exp_type<- factor(data$exp_type,levels=c("single_OR","co-exp_RT","co-exp_MP"))

ggboxplot(data, x="exp_type", y="ratio", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons, method='t.test')+
theme(legend.position="none")+ylab("The number of AAUAAA/intergenic length")
dev.off()




Number of “AAUAAA”/intergenic length

