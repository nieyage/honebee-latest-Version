# Fig5A:

# step1: plot the OR ggtree 
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
pub<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep_tree.aa")

tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")

#OR2 is placed in the last column;
all_OR_gene_fasta<- c(pub,supply_fasta[c(1,4,5)])

all_OR_gene_fasta<- all_OR_gene_fasta[which(names(all_OR_gene_fasta)%in% tree_anno$OR_gene),]
names(all_OR_gene_fasta) <- tree_anno$last_name[match(names(all_OR_gene_fasta),tree_anno$OR_gene)]

aln <- muscle::muscle(all_OR_gene_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="ward.D")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)

library(ggtreeExtra)
library(ggplot2)
library(ggtree)
library(treeio)
library(ggnewscale)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', 
                    '#E0D4CA', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', 
                    '#6778AE', '#B53E2B', '#DCC1DD', '#CCE0F5', '#625D9E', 
                    '#68A180', '#968175', '#FF9999', '#344CB7', '#FFCC1D', 
                    '#24A19C', '#FF9A76',"#BC80BD", "#CCEBC5", "#FFED6F", 
                    "#E95C59", "#476D87",'#9FA3A8')


tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")
OR_target<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_iORdb_pubChem.csv")

# seqnames
tree_anno$Pherobase<- OR_target$Pherobase[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$origin<- OR_target$origin[match(tree_anno$OR_gene,OR_target$OR_gene)]
tree_anno$origin[is.na(tree_anno$origin)]<-"undefine"
tree_anno$Pherobase[is.na(tree_anno$Pherobase)]<-"undefine"
tree_anno$Pherobase[tree_anno$Pherobase=="-"]<-"undefine"

Group =    c("Group1"=myUmapcolors[1],    "Group2"=myUmapcolors[2],    "Group3"=myUmapcolors[3],    "Group4"=myUmapcolors[4],    "Group5"=myUmapcolors[5],    "Group6"=myUmapcolors[6],    "Group7"=myUmapcolors[7],    "Group8"=myUmapcolors[8],    "Group9"=myUmapcolors[9],    "Group10"=myUmapcolors[10],
    "Group11"=myUmapcolors[11],    "Group12"=myUmapcolors[12],    "Group13"=myUmapcolors[13],    "Group14"=myUmapcolors[14],    "Group15"=myUmapcolors[15],    "Group16"=myUmapcolors[16],    "GroupUN3"=myUmapcolors[17],
    "GroupUN226"=myUmapcolors[18],    "GroupUN243"=myUmapcolors[19],    "GroupUN248"=myUmapcolors[20])
rownames(tree_anno)<- tree_anno$last_name
data<-  as.data.frame(tree_anno[,c(5,6,9,15)])
data_long<- melt(data, id.vars = c("last_name"), #需保留的不参与聚合的变量列名
                  measure.vars = c('seqnames','exp_pattern3','Pherobase'),
                  variable.name = c('POS'),#聚合变量的新列名
                  value.name = 'value')
order<- c("Group1","Group2" ,"Group4","Group5","Group7","Group9","Group10",
"Group11","Group12","Group13","Group14","Group15","Group16","GroupUN243","GroupUN248",
"single_OR","MP","SP",unique(tree_anno$Pherobase)[-2],
"undefine")
data_long$value<- factor(data_long$value,levels=order)
pdf("./00_Figure/Fig5/Fig5A-OR_sequence_protein_similarity-tree_add_anno.pdf",width=15,height=15)
clust <- hclust(sdist,method="ward.D")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)
p1<- ggtree(tree, layout="fan",size=0.5) + geom_tiplab(size=3)
p2<- p1+ new_scale_fill() + 
      geom_fruit(
          data=data_long,
          geom=geom_tile,
          mapping=aes(y=last_name, x=POS,fill=value),
          offset=0.2,   # The distance between external layers, default is 0.03 times of x range of tree.
          pwidth=0.1 # width of the external layer, default is 0.2 times of x range of tree.
      ) +
      scale_fill_manual(
          values=c(myUmapcolors[1:15],"#129FAF", "#DE7C5B","#FBD277",myUmapcolors[16:23],"grey"),
          guide=guide_legend(keywidth=1, keyheight=1, order=3)
      ) 
p2
dev.off()


# gene type for dotplot OR 
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
dotplot_data<-read.csv("../05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
coexp_OR <- unique(dotplot_data[dotplot_data$id%in% multiOR_cluster,]$features.plot)
data<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/13_TES_signal/dotplot_feature_TES_type.csv")
RT_OR<- data[data$type=="RT",2]
dotplot_data$gene_type <- ifelse(dotplot_data$features.plot%in% coexp_OR,"coexp(NRT)","single_OR")
dotplot_data[which(dotplot_data$features.plot%in% RT_OR),]$gene_type="coexp(RT)"

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

# for strand "+" OR gene 
strand_p <- strand_p[!duplicated(strand_p$gene_name),]
strand_p <- strand_p[order(strand_p$seqnames,strand_p$start,decreasing=F),]
strand_p_intersections<- data.frame()
for (i in 1:nrow(strand_p)){
  OR1 <- strand_p[i,]
  if(i!=nrow(strand_p)){
  OR2 <- strand_p[i+1,]
  if(OR1$seqnames==OR2$seqnames){
    intersections <- OR2$start-OR1$end
    tmp<- data.frame(OR=OR1$gene_name,intersections=intersections)
    strand_p_intersections<- rbind(strand_p_intersections,tmp)
  }
  }
}
# for strand "-" OR gene 
strand_n <- strand_n[!duplicated(strand_n$gene_name),]
strand_n <- strand_n[order(strand_n$seqnames,strand_n$start,decreasing=T),]
strand_n_intersections<- data.frame()
for (i in 1:nrow(strand_n)){
  OR1 <- strand_n[i,]
  if(i!=nrow(strand_n)){
  OR2 <- strand_n[i+1,]
  if(OR1$seqnames==OR2$seqnames){
    intersections <- -(OR2$start-OR1$end)
    tmp<- data.frame(OR=OR1$gene_name,intersections=intersections)
    strand_n_intersections<- rbind(strand_n_intersections,tmp)
  }
}
}

intersections_last<- rbind(strand_n_intersections,strand_p_intersections)
intersections_last$type <- tree_anno[match(intersections_last$OR,tree_anno$OR_gene),]$exp_pattern3

intersections_last <- na.omit(intersections_last)
intersections_last$OR_symbol<- ORgene_name_trans[match(intersections_last$OR,ORgene_name_trans$gene_name),]$last_name
intersections_last<- intersections_last[intersections_last$type!="undefine",]

intersections_last$type<-factor(intersections_last$type,levels=c("single_OR","SP","MP"))
library(ggpubr)
library(cowplot)
my_comparisons <- list( c("single_OR", "MP"), c("single_OR", "SP"), c("MP", "SP") )

pdf("./00_Figure/Fig5/singlevsNRTvsRT_Intergenic_region_length_distribution.pdf",width=6,height=3)
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,2000)+
stat_compare_means(comparisons = my_comparisons, label.y = c(18, 22, 26))+theme(legend.position="none")+ylab("Intergenic region length")
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,20000)+
stat_compare_means(comparisons = my_comparisons, label.y = c(18, 22, 26))+theme(legend.position="none")+ylab("Intergenic region length")
ggboxplot(intersections_last, x="type", y="intersections", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons)+theme(legend.position="none")+ylab("Intergenic region length")
dev.off()



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
data<- intersections_last
data<- na.omit(data)

data$Polya_Count<- result[match(data$OR_symbol,result$ID),]$Polya_Count


library(ggpubr)
library(cowplot)
my_comparisons <- list( c("single_OR", "MP"), c("single_OR", "SP"), c("MP", "SP") )
data$type<- factor(data$type,levels=c("single_OR","SP","MP"))
pdf("./00_Figure/Fig5/Fig5K_singlevsNRTvsRT_Intergenic_region_AAUAAA.pdf",width=4,height=4)
data$Polya_Count[data$Polya_Count>20]<- 20;
ggboxplot(data, x="type", y="Polya_Count", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons, method='t.test')+
theme(legend.position="none")+ylab("The number of AAUAAA")
dev.off()

# The ration of Polya_Count and Intergenic_region_length
data<- intersections_last

data$Polya_Count<- result[match(data$OR_symbol,result$ID),]$Polya_Count
data$type<- factor(data$type,levels=c("single_OR","SP","MP"))
data$ratio<- data$Polya_Count/abs(data$intersections)
data<- na.omit(data)
pdf("./00_Figure/Fig5/singlevsNRTvsRT_ratio_of_AAUAAA_and_Intergenic_region.pdf",width=4,height=4)
#data$Polya_Count[data$Polya_Count>20]<- 20;
ggboxplot(data, x="type", y="ratio", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons, method='t.test')+
theme(legend.position="none")+ylab("The number of AAUAAA/intergenic length")
dev.off()







