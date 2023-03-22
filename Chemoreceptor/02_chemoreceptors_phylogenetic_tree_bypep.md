## phylogenetic tree only honeybee 

1. OR 
# New gtf OR sequence similarity:

```
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(OR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
chemoreceptor_info_data<-read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=T)
chemoreceptor_info_data$gene_name <- make.unique(chemoreceptor_info_data$gene_name, sep = "_")
Group_info <- unique(chemoreceptor_info_data[,c(1,12)])
rownames(Group_info)<-Group_info$gene_name
groupInfo <- split(row.names(Group_info), Group_info$seqnames)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=16,height=16)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()

# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/all_OR_fasta_aa_similarity.csv");

```

2. GR 

```
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(GR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/GR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=10,height=10)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()

# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/GR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/all_GR_fasta_aa_similarity.csv");

```

3. IR 


```
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
aln <- muscle::muscle(IR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/IR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=10,height=10)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()
# get the similarity rate matrix 
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
dist_percent[dist_percent==100.00]<-0;

data<-data.frame()
for (i in 1:nrow(dist_percent)){
	j=i
	row=i;
	col=j+1;
	if(dist_percent[row,col]>50){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
write.csv(data,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/IR_fasta_aa_similarity_50.csv");
write.csv(dist_percent,"/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/all_IR_fasta_aa_similarity.csv");

```

## phylogenetic tree honeybee, fly and  

1. get sequence 

```
sed -i '/#/d;' AaegL5-ChemoreceptorPeptides.fasta 
```



```
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
AaegL5_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/AaegL5-ChemoreceptorPeptides.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

#select names from your XStringSet
names_in_myXStringSet <- names(AaegL5_fasta)
ORnames_in_myXStringSet <- as.data.frame(grep("^AgamOr",names_in_myXStringSet))
IRnames_in_myXStringSet <- as.data.frame(grep("^AgamIr",names_in_myXStringSet))
GRnames_in_myXStringSet <- as.data.frame(grep("^AgamGr",names_in_myXStringSet))

colnames(ORnames_in_myXStringSet ) <- "Names"
colnames(IRnames_in_myXStringSet ) <- "Names"
colnames(GRnames_in_myXStringSet ) <- "Names"

#Extract complete sequence names and make a new list
final.list <- NULL
for (i in 1:nrow(ORnames_in_myXStringSet)) {
  sel0 <- ORnames_in_myXStringSet$name[[i]] #assuming your list has "name" as the column name
  sel1 <- dplyr::filter(names_in_myXStringSet , grepl(sel0, Names))
  final.list <- rbind(final.list, sel1)
}



#get fly OR sequence;
library(biomaRt)
mart <- useMart("ensembl","dmelanogaster_gene_ensembl")
fly<-readRDS("/md01/liyh526/project/Fly-cooperate/7.18run/WNN_ORN_integrated_antenna.rds")
allgene<-rownames(fly)
OR<-grep("^Or[0-9]",allgene,value=T)
IR<-grep("^Ir[0-9]",allgene,value=T)
GR<-grep("^Gr[0-9]",allgene,value=T)

all<-c(OR,IR,GR,"Orco")
all_ppseqs <- getSequence(id = all,
                          type="external_gene_name",
                          seqType="peptide",
                          mart = mart) 
allppseqs <- all_ppseqs[all_ppseqs$peptide!="Sequence unavailable",]
index <- duplicated(allppseqs$external_gene_name)
i <- 1
for(j in 1:nrow(allppseqs)){
if(index[j]==FALSE){
  i = 1
  allppseqs[j,2] <- paste0(allppseqs[j,2], ".",i)
  }
  else {
    i = i+1
    allppseqs[j,2] <- paste0(allppseqs[j,2], ".",i)
  } 
}
exportFASTA(allppseqs,"./ORN/sequence_tree/fly/fly_OR.fasta")

All_Amel_fasta<- c(OR_fasta.GR_fasta,IR_fasta)
fly_OR_fasta<-readAAStringSet("./ORN/sequence_tree/fly/fly_OR.fasta", format="fasta",
                          nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)




aln <- muscle::muscle(fly_OR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)



aln <- muscle::muscle(IR_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
tree <- groupOTU(tree, groupInfo)

```