## get infomation from the iorbase
1. get the Apis OR sequence and mapping with our dotplot feature 

```
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
dotplot_data<-read.csv("~/project/honeybee/honebee-latest-Version//05_ORN_cluster/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));

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
dotplot_feature_fa<- all_receptor_gene_fasta[dotplot_feature_OR2,]
iOR_fasta<-readAAStringSet("/data/R02/nieyg/project/honeybee/honebee-latest-Version/06_iOR_database/iORbase_ApMe_seq.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
total<- c(dotplot_feature_fa,iOR_fasta)

aln <- muscle::muscle(total)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist, method = "single")
tree<-as.phylo(clust)
pdf("./iOR_protein_mapping_our_Data.pdf",width=16,height=16)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none") + geom_tiplab2(size=3) + theme(legend.position = "right")
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
	if(dist_percent[row,col]>5){
		row_name<-rownames(dist_percent)[row];
		col_name<-rownames(dist_percent)[col];
		if(row_name!=col_name){
			data_subset<-data.frame(gene1=row_name,gene2=col_name,rate=dist_percent[row,col]);
			data<- rbind(data,data_subset)
			}
	}
}
data<-data[order(data$rate,decreasing=TRUE),]
```

2. get the Apis OR pdb file 

```
#!/bin/bash

cat iORbase_id.txt | while read line
do
#echo $line
arr=($line)
id=${arr[0]}
echo  https://www.iorbase.com/protein_download/?change_name=$id -O $id.pdb >> url.txt
done

```
# %s/\r//g modify url.txt
sed -i -e 's/^/wget /' url.txt
bash url.txt

3. calculate the structural similarity by RMSD 
# You have molecule A and B and want to calculate the structural difference between those two. 
# If you just calculate the RMSD straight-forward you might get a too big of a value as seen below. 
# You would need to first recenter the two molecules and then rotate them unto each other to get the true minimal RMSD. 
# This is what this script does.
# 1) the basical model  
/md01/nieyg/software/calculate_rmsd --no-hydrogen --print ./pdb_file/ApMe_9.pdb ./pdb_file/ApMe_91.pdb
# 2) ignore all hydrogens model: (useful for larger molecules where hydrogens move around indistinguishable) print the rotated structure for visual comparison.
calculate_rmsd --no-hydrogen --print tests/ethane.xyz tests/ethane_mini.xyz
# 3) If the atoms are scrambled and not aligned you can use the --reorder argument which will align the atoms from structure B unto A. Use --reorder-method to select what method for reordering. Choose between Hungarian (default), distance (very approximate) and brute force (slow).
calculate_rmsd --reorder tests/water_16.xyz tests/water_16_idx.xyz



