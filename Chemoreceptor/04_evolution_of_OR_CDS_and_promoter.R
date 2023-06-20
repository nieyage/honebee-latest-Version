# Step1: calculate the matrix of all OR CDS

# the matrix in OR CDS region 
# a. the distance between OR 
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
all_receptor_gene_fasta<- c(OR_fasta,GR_fasta,IR_fasta,supply_fasta)
# corelation heatmap tree
dotplot_feature_fa<- all_receptor_gene_fasta
aln <- muscle::muscle(dotplot_feature_fa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)
sequence_dist<- as.matrix(sdist)
sequence_data<- melt(sequence_dist);
sequence_data$OR_pair<- paste(sequence_data$Var1,sequence_data$Var2,sep="-");

# b. dN matrix 
# c. dS matrix
# get the CDS sequence 
cd /md01/nieyg/ref/10X/Amel_HAv3.1
seqkit grep -f OR_transcript.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa -o OR_transcript_cds.fa 
seqkit grep -f GR_transcript.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa -o GR_transcript_cds.fa 
seqkit grep -f IR_transcript.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa -o IR_transcript_cds.fa 
seqkit grep -f supply.txt GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid_cds.fa -o supply_cds.fa 
#ID trans 
```
nLine=`wc -l chemoreceptor_IDtrans.list | cut -d " " -f 1`
for i in `seq 1 ${nLine}`;
do
        transcript=`awk '{print $1}' chemoreceptor_IDtrans.list | sed -n "${i}p"`
        gene_name=`awk '{print $2}' chemoreceptor_IDtrans.list | sed -n "${i}p"`
        sed -i "s/${transcript}/${gene_name}/g" OR_transcript_cds.fa 
        sed -i "s/${transcript}/${gene_name}/g" GR_transcript_cds.fa 
        sed -i "s/${transcript}/${gene_name}/g" IR_transcript_cds.fa 
        sed -i "s/${transcript}/${gene_name}/g" supply_cds.fa 
done
``` 
cat OR_transcript_pep.aa GR_transcript_pep.aa IR_transcript_pep.aa supply.aa >all_receptor_gene_pep.aa
cat OR_transcript_cds.fa GR_transcript_cds.fa IR_transcript_cds.fa supply_cds.fa >all_receptor_gene_cds.fa

# parpare the OR pair 
library(seqinr)
cds <- read.alignment(file ="/md01/nieyg/ref/10X/Amel_HAv3.1/all_receptor_gene_cds.fa" ,format = "fasta")
OR_gene<- cds$nam
OR_pair<- data.frame()
for(i in 1:length(OR_gene)){
	for (j in 1:length(OR_gene)) {
		if(i!=j){
			data_subset<- data.frame(gene1=OR_gene[i],gene2=OR_gene[j])
			OR_pair<- rbind(OR_pair,data_subset)
		}
	}
}
write.csv(OR_pair,"/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/OR_pair.homologs",row.names=F,col.names=F)

sed -i 's/"//g' OR_pair.homologs
conda activate kaks_calculate

ParaAT.pl -h /md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/OR_pair.homologs \
 -n /md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/all_receptor_gene_cds.fa \
 -a /md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/all_receptor_gene_pep.aa \
 -p proc -m mafft -f axt -g -k -o /md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/result_dir

awk '!(NR%2)' merge.kaks > OR_pair_kaks.out

kaks_data<- read.table("/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/result_dir3/kaks/OR_pair_kaks.out")
colnames(kaks_data)<- c("OR_pair","method","KA","Ks","Ka/Ks","p_value","length","S-sites","N-sites","Folo-sites","Substitutions","S-sub","N-sub","Folo-S-sub","Folo-N-sub","Divergence-Time","Substitution-Rate-Ratio","GC","ML_score","AICc","Akaike-Weight","Model")
dN_data <- kaks_data[,c(1,3)]
dS_data <- kaks_data[,c(1,4)]


# d. structure matrix 
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
RMSD_matrix<- data;
RMSD_data<- melt(RMSD_matrix);
RMSD_data$OR_pair<- paste(RMSD_data$Var1,RMSD_data$Var2,sep="-");

# Step2: calculate the martrix of all OR promoter 
# A. Kmer similarity (k=5)
library(stringr)
a=read.csv('/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_fa/matrix.count',header=F,sep="\t")
colnames(a)=c('OR','kmer','count')
library(reshape2)
counts=dcast(a,formula=kmer~OR)
counts[is.na(counts)]=0
rownames(counts)<- counts[,1]
counts=counts[,-1]
# k=5 
k_5_data <- counts[which(str_length(rownames(counts))==5),]
# correlation coefficient 
OR_kmer_cor <- cor(k_5_data,method="pearson")
# distance 
OR_kmer_dist <- as.matrix(dist(t(k_5_data), method = "euclidean"))

# B. motif similarity 
motif_data<- read.csv("/data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/TOBIAS_motif/Or_motif_bindscore.csv")
motif_data<- dcast(motif_data, MotifID ~ symbol, value.var="Bindscore")
motif_data[is.na(motif_data)]= -300

# cor
OR_motif_cor = as.matrix(cor(motif_data[,-1]))

# OR_motif_braycurtis
library(vegan)
OR_motif_braycurtis = as.matrix(vegdist(t(motif_data[,-1]), "bray"))

# Step3: the relationship between the CDS and the promoter 

# A. all OR pair 

# pairwised compared -->R squard
library(reshape2)
colnames(dN_data)<- c("OR_pair","value")
colnames(dS_data)<- c("OR_pair","value")

CDS_list <- list(sequence_data,dN_data,dS_data,RMSD_data)
promoter_list <- list(OR_kmer_cor,OR_kmer_dist,OR_motif_cor,OR_motif_braycurtis)
names(CDS_list)<- c("distance","dN","dS","RMSD")
names(promoter_list)<- c("kmer_cor","kmer_dist","motif_cor","motif_dist")
r2_data<- matrix(ncol=4,nrow=4)
rownames(r2_data)<- c("distance","dN","dS","RMSD")
colnames(r2_data)<- c("kmer_cor","kmer_dist","motif_cor","motif_dist")

library(ggpmisc)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/All_OR_pair_cor.pdf",width=6,height=4)
for (CDS_ID in 1:4){
	for (promoter_ID in 1:4){
		CDS_matrix<- CDS_list[[CDS_ID]];
		promoter_matrix<- promoter_list[[promoter_ID]];
		#CDS_matrix<- melt(CDS_matrix)
		promoter_matrix<- melt(promoter_matrix);
		#CDS_matrix$OR_pair<- paste(CDS_matrix$Var1,CDS_matrix$Var2,sep="-")
		promoter_matrix$OR_pair<- paste(promoter_matrix$Var1,promoter_matrix$Var2,sep="-");
		data<- merge(CDS_matrix,promoter_matrix,by="OR_pair")
		data<- data[,c("OR_pair","value.x","value.y")]
        colnames(data)<- c("OR_pair","CDS","promoter")
        data<- data[!duplicated(data),]
        r2<- cor(data$CDS,data$promoter)
        r2_data[CDS_ID,promoter_ID]<- r2;
        p<- ggplot(data, aes(CDS,promoter))+ 
            geom_point()+ 
            geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
            theme_classic()+ggtitle(paste0(names(CDS_list)[CDS_ID],"_",names(promoter_list)[promoter_ID]))+
            stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
            print(p)
	}
}
dev.off()
library(pheatmap)
pdf("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/CDS_promoter_cor_heatmap.pdf",width=8,height=8)
pheatmap(r2_data,
      cluster_cols = T,
      cluster_rows = T,
      color = colorRampPalette(c("#0E2954","#9BCDD2","#FAF0E4"))(100),
      cellwidth=10,cellheight=10,
      show_rownames=T,
      show_colnames=T
 )
dev.off()

# B. single exp OR pair and co-exp OR

dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
Freq<- as.data.frame(table(dotplot_data$id));
multiple_id <- as.character(Freq[Freq$Freq>1,1])
single_id <- as.character(Freq[Freq$Freq==1,1])
single_OR<- dotplot_data[dotplot_data$id%in%single_id,]$features.plot;
coexp_OR<- dotplot_data[dotplot_data$id%in%multiple_id,]$features.plot;

# single OR 
# pairwised compared -->R squard
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/Single_OR_pair_cor.pdf",width=6,height=4)
for (CDS_ID in 1:4){
	for (promoter_ID in 1:4){
		CDS_matrix<- CDS_list[[CDS_ID]];
		promoter_matrix<- promoter_list[[promoter_ID]];
		#CDS_matrix<- melt(CDS_matrix)
		promoter_matrix<- melt(promoter_matrix);
		#CDS_matrix$OR_pair<- paste(CDS_matrix$Var1,CDS_matrix$Var2,sep="-")
		promoter_matrix$OR_pair<- paste(promoter_matrix$Var1,promoter_matrix$Var2,sep="-");
		data<- merge(CDS_matrix,promoter_matrix,by="OR_pair")
		data<- data[data$Var1%in% single_OR,]
		data<- data[data$Var2%in% single_OR,]
		data<- data[,c("OR_pair","value.x","value.y")]
        colnames(data)<- c("OR_pair","CDS","promoter")
        data<- data[!duplicated(data),]
        r2<- cor(data$CDS,data$promoter)
        r2_data[CDS_ID,promoter_ID]<- r2;
        p<- ggplot(data, aes(CDS,promoter))+ 
            geom_point()+ 
            geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
            theme_classic()+ggtitle(paste0(names(CDS_list)[CDS_ID],"_",names(promoter_list)[promoter_ID]))+
            stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
        print(p)
	}
}
dev.off();

pdf("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/single_OR_CDS_promoter_cor_heatmap.pdf",width=8,height=8)
pheatmap(r2_data,
      cluster_cols = T,
      cluster_rows = T,
      color = colorRampPalette(c("#0E2954","#9BCDD2","#FAF0E4"))(100),
      cellwidth=10,cellheight=10,
      show_rownames=T,
      show_colnames=T
 )
dev.off()

# single OR 
# pairwised compared -->R squard
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/coexp_OR_pair_cor.pdf",width=6,height=4)
for (CDS_ID in 1:4){
	for (promoter_ID in 1:4){
		CDS_matrix<- CDS_list[[CDS_ID]];
		promoter_matrix<- promoter_list[[promoter_ID]];
		#CDS_matrix<- melt(CDS_matrix)
		promoter_matrix<- melt(promoter_matrix);
		#CDS_matrix$OR_pair<- paste(CDS_matrix$Var1,CDS_matrix$Var2,sep="-")
		promoter_matrix$OR_pair<- paste(promoter_matrix$Var1,promoter_matrix$Var2,sep="-");
		data<- merge(CDS_matrix,promoter_matrix,by="OR_pair")
		data<- data[data$Var1%in% coexp_OR,]
		data<- data[data$Var2%in% coexp_OR,]
		data<- data[,c("OR_pair","value.x","value.y")]
        colnames(data)<- c("OR_pair","CDS","promoter")
        data<- data[!duplicated(data),]
        r2<- cor(data$CDS,data$promoter)
        r2_data[CDS_ID,promoter_ID]<- r2;
        p<- ggplot(data, aes(CDS,promoter))+ 
            geom_point()+ 
            geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
            theme_classic()+ggtitle(paste0(names(CDS_list)[CDS_ID],"_",names(promoter_list)[promoter_ID]))+
            stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=2.5)
        print(p)
	}
}
dev.off();
pdf("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/coexp_OR_CDS_promoter_cor_heatmap.pdf",width=8,height=8)
pheatmap(r2_data,
      cluster_cols = T,
      cluster_rows = T,
      color = colorRampPalette(c("#0E2954","#9BCDD2","#FAF0E4"))(100),
      cellwidth=10,cellheight=10,
      show_rownames=T,
      show_colnames=T
 )
dev.off()