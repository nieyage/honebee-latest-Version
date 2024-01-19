# A dot plot of tandem duplication at group2 in the honeybee reference genome detected by CNVnator.
conda activate cnvnator 

# single cell pseudobam file 
single_bam<- "/md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam"


# step1: EXTRACTING READ MAPPING FROM BAM/SAM FILES
cnvnator -root honeybee_singlecell.root -chrom Group2 -tree /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam

# step2: GENERATING A READ DEPTH HISTOGRAM
awk '/^>/{OUT=substr($0, 2) ".fa"} {print > OUT}' genome.fa

cnvnator -root honeybee_singlecell.root -chrom Group2 -his 10 -d genome/

# step3:CALCULATING STATISTICS
cnvnator -root honeybee_singlecell.root -chrom Group2 -stat 10

# step4:RD SIGNAL PARTITIONING
#cnvnator -root honeybee_singlecell.root -chrom Group2 -partition 100 
# Option -ngc specifies not to use GC corrected RD signal. Partitioning is the most time consuming step.
cnvnator -root honeybee_singlecell.root -chrom Group2 -partition 10  

# step5:CNV CALLING
cnvnator -root honeybee_singlecell.root -chrom Group2 -call 10 > honeybee_singlecell_ngc.txt

# step6:REPORTING READ SUPPORT
cnvnator -root honeybee_singlecell.root -pe /md01/nieyg/project/honeybee/honebee-latest-Version/08_ORN_cluster_bw/merged.bam  -qual val(20) -over val(0.8) [-f file]


##### bulk RNAseq ######  

nohup samtools merge honeybee_bulk_merged.bam D1Aligned.sortedByCoord.out.bam  F3Aligned.sortedByCoord.out.bam   NE2Aligned.sortedByCoord.out.bam  N1Aligned.sortedByCoord.out.bam   NE3Aligned.sortedByCoord.out.bam N2Aligned.sortedByCoord.out.bam   Q2Aligned.sortedByCoord.out.bam  F1Aligned.sortedByCoord.out.bam   N3Aligned.sortedByCoord.out.bam   Q3Aligned.sortedByCoord.out.bam  F2Aligned.sortedByCoord.out.bam   NE1Aligned.sortedByCoord.out.bam &

cnvnator -root honeybee_bulk_merged.root -chrom Group2 -tree /md01/nieyg/project/honeybee/bulk-RNAseq/03_STAR_remapping/honeybee_bulk_merged.bam

# step2: GENERATING A READ DEPTH HISTOGRAM
cnvnator -root honeybee_bulk_merged.root -chrom Group2 -his 100 -fasta ../genome.fa.gz # -d genome/

# step3:CALCULATING STATISTICS
cnvnator -root honeybee_bulk_merged.root -chrom Group2 -stat 100

# step4:RD SIGNAL PARTITIONING
#cnvnator -root honeybee_bulk_merged.root -chrom Group2 -partition 100 
# Option -ngc specifies not to use GC corrected RD signal. Partitioning is the most time consuming step.
cnvnator -root honeybee_bulk_merged.root -chrom Group2 -partition 100 -ngc  

# step5:CNV CALLING
cnvnator -root honeybee_bulk_merged.root -ngc -chrom Group2 -call 100 > honeybee_bulk_ngc.txt

# step6:REPORTING READ SUPPORT
cnvnator -root honeybee_bulk_merged.root -pe /md01/nieyg/project/honeybee/bulk-RNAseq/03_STAR_remapping/honeybee_bulk_merged.bam -qual val(20) -over val(0.8) 


Group2:8820521-14609630


cnvnator -root honeybee_bulk_merged.root -genotype 100 

cnvnator -root honeybee_bulk_merged.root  -view 100 

cnvnator plotbaf.py honeybee_bulk_merged.root Group2:8820521-14609630


# part2: The dN/dS in honeybee OR 

# Step1: parpare the OR gene cds and pep file (honeybee and ant)

honeybee_cds=/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/all_receptor_gene_cds.fa
honeybee_pep=/md01/nieyg/ref/10X/Amel_HAv3.1/kaks_calculate/all_receptor_gene_pep.aa
cd /md01/nieyg/project/honeybee/honebee-latest-Version/14_OR_evolution/select_pressure/06_honeybee_OR
cp $honeybee_pep ./
cp $honeybee_cds ./

awk '!/^>.*_[0-9]+$/' all_receptor_gene_pep.aa > all_receptor_gene_pep_uniq.aa
awk '!/^>.*_[0-9]+$/' all_receptor_gene_cds.fa > all_receptor_gene_cds_uniq.fa

makeblastdb -in all_receptor_gene_pep_uniq.aa -dbtype prot
###使用blast进行alignment，得到一个表格，输出格式为m6
blastp -query all_receptor_gene_pep_uniq.aa -db all_receptor_gene_pep_uniq.aa -evalue 1e-5 -max_target_seqs 2 -num_threads 2 -out honeybee_blastp_out.m6 -outfmt 6
###而后使用shell指令进行简单的排序跟取值，如果存在多个物种比对，则需要写更为复杂的脚本去取双向最优序列
awk '$3 != 100.000' honeybee_blastp_out.m6 > honeybee_blastp_out_rmself.m6
cut -f1-2 honeybee_blastp_out_rmself.m6 | sort | uniq > honeybee_OR.homolog

conda activate kaks_calculate

ParaAT.pl -h ./honeybee_OR.homolog \
-n all_receptor_gene_cds_uniq.fa \
-a all_receptor_gene_pep_uniq.aa \
-m muscle -p proc -f axt  -g -k -o ./honeybee_kaks
cat *kaks > merge.kaks
awk '!(NR%2)' merge.kaks > honeybee_kaks_calculate.out

# first col OR type (single coexp_RT coexp_NRT)

# in R 
kaks_data<- read.table("honeybee_kaks_calculate.out")
kaks_data$V1<- gsub("-a","_a",kaks_data$V1)
kaks_data$V1<- gsub("-b","_b",kaks_data$V1)
kaks_data$V1<- gsub("-c","_c",kaks_data$V1)
colnames(kaks_data)<- c("OR_pair","method","KA","Ks","dNdS","p_value","length","S-sites","N-sites","Folo-sites","Substitutions","S-sub","N-sub","Folo-S-sub","Folo-N-sub","Divergence-Time","Substitution-Rate-Ratio","GC","ML_score","AICc","Akaike-Weight","Model")
library(tidyr)
df <- kaks_data %>%
  separate(OR_pair, into = c("OR1", "OR2"), sep = "-")
tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")
df$OR1<- gsub("_","-",df$OR1)
df$exp_pattern<- tree_anno[match(df$OR1,tree_anno$OR_gene),]$exp_pattern
plot_data<- df[which(df$exp_pattern%in%c("single_OR","coexp(RT)","coexp(NRT)")),]

my_comparisons <- list( c("coexp(NRT)", "coexp(RT)"),c("coexp(NRT)", "single_OR"),c("coexp(RT)", "single_OR") )

pdf("/md01/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig7/dNdS_compare_split_RTandNRT_honeybee.pdf",width=5,height=5)
ggboxplot(plot_data, x="exp_pattern", y="dNdS", 
  bxp.errorbar = T,width=0.6, notch = F)+ylim(0,1)+
stat_compare_means(comparisons = my_comparisons,method = "t.test")+theme(legend.position="none")+ylab("dN/dS")
dev.off()

# part3: AaegL OR 
# in R 
# the intergenic region length (genomic distance )
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/AaegL5.0/NCBI_annotation_UPDATED_chemoreceptors_CORRECTED.gtf')

Or_gtf<- gtf[grep("Or",gtf$gene_id),]
Or_gtf<- Or_gtf[Or_gtf$type=="CDS",]

Or_list<- unique(Or_gtf$gene_id)

#the distance of OR pair in the same cluster and random OR pair 
coexp_OR_pair<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/17_AaegL_OR/mosquito_OR.csv")
coexp_OR_pair<- coexp_OR_pair[which(coexp_OR_pair$OR1%in%Or_list&coexp_OR_pair$OR2%in%Or_list),]
coexp_OR_pair$length<-NA
for(i in 1:nrow(coexp_OR_pair)){
	OR_pair<- coexp_OR_pair[i,]
	OR1<- OR_pair[,1];
	OR2<- OR_pair[,2];
	OR1_gtf<- as.data.frame(Or_gtf[which(Or_gtf$gene_id%in%c(OR1)),])
    OR2_gtf<- as.data.frame(Or_gtf[which(Or_gtf$gene_id%in%c(OR2)),])
    OR1_seq<- unique(OR1_gtf$seqnames);
    OR2_seq<- unique(OR2_gtf$seqnames);
    if(OR1_seq==OR2_seq){
    	if(unique(OR1_gtf$strand)=="+"){
    		OR1_start<- min(OR1_gtf$start)
    		OR1_end<- max(OR1_gtf$start)
    		OR2_start<- min(OR2_gtf$start)
    		OR2_end<- max(OR2_gtf$start)
    		tmp<- c(OR1_start,OR1_end,OR2_start,OR2_end)
    		tmp<- tmp[order(tmp)]
    		coexp_OR_pair$length[i]<- tmp[3]-tmp[2]
    	}else{
    		OR1_start<- max(OR1_gtf$start)
    		OR1_end<- min(OR1_gtf$start)
    		OR2_start<- max(OR2_gtf$start)
    		OR2_end<- min(OR2_gtf$start)
    		tmp<- c(OR1_start,OR1_end,OR2_start,OR2_end)
    		tmp<- tmp[order(tmp)]
    		coexp_OR_pair$length[i]<- tmp[3]-tmp[2]
    	}
    }
}



probability<-data.frame(OR_pair_type=c("not same chromosome","same chromosome"),
    count=c(24,44));

pdf("OR_pair_appearing_on_same_chromosome.pdf",width=3,height=3)
p<-ggplot(data = probability, aes_string(x = "OR_pair_type", y = "count",color="OR_pair_type")) +  
        xlab("OR pair type") +
        ylab("The count on the same chromosome") + 
        geom_bar(stat = "identity", width = 0.6,fill="white") +
        theme_classic()+guides(color='none')+
        theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
#add number in plot 
a<- p+geom_text(aes(label = count), size = 3, hjust = 0.5, vjust = 3)
a
p
dev.off()

length<- as.numeric(na.omit(coexp_OR_pair$length))
length<- length[order(length)]

# Fig2K: OR transcriptomic distance 
# the transcript distance between coexp and uncoexp
library(pheatmap)
DefaultAssay(ORN)<-"SCT"
matrix<- ORN@assays$SCT[dotplot_feature,]
cor_data<-cor(as.data.frame(t(matrix)))
cosine_dist <- (1-cosine(cor_data));





