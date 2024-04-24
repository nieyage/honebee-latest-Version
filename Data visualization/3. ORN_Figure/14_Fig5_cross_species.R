# fly, mosquito and honeybee statistics

# functional coexp vs single OR;
coexp_OR<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/coexp_OR.csv")
# Part1:  
# Number of OR gene expressed per cell 
# get the count matrix and OR list 
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
set.seed(1234)
# honeybee 
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
honeybee_OR<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
ORN <- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_latest.rds")
ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")
honeybee_OR<- ORgene_name_trans[match(honeybee_OR,ORgene_name_trans$OR_gene),]$last_name
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
mosquito_neuron<-readRDS("/data/R02/nieyg/project/honeybee/data/publish_data/mosquito/SeuratObject2_Antenna_mergedBatches_Neurons.rds")
Fly_ORN<-readRDS("/md01/liyh526/project/Fly-cooperate/7.18run/WNN_ORN_integrated_antenna.rds")
# get the OR list 
mosquito_OR<- grep("^Or[1-9].*",rownames(mosquito_neuron),value=TRUE)
fly_OR<- grep("^Or[1-9].*",rownames(Fly_ORN),value=TRUE)
honeybee_OR<- honeybee_OR[-which(honeybee_OR%in%Orco)]

# parpare the OR pair
coexp_OR$type<- "coexp";
honeybee_coexp_OR<- unique(c(coexp_OR[coexp_OR$species=="honeybee",]$OR1,coexp_OR[coexp_OR$species=="honeybee",]$OR2))
honeybee_single_OR<- honeybee_OR[-which(honeybee_OR%in%honeybee_coexp_OR)]
fly_coexp_OR<- unique(c(coexp_OR[coexp_OR$species=="fly",]$OR1,coexp_OR[coexp_OR$species=="fly",]$OR2))
fly_single_OR<- fly_OR[-which(fly_OR%in%fly_coexp_OR)]
mosquito_coexp_OR<- unique(c(coexp_OR[coexp_OR$species=="mosquito",]$OR1,coexp_OR[coexp_OR$species=="mosquito",]$OR2))
mosquito_single_OR<- mosquito_OR[-which(mosquito_OR%in%mosquito_coexp_OR)]
mosquito_single_OR<- mosquito_single_OR[which(mosquito_single_OR %in% unique(AaegL_gtf$gene_id))]
single_OR<- data.frame()

# in honeybee
remaining_gene<- honeybee_single_OR
for(i in 1:length(honeybee_single_OR)){
    OR1<- honeybee_single_OR[i];
    remaining_gene<- remaining_gene[-which(remaining_gene%in%OR1)]
    for(j in 1:length(remaining_gene)){
        OR2<- remaining_gene[j];
        tmp<- data.frame(OR1=OR1,OR2=OR2,species="honeybee",NB=1,type="non-coexp");
        single_OR<- rbind(single_OR,tmp)
    }
}

# in fly
remaining_gene<- fly_single_OR
for(i in 1:length(fly_single_OR)){
    OR1<- fly_single_OR[i];
    remaining_gene<- remaining_gene[-which(remaining_gene%in%OR1)]
    for(j in 1:length(remaining_gene)){
        OR2<- remaining_gene[j];
        tmp<- data.frame(OR1=OR1,OR2=OR2,species="fly",NB=1,type="non-coexp");
        single_OR<- rbind(single_OR,tmp)
    }
}

# in mosquito
remaining_gene<- mosquito_single_OR
for(i in 1:length(mosquito_single_OR)){
    OR1<- mosquito_single_OR[i];
    remaining_gene<- remaining_gene[-which(remaining_gene%in%OR1)]
    for(j in 1:length(remaining_gene)){
        OR2<- remaining_gene[j];
        tmp<- data.frame(OR1=OR1,OR2=OR2,species="mosquito",NB=1,type="non-coexp");
        single_OR<- rbind(single_OR,tmp)
    }
}

single_OR<- na.omit(single_OR)
all_OR_pair<- rbind(single_OR,coexp_OR)

# parpare the OR region file (from gtf)
chemoreceptor$gene_name <- ORgene_name_trans[match(chemoreceptor$gene_name,ORgene_name_trans$OR_gene),]$last_name
AaegL_gtf <- rtracklayer::import('/md01/nieyg/ref/10X/AaegL5.0/NCBI_annotation_UPDATED_chemoreceptors_CORRECTED.gtf')
AaegL_gtf<- AaegL_gtf[grep("Or",AaegL_gtf$gene_id),]
AaegL_gtf<- AaegL_gtf[AaegL_gtf$type=="CDS",]
AaegL_gtf<- as.data.frame(AaegL_gtf)
fly_gtf <- rtracklayer::import('/md01/nieyg/ref/10X/fly/Drosophila_melanogaster/genes/genes.gtf.gz')
fly_gtf<- fly_gtf[grep("^Or",fly_gtf$gene_name),]
fly_gtf<- fly_gtf[fly_gtf$type=="CDS",]
fly_gtf<- as.data.frame(fly_gtf)

all_OR_pair$same_Chrom<- NA
all_OR_pair$genomic_dist<- NA
for (i in 1:nrow(all_OR_pair)){
    tmp<- all_OR_pair[i,]
    if(tmp$species=="honeybee"){
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        OR1_seqnames<- chemoreceptor[chemoreceptor$gene_name==OR1,1][1];
        OR2_seqnames<- chemoreceptor[chemoreceptor$gene_name==OR2,1][1];
        if(OR1_seqnames==OR2_seqnames){
            all_OR_pair[i,]$same_Chrom="Yes";
            if(unique(chemoreceptor[chemoreceptor$gene_name==OR1,]$strand)=="+"){
                OR1_start<- chemoreceptor[chemoreceptor$gene_name==OR1,]$start[1];
                OR2_start<- chemoreceptor[chemoreceptor$gene_name==OR2,]$start[1];
                all_OR_pair[i,]$genomic_dist=abs(OR1_start-OR2_start)
            }else{
                OR1_start<- chemoreceptor[chemoreceptor$gene_name==OR1,]$end[1];
                OR2_start<- chemoreceptor[chemoreceptor$gene_name==OR2,]$end[1];
                all_OR_pair[i,]$genomic_dist=abs(OR1_start-OR2_start)
            }
        }else{all_OR_pair[i,]$same_Chrom="No"};

    }
    if(tmp$species=="fly"){
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        if(OR1%in%unique(fly_gtf$gene_name) &OR2%in%unique(fly_gtf$gene_name) ){
        OR1_seqnames<- fly_gtf[fly_gtf$gene_name==OR1,1][1];
        OR2_seqnames<- fly_gtf[fly_gtf$gene_name==OR2,1][1];
        if(OR1_seqnames==OR2_seqnames){
            all_OR_pair[i,]$same_Chrom="Yes";
            if(unique(fly_gtf[fly_gtf$gene_name==OR1,]$strand)=="+"){
                OR1_start<- fly_gtf[fly_gtf$gene_name==OR1,]$start[1];
                OR2_start<- fly_gtf[fly_gtf$gene_name==OR2,]$start[1];
                all_OR_pair[i,]$genomic_dist=abs(OR1_start-OR2_start)
            }else{
                OR1_start<- fly_gtf[fly_gtf$gene_name==OR1,]$end[1];
                OR2_start<- fly_gtf[fly_gtf$gene_name==OR2,]$end[1];
                all_OR_pair[i,]$genomic_dist=abs(OR1_start-OR2_start)
            }
        }else{all_OR_pair[i,]$same_Chrom="No"};}
    }
    if(tmp$species=="mosquito"){
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        if(OR1%in%unique(AaegL_gtf$gene_id) &OR2%in%unique(AaegL_gtf$gene_id) ){
        OR1_seqnames<- AaegL_gtf[AaegL_gtf$gene_id==OR1,1][1];
        OR2_seqnames<- AaegL_gtf[AaegL_gtf$gene_id==OR2,1][1];
        if(OR1_seqnames==OR2_seqnames){
            all_OR_pair[i,]$same_Chrom="Yes";
            if(unique(AaegL_gtf[AaegL_gtf$gene_id==OR1,]$strand)=="+"){
                OR1_start<- AaegL_gtf[AaegL_gtf$gene_id==OR1,]$start[1];
                OR2_start<- AaegL_gtf[AaegL_gtf$gene_id==OR2,]$start[1];
                all_OR_pair[i,]$genomic_dist=abs(OR1_start-OR2_start)
            }else{
                OR1_start<- AaegL_gtf[AaegL_gtf$gene_id==OR1,]$end[1];
                OR2_start<- AaegL_gtf[AaegL_gtf$gene_id==OR2,]$end[1];
                all_OR_pair[i,]$genomic_dist=abs(OR1_start-OR2_start)
            }
        }else{all_OR_pair[i,]$same_Chrom="No"};
    }}
}



all_OR_pair<- all_OR_pair[-which(is.na(all_OR_pair$same_Chrom)),]

table(all_OR_pair$species,all_OR_pair$type);
table(all_OR_pair$same_Chrom,all_OR_pair$type,all_OR_pair$species);

chrom_data<- data.frame(species=c("honeybee","honeybee","fly","fly","mosquito","mosquito"),
    type=rep(c("non-coexp","coexp"),3),
    total=c(9870,4,1326,5,1378,65),
    same=c(2655,4,281,5,466,41)
    )

chrom_data$ratio<- chrom_data$same/chrom_data$total
chrom_data$type<- factor(chrom_data$type,levels=c("non-coexp","coexp"))
# barplot by group 
#ggplot2 分组柱状图
p <- ggplot(chrom_data, aes(species, ratio, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6)  +
  scale_fill_manual(values = c("#3C5488B2", "#F39B7FB2")) +
  labs(title = NULL, x = NULL, y = 'The probability on the same chromosome', fill = NULL) +
  theme_classic()+guides(color='none')+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
pdf("./00_Figure/Fig5/Fig5C-OR_pair_appearing_on_same_chromosome.pdf",width=5,height=3)
p;
dev.off()


# Genomic distance 
Genomic_data<- na.omit(all_OR_pair)

bin.levels = c(as.character(1:10));

honeybee_data<- Genomic_data[Genomic_data$species=="honeybee",]
honeybee_data$genomic_dist[honeybee_data$genomic_dist>100000]<- 100000
honeybee.comp.bins <- as.data.frame(honeybee_data) %>%
  mutate(bin.genome.10th = as.numeric(cut(genomic_dist, breaks = 10)),
         #bin.genome.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
         bin.genome.10th = factor(bin.genome.10th, levels = bin.levels),
         coexp=ifelse(type == "coexp", 1, 0),noncoexp=ifelse(type == "non-coexp", 1, 0)) 
honeybee.comp.stat <- honeybee.comp.bins %>% group_by(species,bin.genome.10th,type) %>%
       summarise(coexp_count = sum(coexp),noncoexp_count = sum(noncoexp))

p <- ggplot(honeybee.comp.bins, aes(bin.genome.10th, type, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6)  +
  scale_fill_manual(values = c("#0576A4", "#EE7A07","#D02D36")) +
  labs(title = NULL, x = 'genome distance (bin #)', y = 'Number of coexp_OR', fill = NULL) +
  theme_classic()+guides(color='none')+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
pdf("./00_Figure/Fig5/Fig5C-coexOR_pair_appearing_on_genomic_bin.pdf",width=6,height=3)
p;
dev.off()



fly_data<- Genomic_data[Genomic_data$species=="fly",]
fly_data$genomic_dist[fly_data$genomic_dist>2000000]<- 2000000
fly.comp.bins <- as.data.frame(fly_data) %>%
  mutate(bin.genome.10th = as.numeric(cut(genomic_dist, breaks = 10)),
         #bin.genome.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
         bin.genome.10th = factor(bin.genome.10th, levels = bin.levels),
         coexp=ifelse(type == "coexp", 1, 0),noncoexp=ifelse(type == "non-coexp", 1, 0)) 


mosquito_data<- Genomic_data[Genomic_data$species=="mosquito",]
mosquito_data$genomic_dist[mosquito_data$genomic_dist>80000000]<- 80000000
mosquito.comp.bins <- as.data.frame(mosquito_data) %>%
  mutate(bin.genome.10th = as.numeric(cut(genomic_dist, breaks = 10)),
         #bin.genome.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
         bin.genome.10th = factor(bin.genome.10th, levels = bin.levels),
         coexp=ifelse(type == "coexp", 1, 0),noncoexp=ifelse(type == "non-coexp", 1, 0)) 

comp.bins<- rbind(honeybee.comp.bins,fly.comp.bins,mosquito.comp.bins)


pairwise_distances.comp.stat <- comp.bins %>% group_by(species,bin.genome.10th) %>%
       summarise(coexp_count = sum(coexp),noncoexp_count = sum(noncoexp))
pairwise_distances.comp.stat$ratio<- pairwise_distances.comp.stat$coexp_count/(pairwise_distances.comp.stat$coexp_count+pairwise_distances.comp.stat$noncoexp_count)


p <- ggplot(pairwise_distances.comp.stat, aes(bin.genome.10th, coexp_count, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6)  +
  scale_fill_manual(values = c("#0576A4", "#EE7A07","#D02D36")) +
  labs(title = NULL, x = 'genome distance (bin #)', y = 'Number of coexp_OR', fill = NULL) +
  theme_classic()+guides(color='none')+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
pdf("./00_Figure/Fig5/Fig5C-coexOR_pair_appearing_on_genomic_bin.pdf",width=6,height=3)
p;
dev.off()

# parpare the sequence data

# Sequence similarity
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
honeybee_aa <-c(OR_fasta,supply_fasta)
names(honeybee_aa) <- ORgene_name_trans[match(names(honeybee_aa),ORgene_name_trans$OR_gene),]$last_name

AaegL5_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/AaegL5-ChemoreceptorPeptides.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
names_in_myXStringSet <- names(AaegL5_fasta)
#select names from your XStringSet
Aaeg_OR <- AaegL5_fasta[names_in_myXStringSet[grep("AaegOr",names_in_myXStringSet)],]
names(Aaeg_OR) <- gsub("Aaeg","",names(Aaeg_OR) )
names(Aaeg_OR) <- gsub(" new","",names(Aaeg_OR) )

#get fly OR sequence;
library(biomaRt)
mart <- useMart("ensembl","dmelanogaster_gene_ensembl")
OR_ppseqs <- getSequence(id = fly_OR,
                          type="external_gene_name",
                          seqType="peptide",
                          mart = mart)                               
OR_ppseqs <- OR_ppseqs[OR_ppseqs$peptide!="Sequence unavaliable",]
OR_ppseqs$external_gene_name <- make.unique(OR_ppseqs$external_gene_name)
exportFASTA(OR_ppseqs,"./fly_OR_pep.fasta")
fly_aa<-readAAStringSet("./fly_OR_pep.fasta", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)

dot_aa<- honeybee_aa[which(names(honeybee_aa)%in% dotplot_feature),]

aln <- muscle::muscle(honeybee_aa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
dist_matrix <- as.matrix(sdist)
honeybee_dist<-dist_matrix
tree<-as.phylo(clust)
pdf("/data/R02/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotOR_sequence_protein_similarity-tree_add_groupinfo.pdf",width=16,height=16)
ggtree(tree, layout="fan", ladderize = FALSE, branch.length = "none") + geom_tiplab2(size=3) + theme(legend.position = "right")
dev.off()


aln <- muscle::muscle(Aaeg_OR)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
dist_matrix <- as.matrix(sdist)
mosquito_dist<-dist_matrix
aln <- muscle::muscle(fly_aa)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming") 
clust <- hclust(sdist)
dist_matrix <- as.matrix(sdist)
fly_dist<-dist_matrix

coexp_OR_pair<- Genomic_data
coexp_OR_pair$Sequence_dist<-NA
for (i in 1:nrow(coexp_OR_pair)){
    tmp<- coexp_OR_pair[i,]
    if(tmp$species=="honeybee"){
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        if(OR1%in% rownames(honeybee_dist) & OR2%in% colnames(honeybee_dist)){
            coexp_OR_pair$Sequence_dist[i]<- honeybee_dist[OR1,OR2]}
    }
    if(tmp$species=="fly"){
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        if(OR1%in% rownames(fly_dist) & OR2%in% colnames(fly_dist)){
            coexp_OR_pair$Sequence_dist[i]<- fly_dist[OR1,OR2]}
    }
    if(tmp$species=="mosquito"){
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        if(OR1%in% rownames(mosquito_dist) & OR2%in% colnames(mosquito_dist)){
            coexp_OR_pair$Sequence_dist[i]<- mosquito_dist[OR1,OR2]}
    }
}

# only keep the NB pair 
OR_pair_NB<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/ORpairtype.csv")
OR_pair_NB$Sequence_dist<-NA
for (i in 1:nrow(OR_pair_NB)){
    tmp<- OR_pair_NB[i,]
        OR1<- tmp$OR1;
        OR2<- tmp$OR2;
        if(OR1%in% rownames(honeybee_dist) & OR2%in% colnames(honeybee_dist)){
            OR_pair_NB$Sequence_dist[i]<- honeybee_dist[OR1,OR2]}
}
OR_pair_NB<- na.omit(OR_pair_NB)

bin.levels = c(as.character(1:6),"7-10");
honeybee.comp.bins <- as.data.frame(OR_pair_NB) %>%
   mutate(bin.seq.10th = as.numeric(cut(Sequence_dist, breaks = 10)),
          bin.seq.10th = ifelse(bin.seq.10th >= 7, "7-10", as.character(bin.seq.10th)),
          bin.seq.10th = factor(bin.seq.10th, levels = bin.levels)) 

data<- as.data.frame(table(honeybee.comp.bins$bin.seq.10th,honeybee.comp.bins$exp_type))
data$Var2<- factor(data$Var2,levels=c("single","MP","RT"))
p <- ggplot(data, aes(Var1, Freq, fill = Var2)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6)  +
  scale_fill_manual(values = c("#129FAF", "#DE7C5B","#FBD277")) +
  labs(title = NULL, x = 'sequence distance (bin #)', y = 'Number of OR pair', fill = NULL) +
  theme_classic()+guides(color='none')+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
pdf("./00_Figure/Fig5/Fig5E-NB_coexOR_pair_appearing_on_sequence_distance_bin.pdf",width=6,height=3)
p;
dev.off()

#coexp_OR_pair<- coexp_OR_pair[coexp_OR_pair$type=="coexp",]
bin.levels = c(as.character(1:10));

 honeybee_data<- coexp_OR_pair[coexp_OR_pair$species=="honeybee",]
 honeybee_data$Sequence_dist[honeybee_data$Sequence_dist>300]<- 600
 honeybee.comp.bins <- as.data.frame(honeybee_data) %>%
   mutate(bin.seq.10th = as.numeric(cut(Sequence_dist, breaks = 10)),
          #bin.seq.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
          bin.seq.10th = factor(bin.seq.10th, levels = bin.levels),
          coexp=ifelse(type == "coexp", 1, 0),noncoexp=ifelse(type == "non-coexp", 1, 0)) 
 

 fly_data<- coexp_OR_pair[coexp_OR_pair$species=="fly",]
 fly_data$Sequence_dist[fly_data$Sequence_dist>300]<- 600
 fly.comp.bins <- as.data.frame(fly_data) %>%
   mutate(bin.seq.10th = as.numeric(cut(Sequence_dist, breaks = 10)),
          #bin.genome.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
          bin.seq.10th = factor(bin.seq.10th, levels = bin.levels),
          coexp=ifelse(type == "coexp", 1, 0),noncoexp=ifelse(type == "non-coexp", 1, 0)) 
 
 mosquito_data<- coexp_OR_pair[coexp_OR_pair$species=="mosquito",]
 #mosquito_data$Sequence_dist[mosquito_data$Sequence_dist>300]<- 600
 mosquito.comp.bins <- as.data.frame(mosquito_data) %>%
   mutate(bin.seq.10th = as.numeric(cut(Sequence_dist, breaks = 10)),
          #bin.genome.10th = ifelse(bin.genome.10th >= 7, "7-10", as.character(bin.genome.10th)),
          bin.seq.10th = factor(bin.seq.10th, levels = bin.levels),
          coexp=ifelse(type == "coexp", 1, 0),noncoexp=ifelse(type == "non-coexp", 1, 0)) 

comp.bins<- rbind(honeybee.comp.bins,fly.comp.bins,mosquito.comp.bins)
comp.bins<- na.omit(comp.bins)

pairwise_distances.comp.stat <- comp.bins %>% group_by(species,bin.seq.10th,.drop = FALSE) %>%
       summarise(coexp_count = sum(coexp),noncoexp_count = sum(noncoexp));
pairwise_distances.comp.stat$ratio<- pairwise_distances.comp.stat$coexp_count/(pairwise_distances.comp.stat$coexp_count+pairwise_distances.comp.stat$noncoexp_count)*100

p <- ggplot(pairwise_distances.comp.stat, aes(bin.seq.10th, coexp_count, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6)  +
  scale_fill_manual(values = c("#0576A4", "#EE7A07","#D02D36")) +
  labs(title = NULL, x = 'sequence distance (bin #)', y = 'Number of coexp_OR', fill = NULL) +
  theme_classic()+guides(color='none')+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
pdf("./00_Figure/Fig5/Fig5D-coexOR_pair_appearing_on_sequence_distance_bin.pdf",width=6,height=3)
p;
dev.off()

p <- ggplot(pairwise_distances.comp.stat, aes(bin.seq.10th, ratio, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6)  +
  scale_fill_manual(values = c("#0576A4", "#EE7A07","#D02D36")) +
  labs(title = NULL, x = 'sequence distance (bin #)', y = 'coexp_OR/total', fill = NULL) +
  theme_classic()+guides(color='none')+
  theme(axis.text.x = element_text(angle = 0,vjust = 0.5,hjust = 0.5))
pdf("./00_Figure/Fig5/Fig5D-ratio_coexOR_pair_appearing_on_sequence_distance_bin.pdf",width=6,height=3)
p;
dev.off()

3 
5
cross_species_glomeruli<-data.frame(species=c("Apis mellifera","D. melanogaster"),
  MP=c(5,3));
cross_species_glomeruli$species<- factor(cross_species_glomeruli$species,levels=c("D. melanogaster","Apis mellifera"))
pdf("./00_Figure/Fig5/Fig5C-cross_species_MP.pdf",width=3,height=4)
p<-ggplot(data = cross_species_glomeruli, aes_string(x = "species", y = "MP", 
        fill = "species")) +  xlab(" ") + ylab("# of coexp_OR group") + 
        geom_bar( stat = "identity",width=0.6,color = 'black', fill='grey') +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = MP), size = 3, hjust = 0.5, vjust = 3) 
dev.off();




