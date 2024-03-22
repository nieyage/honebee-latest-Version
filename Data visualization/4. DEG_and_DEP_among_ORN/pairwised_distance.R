# cluster DEG and DEP 
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
require(tidyverse)
require(ggpubr)
require(patchwork)
require(lsa)
require(ggtree)
require(ggtreeExtra)
require(ggnewscale)
require(viridis)
require(fields)
require(ggrastr)
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
DefaultAssay(ORN)<- "raw_RNA"
OR_pair<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/19_without_OR/OR_pair_genomic_dist.csv")
OR_data<- ORN@assays$RNA
OR_data<-OR_data[which(rownames(OR_data)%in%all_receptor_gene[-which(all_receptor_gene%in%Orco)]),]
OR_data<-as.matrix(OR_data);
cell <- names(which(colSums(OR_data)!=0))
OR_data<- OR_data[,cell]
row_max <- apply(OR_data, 2, which.max)
# 根据行号获取行名
row_names <- rownames(OR_data)[row_max]
ORN<- subset(ORN,cells=cell)
names(row_names)<- cell
row_names<- row_names[colnames(ORN)]
ORN$OR_max<- row_names;


library(Matrix)
library(patchwork)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(effsize)
# Compute the euclidean distance 
pca.coord.without.or <- ORN@reductions$pca_withOR@cell.embeddings[,1:20]
dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)
colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)
osn.idents <- ORN$OR_max

dist.mat.without.or[lower.tri(dist.mat.without.or, diag = TRUE)] <- NA
 
dist.mat.without.or <- melt(data = dist.mat.without.or,
                            na.rm = TRUE,
                            varnames = c("OSN_1", "OSN_2"),
                            value.name = "euclidean_distance")

dist.mat.without.or <- data.table::data.table(dist.mat.without.or)
dist.mat.without.or$OSN_1 <- as.character(x = dist.mat.without.or$OSN_1)
dist.mat.without.or$OSN_2 <- as.character(x = dist.mat.without.or$OSN_2)

dist.mat.without.or$OSN_1_ident <- osn.idents[dist.mat.without.or$OSN_1]
dist.mat.without.or$OSN_2_ident <- osn.idents[dist.mat.without.or$OSN_2]

dist.mat.without.or$osn_pairs <- ifelse(test = dist.mat.without.or$OSN_1_ident == dist.mat.without.or$OSN_2_ident,
                                        yes = "same OSN pop.",
                                        no = "different OSN pop.")

# OR phylogeny
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
all_OR_gene_fasta<- c(OR_fasta,supply_fasta)
aln <- muscle::muscle(all_OR_gene_fasta)
auto <- maskGaps(aln, min.fraction=0.5, min.block.width=4)
sdist <- stringDist(as(auto,"AAStringSet"), method="hamming")
clust <- hclust(sdist,method="complete")#"ward.D"’, ‘"ward.D2"’,‘"single"’, ‘"complete"’, ‘"average"’ (= UPGMA), ‘"mcquitty"’, ‘"median"’ or ‘"centroid"’ (= UPGMC)
tree <- as.phylo(clust)
dist_matrix <- as.matrix(sdist)
dist_percent<-(max(dist_matrix)-dist_matrix)/max(dist_matrix)*100
OR_pair$Sequence_similarity<- NA
for (i in 1:nrow(OR_pair)){
    gene1<-OR_pair[i,]$OR1;
    gene2<-OR_pair[i,]$OR2;
    if(gene1 %in% names(all_OR_gene_fasta)){
    	if(gene2 %in% names(all_OR_gene_fasta)){
    		    OR_pair$Sequence_similarity[i]<-dist_percent[gene1,gene2];
    	}
    }
}


PCA  <- ORN@reductions$pca_withOR@cell.embeddings
osn_ident <- ORN$OR_max

# fetch values of a given set of PCs
  osn_pca <- PCA[,1:20] %>%  
    as.data.frame(.) %>%
    rownames_to_column("cell.ident") %>%
    mutate(olfr.ident = osn_ident[cell.ident])
  
  # split by OSN population
  osn_pca.list <- split(osn_pca, osn_pca$olfr.ident)
  

  # calculate centroids for OSN populations
  osn_pca.centroids.matrix <- lapply(osn_pca.list, function(mat) {
    mat <- mat[,which(startsWith(colnames(mat), "pca_withor_"))]
    
    # optionally: outlier identification using mean distance of k-nearest
    # neighbors
    find_centroid(mat) }) %>%
    do.call(what = rbind, args = .)
  
  # calculate distances between centroids and make the pairwise_distances dataset
  #osn_dist <- rdist(as.matrix(osn_pca.centroids.matrix))
  osn_dist <- 1-cosine(t(as.matrix(osn_pca.centroids.matrix)))

  rownames(osn_dist) <- rownames(osn_pca.centroids.matrix)
  colnames(osn_dist) <- rownames(osn_pca.centroids.matrix)

  pairwise_distances <- as.data.frame(osn_dist) %>%
    rownames_to_column("gene_name.1") %>%
    pivot_longer(cols = where(is.numeric), values_to = "dist.trans", names_to = "gene_name.2") %>%
    mutate(pair = unlist(mapply(
      a = gene_name.1,
      b = gene_name.2,
      SIMPLIFY = F,
      FUN = function(a,b) {
        x <- sort(c(a,b))
        return(paste(x[1], x[2], sep = ":"))
      }
    ))) 

summary(na.omit(pairwise_distances$dist.trans))
# add metadata and genomic distance
pairwise_distances <- left_join(pairwise_distances, 
                                select(OR_pair,OR1, OR2, genomic_dist, Sequence_similarity), 
                                by = c("gene_name.1" = "OR1"));
pairwise_distances <- pairwise_distances[pairwise_distances$OR2 == pairwise_distances$gene_name.2, ]

# unique comparisons
pairwise_distances.comp <- arrange(pairwise_distances, gene_name.1, gene_name.2) %>%
  # remove self-comparisons
  filter(gene_name.1 != gene_name.2) %>%
  distinct()

# B) Distribution of pairwise amino acid difference across bins of 
#    transcriptomic distances
pairwise_distances.comp<- na.omit(pairwise_distances.comp)
# Almost identical sequence definition 
dist.aa.limits <- filter(pairwise_distances.comp) %>%
  pivot_longer(cols = "Sequence_similarity", 
               names_to = "seq.range", 
               values_to = "dist.aa") %>%
  group_by(seq.range) %>%
  summarise(almost_ident.thresh = quantile(dist.aa, probs = 0.60))


dist.aa.limits.names <- dist.aa.limits$seq.range
dist.aa.limits <- pull(dist.aa.limits, almost_ident.thresh)
names(dist.aa.limits) <- dist.aa.limits.names

dist.aa.limits["Sequence_similarity"] <- 70

# 95th percentile of the distribution of intergenic distance
intergenic.dist<-pairwise_distances.comp$genomic_dist
intergenic.dist.mean <- mean(intergenic.dist)
intergenic.dist.90th <- quantile(intergenic.dist, probs = 0.9)
dist.genome.limits <- intergenic.dist.90th
names(dist.genome.limits) <- "dist.genome"

dist.limits <- c(dist.aa.limits, dist.genome.limits)


# definintion of bins of transcriptomic distances values for all pairwise
# distances within clusters
bin.levels = c(as.character(1:6), "7-10")

pairwise_distances.comp.bins <- pairwise_distances.comp %>%
  mutate(bin.trans.10th = as.numeric(cut(dist.trans, breaks = 10)),
         bin.trans.10th = ifelse(bin.trans.10th >= 7, "7-10", as.character(bin.trans.10th)),
         bin.trans.10th = factor(bin.trans.10th, levels = bin.levels),
         bin.genome.10th = as.numeric(genomic_dist),
         bin.genome.10th = floor(bin.genome.10th/1000),
         bin.genome.10th = ifelse(bin.genome.10th >= 10, ">10", as.character(bin.genome.10th)),
         bin.genome.4genes = ceiling(log2(genomic_dist)-log2(intergenic.dist.mean*4)),
         bin.genome.4genes = as.character(ifelse(bin.genome.4genes < 0, 0, bin.genome.4genes) + 1),
         bin.genome.4genes = factor(bin.genome.4genes, levels = as.character(1:7)),
         bin.aa.10th = as.numeric(cut(Sequence_similarity, breaks = 10)),
         bin.aa.10th = as.character(bin.aa.10th),
         bin.aa.10th = factor(bin.aa.10th, levels = as.character(1:10))) %>%
  select(pair, genomic_dist, dist.trans, 
         Sequence_similarity, bin.trans.10th, bin.genome.10th, bin.genome.4genes, bin.aa.10th) 

# statistics on pairwise_distances.comp.bins combuted with transcriptomic 
# distance bins
pairwise_distances.comp.transbins.stat <- group_by(pairwise_distances.comp.bins, bin.trans.10th) %>%
  summarise(n = n(),
            almost_ident.prop = length(which(Sequence_similarity > 70))/n,
            nearby.prop = length(which(genomic_dist < 3000))/n,
            cor.aa.pearson = cor(dist.trans,Sequence_similarity, method = "pearson"),
            cor.aa.spearman = cor(dist.trans, Sequence_similarity, method = "spearman"),
            cor.genome.pearson = cor(dist.trans, genomic_dist, method = "pearson"),
            cor.genome.spearman = cor(dist.trans, genomic_dist, method = "spearman"))

pairwise_distances.comp.noconfound <- bind_rows(
  # distant & homologous
  filter(pairwise_distances.comp, 
       Sequence_similarity < 70) %>%
    mutate(cat = "D-H",
           genome = "distant",
           aa = "homologous"),
  
  # close & non-homologous
  filter(pairwise_distances.comp, 
         Sequence_similarity < 70,
         genomic_dist < 3000) %>%
    mutate(cat = "C-N",
           genome = "close",
           aa = "non-homologous"),
  
  # distant & non-homologous
  filter(pairwise_distances.comp, 
         Sequence_similarity > 70) %>%
    mutate(cat = "D-N",
           genome = "distant",
           aa = "non-homologous"),
  
  # close & homologous
  filter(pairwise_distances.comp, 
         Sequence_similarity > 70,
         genomic_dist < 3000) %>%
    mutate(cat = "C-H",
           genome = "close",
           aa = "homologous")
) %>%
  mutate(cat = factor(cat, levels = c("C-H", "C-N", "D-H", "D-N")))

pairwise_distances.comp.noconfound.stat <- select(pairwise_distances.comp.noconfound, cat) %>%
  group_by(cat) %>%
  summarise(n = n())

# common layout
common_layout <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.7, "lines")
  )

p.close_distant <- ggplot(pairwise_distances.comp.noconfound, aes(x=cat, y=dist.trans)) +
  geom_violin(color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  stat_compare_means(aes(label = ..p.signif..), 
                     method = "wilcox.test", comparisons = list(
                       c("C-H", "C-N"),
                       c("D-H", "D-N"), 
                       c("C-H", "D-H")))+
                     #label.y = c(28, 32, 34)) +
  common_layout +  ylab("transcriptomic distance")
  #scale_y_continuous(limits = c(0, 38), expand = c(0,0), breaks = c(10, 20, 30)) +


# 
p.binning.prop <- ggplot(pairwise_distances.comp.transbins.stat, 
                         aes(x=bin.trans.10th)) +
  geom_col(aes(y=almost_ident.prop), fill = "#8D85BE", color = "#8D85BE", width = 0.4, position = position_nudge(-0.25)) +
  stat_cor(aes(x=as.numeric(bin.trans.10th), y=almost_ident.prop), method = "spearman",label.y = 1, label.x = 1,   color = "#8D85BE") +
  geom_col(aes(y=nearby.prop), fill = "#66C2A4", color = "#66C2A4", width = 0.4, position = position_nudge(+0.25)) +
  stat_cor(aes(x=as.numeric(bin.trans.10th), y=nearby.prop), method = "spearman",label.y = 1, label.x = 4,  color = "#66C2A4") +
  scale_x_discrete(limits = c(as.character(1:6), "7-10")) +
  #scale_y_continuous(limits = c(0, 0.33), expand = c(0,0)) +
  common_layout +
  xlab("transcriptomic distance bin") +
  ylab("proportion of pairs")

dist.trans.cuts <- levels(cut(pairwise_distances.comp.bins$dist.trans, breaks = 10)) %>%
  str_remove(., ".+,") %>%
  str_remove(., "\\]") %>%
  as.numeric()

dist.aa.cuts <- levels(cut(pairwise_distances.comp.bins$Sequence_similarity, breaks = 10)) %>%
  str_remove(., ".+,") %>%
  str_remove(., "\\]") %>%
  as.numeric()

dist.genome.cuts <- levels(cut(pairwise_distances.comp.bins$genomic_dist, breaks = 10)) %>%
  str_remove(., ".+,") %>%
  str_remove(., "\\]") %>%
  as.numeric()

dist.genome.cuts.4genes <- intergenic.dist.mean*4*2**c(0, seq_along(unique(as.numeric(pairwise_distances.comp.bins$bin.genome.4genes))))

pairwise_distances.comp.bins$bin.genome.10th<- factor(pairwise_distances.comp.bins$bin.genome.10th,levels=c("1","2","3","4","5","6","7","8","9",">10"))
# binning and violin plots of transcriptomic distance vs genomic distance
p.binning.genome <- ggplot(pairwise_distances.comp.bins, 
                       aes(x=bin.genome.10th, y=dist.trans)) +
  geom_violin(aes(group=bin.genome.10th), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  #scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance (*1000)") +
  ylab("transcriptomic distance")

# scatterplot of the same dataset
p.scatterplot.genome <- ggplot(pairwise_distances.comp.bins, 
                               aes(x=dist.genome, y=dist.trans)) +
  rasterize(geom_point(alpha=0.3, size = 0.2, color="dodgerblue4", shape = 19), dpi = 600) +
  stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 10, alpha = 0.4, show.legend = F) +
  #geom_vline(xintercept = dist.genome.cuts, linetype = "dashed") +
  #geom_vline(xintercept = dist.genome.cuts.4genes[1:6], linetype = "dashed", color = "red") +
  #stat_smooth(method = "lm", se = F) +
  #stat_cor(method = "spearman") +
  #stat_smooth(data = filter(pairwise_distances.comp.bins, as.numeric(bin.genome.10th) <= 3),
  #            method = "lm", se = F, color = "red") +
  #scale_x_continuous(limits = c(0, 5.1e+6), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance(*1000)") +
  ylab("transcriptomic distance")


p.binning.genome.2 <- ggplot(pairwise_distances.comp.bins, 
                           aes(x=bin.genome.4genes, y=dist.trans)) +
  geom_violin(aes(group=bin.genome.4genes), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
 # scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("genomic distance (bin #)") +
  ylab("transcriptomic distance")

x <- filter(pairwise_distances.comp.bins, bin.genome.10th %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.genome.10th %in% c("1", "2", "3")) %>%
  pull(genomic_dist)
p.corr.genome <- cor.test(x, y, method = "spearman")
y.aa <- filter(pairwise_distances.comp.bins, bin.genome.10th %in% c("1", "2", "3")) %>%
  pull(Sequence_similarity)
p.corr.genome.but.aa <- cor.test(x, y.aa, method = "spearman")
length(x)

x <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(Sequence_similarity)
p.corr.aa <- cor.test(x, y, method = "spearman")
y.genome <- filter(pairwise_distances.comp.bins, bin.aa.10th %in% c("1", "2", "3")) %>%
  pull(genomic_dist)
p.corr.aa.but.genome <- cor.test(x, y.genome, method = "spearman")
length(x)

x <- filter(pairwise_distances.comp.bins, bin.genome.4genes %in% c("1", "2", "3")) %>%
  pull(dist.trans)
y <- filter(pairwise_distances.comp.bins, bin.genome.4genes %in% c("1", "2", "3")) %>%
  pull(genomic_dist)
p.corr.genome.4genes <- cor.test(x, y, method = "spearman")
length(x)
pairwise_distances.comp.bins$bin.aa.10th<- factor(pairwise_distances.comp.bins$bin.aa.10th,levels=as.character(1:10))
p.binning.aa <- ggplot(pairwise_distances.comp.bins, 
                       aes(x=as.character(bin.aa.10th), y=dist.trans)) +
  geom_violin(aes(group=bin.aa.10th), color = "dodgerblue1", fill = "dodgerblue1", alpha = 0.4, scale = "width") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..), alpha = 0.5) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..), alpha = 0.5) +
  stat_summary(fun = "median", geom = "point", alpha = 0.8) +
  scale_x_discrete(limits = c(as.character(1:10))) +
 # scale_y_continuous(limits = c(0, 38), expand = c(0,0)) +
  common_layout +
  xlab("amino acid similarity (bin #)") +
  ylab("transcriptomic distance")

、
layout <- "
AABB
CCD#
"

p.binning <- wrap_plots(p.binning.genome, p.binning.prop, p.binning.aa, p.close_distant, design = layout)

ggsave(filename = "19_without_OR/fig2efg_binning.pdf", 
       plot = p.binning,
       device = "pdf", 
       units = "cm",
       width = 20, 
       height = 10, 
       useDingbats=FALSE)
