# cluster DEG and DEP 
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
DefaultAssay(ORN)<- "raw_RNA"

# Step1: tau cluster specific DEG 
library(VGAM)
ORN_matrix<-as.matrix(GetAssayData(ORN));
#filter the gene only appear in a few cells 
cell_pct = function(data){
    pct<-length(data[data!=0])/length(data)
    return(pct)
}
gene_pct<-apply(ORN_matrix,1,cell_pct)
gene_pass_pct<-names(gene_pct[gene_pct>0.005])
# gene_pass_pct_rmreceptor<- gene_pass_pct[-which(gene_pass_pct%in%all_receptor_gene)]

ORN<-NormalizeData(ORN)
ORN_avg<-AverageExpression(
       ORN,
       assays = "raw_RNA",
       features = gene_pass_pct,
       return.seurat = FALSE,
       group.by = "cell_group",
       #add.ident = NULL,
       slot = "data")
ORN_avg<-ORN_avg$raw_RNA;
#https://fmicompbio.github.io/swissknife/reference/specificityScore-methods.html#value-1
source("/md01/nieyg/software/swissknife/R/tissue_specificity_score.R")
library(matrixStats)
gene_tau<-specificityScore(
  ORN_avg,
  method = c("tau", "TSI", "counts"),
  #group = ORN$subcluster,
  thresh = 0,
  expr_values = "logcounts",
  na.rm = FALSE
)
names(gene_tau)<-rownames(ORN_avg)
# plot the density plot for gene_tau 

pdf("./05_ORN_cluster2/07_DEG_and_DEP/03_rm_receptor/01_DEG_by_tau/gene_tau_distribution.pdf",width=10,height=5)
data<- as.data.frame(gene_tau)
ggplot(data, aes(x=data[,1])) + xlab("")+
              geom_density(alpha=.25) + theme_classic() 
d <- density(data[,1])
d$x[which.min(abs(diff(d$y)))]
hist(data[,1],prob=TRUE)
lines(d, col="red", lty=2)
#Kmeans
dev.off()

#plot the cluster specific pheatmap
gene_specific<-names(which(gene_tau>0.8))
gene_specific_data<-as.data.frame(ORN_avg[gene_specific,])
data<-data.frame()
for (gene in gene_specific){
    gene_cluster_avg<-gene_specific_data[gene,]
    gene_specific_cluster<-names(gene_cluster_avg[which(gene_cluster_avg==max(gene_cluster_avg))])
    data_subset<-data.frame(gene,cluster=gene_specific_cluster);
    data<-rbind(data,data_subset)
}
data$cluster<-factor(data$cluster,levels=colnames(gene_specific_data))
data<-data[order(data$cluster),]
tau1_gene<-data$gene
write.csv(data,"./05_ORN_cluster2/07_DEG_and_DEP/03_rm_receptor/01_DEG_by_tau/DEG_tau-0.8_cluster_specfic_data.csv")
#gene_specific_data<-gene_specific_data[do.call(order,gene_specific_data),]
# plot the tau cluster specfic genes heatmap 
#modify_tau_gene<-names(which(gene_tau>0.945))
tau_Avg <- AverageExpression(ORN,features=tau1_gene,assays = "raw_RNA")
library(pheatmap)
count=t(scale(t(tau_Avg$raw_RNA),scale = T,center = T))
pdf("./19_without_OR/DEG_tau_heatmap_withOR.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows /= F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
dev.off();
cluster_DEG_withOR<- tau1_gene
cluster_DEG_withoutOR<-tau1_gene[-which(tau1_gene%in%all_receptor_gene)]

tau_Avg <- AverageExpression(ORN,features=cluster_DEG_withoutOR,assays = "raw_RNA")
library(pheatmap)
count=t(scale(t(tau_Avg$raw_RNA),scale = T,center = T))

pdf("./19_without_OR/DEG_tau_heatmap_withoutOR.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=T,show_colnames=T)
dev.off();


# Step2: RunPCA by DEG(with OR or without OR )
library(reticulate)
py_install("kneed")
py_install("matplotlib")

DefaultAssay(ORN) <- "integratedRNA_onecluster";
ORN<- ScaleData(ORN,features=cluster_DEG_withOR)
ORN <-  RunPCA(object = ORN,features=cluster_DEG_withOR,reduction.name = "pca_withOR")

# Select PCA
pdf("./19_without_OR/ElbowPlot_withOR_DEG_PCA.pdf")
ElbowPlot(ORN,ndims = 50, reduction = "pca_withOR")
dev.off()

DefaultAssay(ORN) <- "integratedRNA_onecluster";
ORN<- ScaleData(ORN,features=cluster_DEG_withoutOR)
ORN <-  RunPCA(object = ORN,features=cluster_DEG_withoutOR,reduction.name = "pca_withoutOR")

# Select PCA
pdf("./19_without_OR/ElbowPlot_withoutOR_DEG_PCA.pdf")
ElbowPlot(ORN,ndims = 50, reduction = "pca_withoutOR")
dev.off()

# plot the transcript dist heatmap 

# cell cosine simility heatmap 
library(lsa)
embeddings <- Embeddings(object = ORN, reduction = "pca_withOR")[,1:20]
metadata<- data.frame(cell_group=ORN$cell_group)
rownames(metadata) <- colnames(ORN)

data<- data.frame(barcode= colnames(ORN),cell_group=ORN$cell_group)
data$cell_group<- factor(data$cell_group,levels=levels(ORN))
data<- data[order(data$cell_group),]

embeddings <- embeddings[data$barcode,]
trans_dist <- 1-cosine(t(embeddings))

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )

col<-myUmapcolors[1:length(unique(metadata$cell_group))]
names(col)<-unique(metadata$cell_group)
ann_colors= list(cell_group = col)
pdf("./19_without_OR/transcript_dist_heatmap_withOR.pdf",width=10,height=10)
pheatmap(trans_dist,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = metadata,
         annotation_colors = ann_colors,
         annotation_row = metadata,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
dev.off()

embeddings <- Embeddings(object = ORN, reduction = "pca_withoutOR")[,1:20]
embeddings <- embeddings[data$barcode,]
trans_dist <- 1-cosine(t(embeddings))
col<-myUmapcolors[1:length(unique(metadata$cell_group))]
names(col)<-unique(metadata$cell_group)
ann_colors= list(cell_group = col)
pdf("./19_without_OR/transcript_dist_heatmap_withoutOR.pdf",width=10,height=10)

pheatmap(trans_dist,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
         annotation_col = metadata,
         annotation_colors = ann_colors,
         annotation_row = metadata,
         annotation_legend = TRUE,
         show_rownames=F,
         show_colnames=F
  )
dev.off()

##################################
#                                #        
#   ED all cells without OR      #
#                                #
##################################

library(Matrix)
library(patchwork)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(effsize)
# Compute the euclidean distance 
pca.coord.without.or <- ORN@reductions$pca_withoutOR@cell.embeddings[,1:20]

dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)

colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)

osn.idents <- ORN$cell_group

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

saveRDS(dist.mat.without.or,"eucl_dist_matrix_without_or.rds")


## Permutations

# Create a matrix containing the first 15 PC.
pca.coords <- ORN@reductions$pca@cell.embeddings[,1:15]
osn.idents <- ORN$cell_group

# Determine a number of permutations
permutations <- paste("permutation_",
                      1:1000,
                      sep = "")
names(x = permutations) <- permutations

permuted.distances <- lapply(X = permutations,
                             FUN = function(permutation,
                                            pca.coords,
                                            osn.idents) {
                               print(x = permutation)
                               # Compute euclidean distances of PCA eigenvectors after permuting OSN identities
                               rownames(x = pca.coords) <- sample(x = rownames(x = pca.coords),
                                                                  size = nrow(x = pca.coords),
                                                                  replace = FALSE)
                               
                               dist.mat <- fields::rdist(x1 = pca.coords,
                                                         x2 = NULL)
                               
                               colnames(dist.mat) <- rownames(pca.coords)
                               rownames(dist.mat) <- rownames(pca.coords)
                               
                               dist.mat[lower.tri(dist.mat, diag = TRUE)] <- NA
                               
                               dist.mat <- melt(data = dist.mat,
                                                na.rm = TRUE,
                                                varnames = c("OSN_1", "OSN_2"),
                                                value.name = "euclidean_distance")
                               
                               dist.mat <- data.table::data.table(dist.mat)
                               dist.mat$OSN_1 <- as.character(x = dist.mat$OSN_1)
                               dist.mat$OSN_2 <- as.character(x = dist.mat$OSN_2)
                               
                               dist.mat$OSN_1_ident <- osn.idents[dist.mat$OSN_1]
                               dist.mat$OSN_2_ident <- osn.idents[dist.mat$OSN_2]
                               
                               dist.mat$osn_pairs <- ifelse(test = dist.mat$OSN_1_ident == dist.mat$OSN_2_ident,
                                                            yes = "same OSN pop.",
                                                            no = "different OSN pop.")
                               
                               dist.mat <- dist.mat[dist.mat$osn_pairs == "same OSN pop.",]
                               
                               return(dist.mat)
                             },
                             pca.coords = pca.coords,
                             osn.idents = osn.idents)

#%%% Checkpoint %%%#

saveRDS(permuted.distances,file = "eucl_dist_permuations_without_or.rds")

# Read the data
permuted.distances <- read_rds("eucl_dist_permuations_without_or.rds")
dist.mat.without.or <- read_rds("eucl_dist_matrix_without_or.rds")

# Cohen's D
same.osn.dist <- dist.mat.without.or$euclidean_distance[dist.mat.without.or$osn_pairs == "same OSN pop."]
different.osn.dist <- unlist(dist.mat.without.or$euclidean_distance[dist.mat.without.or$osn_pairs == "different OSN pop."])
permutation <- lapply(permuted.distances, '[[', 3)

intra.inter <- cohen.d(same.osn.dist,sample(different.osn.dist))


density.cohend <- lapply(c(1:1000), FUN = function(x){
  
  p <- unlist(permutation[x])
  
  sp <- cohen.d(same.osn.dist,p)
  sp <- sp$estimate
  
  dp <- cohen.d(different.osn.dist,p)
  dp <- dp$estimate
  
  c(sp,dp)
  
})


## Checkpoint ##

saveRDS(density.cohend,file = "cohendsD_density_without_or.rds")

density.cohend <- readRDS("cohendsD_density_without_or.rds")

density.cohend.dframe <- data.frame(t(data.frame(density.cohend)))
a <- data.frame(dist = density.cohend.dframe$X1, type = "intra_perm")
b <- data.frame(dist = density.cohend.dframe$X2, type = "inter_perm")

density.cohend.dframe <- rbind(a,b)
density.cohend.dframe$dist <- density.cohend.dframe$dist*c(-1)

stat.test <- ks.test(density.cohend.dframe$dist[density.cohend.dframe$type == "intra_perm"],
                     density.cohend.dframe$dist[density.cohend.dframe$type == "inter_perm"])

common_layout <- theme_classic(base_size = 10, base_family = "ArialMT") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.7, "lines")
  )

cohend_all_clusters <- ggplot(density.cohend.dframe,aes(y = dist,
                                                        x= type))+ 
  geom_violin(scale = "width")+ 
  geom_hline(yintercept = -1*intra.inter$estimate,
             linetype="dotted",color="red")+
  ylim(-0.3,1.7)+
  ylab("Cohen's D")+
  annotate(geom = "text",label ="***",y=1.7,x=1.5)+
  common_layout+
  theme(axis.title.x = element_blank())


# Combine the permuted matrix
permuted.distances <- dplyr::bind_rows(permuted.distances)
permuted.distances$osn_pairs <- "permutation"


# Combine the distance matrix to the permuted distance matrix
dist.mat.without.or.per <- dplyr::bind_rows(dist.mat.without.or, permuted.distances)
dist.mat.without.or.per$osn_pairs <- factor(dist.mat.without.or.per$osn_pairs,
                                            levels = c("same OSN pop.",
                                                       "different OSN pop.",
                                                       "permutation"))

levels(dist.mat.without.or.per$osn_pairs) <- c("intra.","inter.","perm.")

euclidean.dist.noOR.plot <- ggplot(data = dist.mat.without.or.per,
                                   mapping = aes(y = euclidean_distance,
                                                 x = osn_pairs,
                                                 fill = osn_pairs,
                                                 group = osn_pairs)) +
  geom_violin(alpha = 0.5,scale = "width",color="gray50",fill = "gray50") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "eucl. dist.") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))


##################################
#                                #        
#   ED all cells with OR         #
#                                #
##################################

# Compute the euclidean distance 
pca.coord.with.or <- ORN@reductions$pca_withOR@cell.embeddings[,1:20]

dist.mat.with.or <- fields::rdist(x1 = pca.coord.with.or, x2 = NULL)

colnames(dist.mat.with.or) <- rownames(pca.coord.with.or)
rownames(dist.mat.with.or) <- rownames(pca.coord.with.or)

osn.idents <- ORN$cell_group

dist.mat.with.or[lower.tri(dist.mat.with.or, diag = TRUE)] <- NA

dist.mat.with.or <- melt(data = dist.mat.with.or,
                         na.rm = TRUE,
                         varnames = c("OSN_1", "OSN_2"),
                         value.name = "euclidean_distance")

dist.mat.with.or <- data.table::data.table(dist.mat.with.or)
dist.mat.with.or$OSN_1 <- as.character(x = dist.mat.with.or$OSN_1)
dist.mat.with.or$OSN_2 <- as.character(x = dist.mat.with.or$OSN_2)

dist.mat.with.or$OSN_1_ident <- osn.idents[dist.mat.with.or$OSN_1]
dist.mat.with.or$OSN_2_ident <- osn.idents[dist.mat.with.or$OSN_2]

dist.mat.with.or$osn_pairs <- ifelse(test = dist.mat.with.or$OSN_1_ident == dist.mat.with.or$OSN_2_ident,
                                     yes = "same OSN pop.",
                                     no = "different OSN pop.")

dist.mat.with.or$osn_pairs <- factor(dist.mat.with.or$osn_pairs,
                                     levels = c("same OSN pop.","different OSN pop."))



# Cohen's D
same.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "same OSN pop."]
different.osn.dist <- dist.mat.with.or$euclidean_distance[dist.mat.with.or$osn_pairs == "different OSN pop."]

same_diff <- cohen.d(same.osn.dist,different.osn.dist)


levels(dist.mat.with.or$osn_pairs) <- c("intra.","inter.")

# Plots
euclidean.dist.or.plot <- ggplot(data = dist.mat.with.or,
                                 mapping = aes(y = euclidean_distance,
                                               x = osn_pairs,
                                               fill = osn_pairs,
                                               group = osn_pairs)) + 
  geom_violin(alpha = 0.5,scale = "width",color="gray50",fill = "gray50") +
  stat_mean(geom="segment", mapping = aes(xend=..x.. + 0.5, yend=..y..)) +
  stat_mean(geom="segment", mapping = aes(xend=..x.. - 0.5, yend=..y..)) +
  stat_summary(geom = "point", fun = "median")+
  labs(y = "eucl. dist.") + 
  ylim(0,70)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))



euclidean.dist.plot.all <- wrap_plots(A = euclidean.dist.noOR.plot,
                                      B = cohend_all_clusters,
                                      C = euclidean.dist.or.plot,
                                      design = "AABC")

ggsave(euclidean.dist.plot.all,
       filename = "fig1_E_eucldist_all.pdf",
       width = 6,
       height = 2)