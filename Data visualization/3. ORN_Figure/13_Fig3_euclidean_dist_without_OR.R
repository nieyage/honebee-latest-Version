# Fig3 the euclidean dist among 4 clusters
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


obj<- readRDS("./00_Figure/Fig4/Fig4-last-data-obj.rds")
Idents(obj)<-obj$group_manully
obj_features<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")

# Find the DEG among the 4 cluster:
DefaultAssay(obj)<- "SCT"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
markers<- markers[!duplicated(rownames(markers)),]$gene
cluster_DEG_withoutOR <- markers[-which(markers%in%obj_features)]

DefaultAssay(obj) <- "integratedRNA_onecluster";
obj<- ScaleData(obj,features=cluster_DEG_withoutOR)
obj <-  RunPCA(object = obj,features=cluster_DEG_withoutOR,reduction.name = "pca_withoutOR")

# Select PCA
pdf("./C1234_ElbowPlot_withoutOR_DEG_PCA.pdf")
ElbowPlot(obj,ndims = 50, reduction = "pca_withoutOR")
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
pca.coord.without.or <- obj@reductions$pca_withoutOR@cell.embeddings[,1:20]

dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)
colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)
osn.idents <- obj$group_manully

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

saveRDS(dist.mat.without.or,"C1234eucl_dist_matrix_without_or.rds")


## Permutations

# Create a matrix containing the first 20 PC.
pca.coords <- obj@reductions$pca@cell.embeddings[,1:20]
osn.idents <- obj$group_manully

# Determine a number of permutations
permutations <- paste("permutation_",
                      1:300,
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

saveRDS(permuted.distances,file = "C1234eucl_dist_permuations_without_or.rds")

# Read the data
permuted.distances <- read_rds("C1234eucl_dist_permuations_without_or.rds")
dist.mat.without.or <- read_rds("C1234eucl_dist_matrix_without_or.rds")

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
  geom_hline(yintercept = 0.5,
             linetype="dotted",color="red")+
  ylim(-0.1,1)+
  ylab("Cohen's D")+
  annotate(geom = "text",label ="***",y=1,x=1.5)+
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
  ylim(0,20)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))


euclidean.dist.plot.all <- wrap_plots(A = euclidean.dist.noOR.plot,
                                      B = cohend_all_clusters,
                                      design = "AAB")

ggsave(euclidean.dist.plot.all,
       filename = "./00_Figure/Fig3_among_cluster_euclidean_dist.pdf",
       width = 4,
       height = 2)


# In cluster 2:
tmp<- subset(obj,idents="C2")
tmp_data<- tmp@assays$raw_RNA[obj_features,]
Cbin<- ifelse(as.matrix(tmp_data) > 0, 1, 0)
exp_gene_number<- colSums(Cbin)
exp_gene_number<- exp_gene_number[colnames(tmp)]
tmp$gene_count <- exp_gene_number


##################################
#                                #        
#   ED all cells without OR      #
#                                #
##################################

# Compute the euclidean distance 
pca.coord.without.or <- tmp@reductions$pca_withoutOR@cell.embeddings[,1:10]

dist.mat.without.or <- fields::rdist(x1 = pca.coord.without.or, x2 = NULL)
colnames(dist.mat.without.or) <- rownames(pca.coord.without.or)
rownames(dist.mat.without.or) <- rownames(pca.coord.without.or)
osn.idents <- tmp$gene_count

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
## Permutations
# Create a matrix containing the first 20 PC.
pca.coords <- tmp@reductions$pca@cell.embeddings[,1:10]
osn.idents <- tmp$gene_count

# Determine a number of permutations
permutations <- paste("permutation_",
                      1:50,
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
# Cohen's D
same.osn.dist <- dist.mat.without.or$euclidean_distance[dist.mat.without.or$osn_pairs == "same OSN pop."]
different.osn.dist <- unlist(dist.mat.without.or$euclidean_distance[dist.mat.without.or$osn_pairs == "different OSN pop."])
permutation <- lapply(permuted.distances, '[[', 3)
intra.inter <- cohen.d(same.osn.dist,sample(different.osn.dist))
density.cohend <- lapply(c(1:50), FUN = function(x){ 
  p <- unlist(permutation[x])
  sp <- cohen.d(same.osn.dist,p)
  sp <- sp$estimate
  dp <- cohen.d(different.osn.dist,p)
  dp <- dp$estimate
  c(sp,dp)
})


## Checkpoint ##
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
  geom_hline(yintercept = 0.5,
             linetype="dotted",color="red")+
  ylim(-0.1,1)+
  ylab("Cohen's D")+
  annotate(geom = "text",label ="ns",y=1,x=1.5)+
  common_layout+
  theme(axis.title.x = element_blank())

# Combine the permuted matrix
permuted.distances <- dplyr::bind_rows(permuted.distances)
permuted.distances$osn_pairs <- "permutation"
# Combine the distance matrix to the permuted distance matrix
 median(dist.mat.without.or$euclidean_distance)
  median(permuted.distances$euclidean_distance)

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
  ylim(0,30)+
  common_layout +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0,20),
        legend.background = element_rect(fill='transparent'))


euclidean.dist.plot.all <- wrap_plots(A = euclidean.dist.noOR.plot,
                                      B = cohend_all_clusters,
                                      design = "AAB")

ggsave(euclidean.dist.plot.all,
       filename = "./00_Figure/Fig3_cluster2_euclidean_dist_withoutOR.pdf",
       width = 4,
       height = 2)






