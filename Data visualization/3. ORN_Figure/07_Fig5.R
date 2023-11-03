# Fig5A plot the tree and label 39-42:

colors_for_exp_pattern<- c("#476D87","#E95C59")
DefaultAssay(ORN) <- "integratedRNA_onecluster"
object<- ORN
DefaultAssay(object)<-"RNA"
obj<-FindVariableFeatures(object, selection.method = "vst")
top <- head(VariableFeatures(obj),500)
DefaultAssay(obj)<-"raw_RNA"
obj<-ScaleData(obj,rownames(obj))
DefaultAssay(obj)<-"integratedRNA_onecluster"
object <- RunPCA(obj,features=c(all_receptor_gene,top),reduction.name="obj_features_pca") 
embeddings <- Embeddings(object = object, reduction = "obj_features_pca")[,1:50]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
ORN<- object
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])

multiOR_cluster<- cell_group[multiOR_cluster]

library(ggtree)
pdf("./00_Figure/Fig2/Fig2C-1-remove_nopower-cluster-ORN-tree-cosine.pdf",width=5,height=14)
tree <- groupOTU(data.tree, .node=multiOR_cluster)
ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
dev.off()

# Fig2C-2: chemoreceptor gene heatmap
m<-ggtree(tree,aes(color=group)) + geom_tiplab()+ geom_treescale() + scale_color_manual (values =colors_for_exp_pattern) 
cluster_order<-na.omit(m$data[order(m$data$y),]$label)
cluster_order<-as.character(cluster_order)
Idents(ORN)<-factor(ORN$cell_group,levels=cluster_order)



