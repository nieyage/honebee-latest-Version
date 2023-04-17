library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(RColorBrewer)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
set.seed(1234)
honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype_Neuron_subcluster.rds")

# find the DEG among all celltypes and get the specifical TF in honeybees
Idents(honeybee)<- honeybee$Annotation_subcluster
DefaultAssay(honeybee)<-"peaks"
da_peaks <- FindAllMarkers(
  object = honeybee,
  test.use = 'LR',
  logfc.threshold = 0.1,
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
#markers <- FindAllMarkers(honeybee, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
table(da_peaks$cluster)
write.csv(da_peaks,"./02_All_celltype/DEP_FindAllMarkers_all_celltype.csv")
#verification
peak2show<- rownames(da_peaks)
honeybee<-ScaleData(honeybee,features=rownames(honeybee))

library(ArchR)
pdf("./02_All_celltype/DEP_FindAllMarkers_all_celltype_heatmap.pdf",width=20,height=10)
DoHeatmap(object = honeybee,features=peak2show,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = honeybee,features=peak2show,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
DoHeatmap(object = honeybee,features=peak2show,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = ArchRPalettes$solarExtra)+NoLegend()
dev.off();

# motif enrichment 
# Get a list of motif position frequency matrices from the JASPAR database

library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'insects', all_versions = FALSE)
)
DefaultAssay(honeybee) <- 'peaks'
# add motif information
honeybee <- AddMotifs(
  object = honeybee,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = pfm
)





pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)



All_motif_info <- data.frame()
for (cluster in levels(honeybee)){
  cluster_peak <- markers[markers$cluster==cluster,]$gene;
  enriched.motifs <- FindMotifs(honeybee,features = cluster_peak);
  enriched.motifs$cluster <- cluster;
  All_motif_info <- rbind(All_motif_info,enriched.motifs)
}
library(dplyr)
top3 <- All_motif_info %>% group_by(cluster) %>% top_n(n = 3, wt = fold.enrichment)
motif2show<-unique(top3$motif.name)

# motif enrich matrix
motif_matrix<-matrix(ncol=length(motif2show),nrow=length(levels(honeybee)))
colnames(motif_matrix)<-motif2show
rownames(motif_matrix)<-levels(honeybee)

last_motif_info<-All_motif_info[which(All_motif_info$motif.name%in%motif2show),]
for (i in 1:nrow(last_motif_info)){
  cluster=last_motif_info[i,]$cluster;
  motif=last_motif_info[i,]$motif.name;
  motif_matrix[cluster,motif]<-last_motif_info[i,]$fold.enrichment
}

library(pheatmap)
#count=t(scale(t(motif_matrix),scale = T,center = T))
pdf("./02_All_celltype/FindAllMarker_DEP_motif_heatmap.pdf",width=20,height=20)
pheatmap(motif_matrix,cluster_cols = T,cluster_rows = F,
              color = colorRampPalette(c("white", "firebrick3"))(100),
              cellwidth = 10, cellheight = 10,
              show_rownames=T,show_colnames=T)
dev.off()

#cluster_closest_open <- ClosestFeature(ORN, cluster_peak)

enriched.motifs <- FindMotifs(
  object = ORN,
  features = cluster_peak
)
MotifPlot(
  object = ORN,
  motifs = head(rownames(enriched.motifs))
)
honeybee <- RunChromVAR(
  object = honeybee,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor
)

p2 <- FeaturePlot(
  object = ORN,
  features = head(rownames(enriched.motifs)),
  min.cutoff = 'q10',
  )

DefaultAssay(ORN) <- 'chromvar'

differential.activity <- FindAllMarkers(
  object = ORN,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = ORN,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)


# Motif footprinting
# Now we can footprint any motif that we have positional information for. By default, this includes every instance of the motif in the genome. We can instead use the in.peaks = TRUE parameter to include only those motifs that fall inside a peak in the assay. The Footprint() function gathers all the required data and stores it in the assay. 
# We can then plot the footprinted motifs using the PlotFootprint() function.
# gather the footprinting information for sets of motifs
DefaultAssay(ORN) <- 'peaks'

# add motif information
ORN <- AddMotifs(
  object = ORN,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = pfm
)

ORN <- Footprint(
  object = ORN,
  motif.name =head(rownames(differential.activity)),
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(ORN, features = head(rownames(differential.activity))

open_peak <- da_peaks[da_peaks$avg_log2FC > 0.25, ]$gene
close_peak <- da_peaks[da_peaks$avg_log2FC < -0.25, ]$gene
closest_open <- ClosestFeature(ORN, unique(open_peak))
closest_close <- ClosestFeature(ORN, close_peak)

# TF exp 

TF<-c("LOC552100",#ovo
  "LOC410499",#prd
  "LOC411079",#grh
  "Usp","Dl",
  "LOC411207",#slbo
  #"LOC10057220",#ro(rough)
  "LOC100576147", "LOC724740"  ,  "LOC725966"  ,  "LOC726165"#fkh
  )
label<-c("ovo",#ovo
  "prd",#prd
  "grh",#grh
  "Usp","Dl",
  "slbo",#slbo
  #"ro",#ro(rough)
  "fkh FD4", "slp2"  ,  "fkh crocodile"  ,  "slp1"#fkh
  )
DefaultAssay(ORN)<-"raw_RNA"
pdf("./ORN/remove_nopower/TF-ORN-exp.pdf",width=6,height=12)
p<-DotPlot(ORN, features = TF,dot.scale = 3)+theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p
p&scale_x_discrete(labels=label)
dev.off()
 




