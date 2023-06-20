#Identify DEP by FindMarkers, Kendall tau,correlation and mutual information 

library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)
set.seed(1234);
library(pheatmap)
# remove nopower cluster,then plot UMAP
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )
# 1. Find All Markers:
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
da_peaks <- FindAllMarkers(
  object = ORN,
  test.use = 'LR',
  #logfc.threshold = 0.1,
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
write.csv(da_peaks,"./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_FindAllMarkers.csv")
#markers <- FindAllMarkers(ORN, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
da_peaks<-da_peaks[da_peaks$p_val_adj<0.05,]
da_peaks<- da_peaks[da_peaks$avg_log2FC>1,]
table(da_peaks$cluster)
#verification
peak2show<- rownames(da_peaks)


# plot the avg heatmap 
peak_Avg <- AverageExpression(ORN,features=peak2show,assays = "peaks")
count=t(scale(t(peak_Avg$peaks),scale = T,center = F))
pdf("./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_FindAllMarkers_ORN_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();

# 2. tau cluster specific:
library(VGAM)
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
ORN_matrix<-as.matrix(GetAssayData(ORN));
#filter the gene only appear in a few cells 
cell_pct = function(data){
    pct<-length(data[data!=0])/length(data)
    return(pct)
}
peak_pct<-apply(ORN_matrix,1,cell_pct)
peak_pass_pct<-names(peak_pct[peak_pct>0.005])
ORN<-NormalizeData(ORN)
ORN_avg<-AverageExpression(
       ORN,
       assays = "peaks_ORN_subcluster",
       features = peak_pass_pct,
       return.seurat = FALSE,
       slot = "counts")
ORN_avg<-ORN_avg$peaks_ORN_subcluster
#https://fmicompbio.github.io/swissknife/reference/specificityScore-methods.html#value-1
source("/md01/nieyg/project/honeybee/add_antenna/swissknife-master/R/tissue_specificity_score.R")
library(matrixStats)
peak_tau<-specificityScore(
  ORN_avg,
  method = c("tau", "TSI", "counts"),
  thresh = 0,
  expr_values = "logcounts",
  na.rm = FALSE
)
names(peak_tau)<-rownames(ORN_avg)
# plot the density plot for peak_tau 
pdf("./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_tau_density.pdf",width=10,height=5)
data<- as.data.frame(peak_tau)
ggplot(data, aes(x=data[,1])) + xlab("")+
              geom_density(alpha=.25) + theme_classic() 
d <- density(data[,1])
d$x[which.min(abs(diff(d$y)))]
hist(data[,1],prob=TRUE)
lines(d, col="red", lty=2)
dev.off()

#plot the cluster specific pheatmap
DefaultAssay(ORN)<- "peaks_ORN_subcluster"
peak_specific<-names(which(peak_tau>0.92))
colnames(ORN_avg) <- levels(ORN)
peak_specific_data<-as.data.frame(ORN_avg[peak_specific,])

data<-data.frame()
for (peak in peak_specific){
    peak_cluster_avg<-peak_specific_data[peak,]
    peak_specific_cluster<-names(peak_cluster_avg[which(peak_cluster_avg==max(peak_cluster_avg))])
    data_subset<-data.frame(peak,cluster=peak_specific_cluster);
    data<-rbind(data,data_subset)
}
data$cluster<-factor(data$cluster,levels=levels(ORN))
data<-data[order(data$cluster),]
tau1_peak<-data$peak
write.csv(data,"./05_ORN_cluster2/07_DEG_and_DEP/02_without_GRN/DEP_tau-0.9_cluster_specfic_data_peak_ORN.csv")
# plot the tau cluster specfic genes heatmap 
# plot the avg heatmap 
peak_Avg <-ORN_avg[tau1_peak,]
count=t(scale(t(peak_Avg),scale = T,center = F))
count<- count[tau1_peak,]
pdf("./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_tau_ORN_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();

# 3. the correlation between ORx and other gene:
# library(dplyr)
# library(BBmisc)
# library(mlr)
# library(infotheo)
# peak_CV<-apply(t_ORN_matrix, 2, cal_cv)
# 
# # density plot 
# data<-as.data.frame(peak_CV);
# range(data$peak_CV)
# #[1]  0.4903281 40.3753042
# library(ggplot2)
# pdf("./ORN/remove_nopower/DEGandDEP/allpeak_cv_ATAC.pdf",width=10,height=4)
# ggplot(data, aes(x=peak_CV)) + xlab("peak_CV")+
#               geom_density(alpha=.25) + theme_classic() 
# dev.off()
# data$peak<-rownames(data)
# #min(data[data$peak%in%peakname,]$peak_CV)
# #[1] 4.769134
# peak_hCV<-data[data$peak_CV>4.5,2]
# 
# # the OR gene promoter region peak 
# library(tidyr)
# peak_filtered<-separate(as.data.frame(peak_hCV),"peak_hCV",c("Group","start","end"),"-")
# # create GRange object 
# gr <- GRanges(
#     seqnames = Rle(peak_filtered$Group),
#     ranges = IRanges(start= as.numeric(peak_filtered$start),
#     	             end = as.numeric(peak_filtered$end), 
#     	             names = peak_hCV)
#     )
# 
# #all_peak_in_ourdata<-granges(ORN)
# library(rtracklayer)
# dotplot_data<-read.csv("./ORN/remove_nopower/dotplot_data_remove_nopower.csv",row.names=1)
# dotplot_feature<- unique(dotplot_data$features.plot)
# tss =  import("/md01/nieyg/ref/10X/honeybee/de_novo_antenna/Amel_antenan/regions/tss.bed")
# OR_ID<-unique(Annotation(ORN)[which(Annotation(ORN)$gene_name%in%dotplot_feature),]$gene_id)
# OR_tss<-tss[tss$name%in%OR_ID,]
# #OR_gene_promoter<-promoters(OR_tss, upstream=1000, downstream=1000, use.names=TRUE)
# gene_id2name<-unique(as.data.frame(Annotation(ORN))[,c(10,12)])
# #OR_gene_promoter$name<-gene_id2name[match(OR_gene_promoter$name,gene_id2name$gene_id),]$gene_name
# OR_tss$name <- gene_id2name[match(OR_tss$name,gene_id2name$gene_id),]$gene_name
# OR_tss_nearest<-gr[nearest(OR_tss, gr),]
# OR_tss_nearest$name<- OR_tss$name
# #OR_tss_nearest_peak<-as.data.frame(OR_tss_nearest)[,1:3]
# #peakname<-do.call(paste, c(OR_tss_nearest_peak[ ,1:3], list(sep = '-')))
# peakname<-names(OR_tss_nearest)
# peak2calculate<-intersect(peak_pass_pct,peak_hCV)
# peak2calculate<-c(peak2calculate,peakname)
# t_ORN_matrix<-t_ORN_matrix[,unique(peak2calculate)]
# 
# # Step: Feature Selection 
# # Correlation Matrix
# cor_data<-cor(t_ORN_matrix)
# Or_peak_cor_data<-cor_data[unique(peakname),]
# # Manage correalation results 
# #cluster_info<-as.data.frame(table(dotplot_data$id))
# #multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])
# #cluster_info$Var1<-factor(cluster_info$Var1,levels=levels(ORN))
# #cluster_info<-cluster_info[order(cluster_info$Var1),]
# Or_cor_data<-data.frame()
# for (i in levels(ORN)){
# 	cluster<-i;
# 	print(cluster);
# 	cluster_OR <- dotplot_data[dotplot_data$id==cluster,]$features.plot
# 	cluster_peak <- unique(names(OR_tss_nearest[OR_tss_nearest$name%in%cluster_OR,]))
#   if(length(cluster_peak)==1){
#     cor_value<-Or_peak_cor_data[cluster_peak,]
#     cor_value<-cor_value[order(cor_value,decreasing=T)];
#     top100_cor<-names(cor_value[1:100]);
#     OR_char<-paste(cluster_OR,collapse="_")
#     Or_cor_data_subset<-data.frame(OR=OR_char,peak=top100_cor,cluster);
#     Or_cor_data<-rbind(Or_cor_data,Or_cor_data_subset);}
#   if(length(cluster_peak)>1){
#     cor_value_all<-Or_peak_cor_data[cluster_peak,]
#     cor_value<-colSums(cor_value_all)
#     cor_value<-cor_value[order(cor_value,decreasing=T)]
#     top100_cor<-names(cor_value[1:100])
#     OR_char<-paste(cluster_OR,collapse="_")
#     Or_cor_data_subset<-data.frame(OR=OR_char,peak=top100_cor,cluster);
#     Or_cor_data<-rbind(Or_cor_data,Or_cor_data_subset);
#   }
#   if(length(cluster_peak)==0){
#   	p<-paste("cluster",cluster,"not have OR promoter peak!")
#   	print(p)
#   }
# }
# 
# write.csv(Or_cor_data,"./ORN/remove_nopower/DEGandDEP/DEP_correlation_top100_data_ATAC.csv")
# ORN<-ScaleData(ORN,features=Or_cor_data$peak)
# pdf("./ORN/remove_nopower/DEGandDEP/DEP-Correalation_heatmap-top100_peak_ORN.pdf",width=20,height=10)
# DoHeatmap(object = ORN,features=Or_cor_data$peak,label=TRUE,size = 3.5,group.colors =myUmapcolors) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
# DoHeatmap(object = ORN,features=Or_cor_data$peak,label=TRUE,size = 3.5,group.colors =myUmapcolors) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
# DoHeatmap(object = ORN,features=Or_cor_data$peak,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = ArchRPalettes$solarExtra)+NoLegend()
# dev.off();

# 4. the Mutinformation between ORx and other gene:
#library(infotheo)
#Or_mi_data<-data.frame()
#for (i in levels(ORN)){
#	cluster<-i;
#	print(cluster);
#	cluster_OR <- dotplot_data[dotplot_data$id==cluster,]$features.plot
#	cluster_peak <- unique(names(OR_tss_nearest[OR_tss_nearest$name%in%cluster_OR,]))
#  if(length(cluster_peak)==1){
#    mi_value<-data.frame() 
#    Or_peak_value<-t_ORN_matrix[,cluster_peak];
#    for (j in 1:ncol(t_ORN_matrix)){
#     mivalue<-mutinformation(Or_peak_value,t_ORN_matrix[,j]);
#     OR_char<-paste(cluster_OR,collapse="_")
#     mi_value_subset<-data.frame(Or_gene=OR_char,mi_value=mivalue,peak=colnames(t_ORN_matrix)[j],cluster)
#     mi_value<-rbind(mi_value,mi_value_subset)
#   };
#   mi_value<-mi_value[order(mi_value$mi_value,decreasing=T),]
#   top100_mi<-mi_value[1:100,];
#   Or_mi_data<-rbind(Or_mi_data,top100_mi);}
#  if(length(cluster_peak)>1){
#    mi_value<-data.frame() 
#    Or_peak_value<-t_ORN_matrix[,cluster_peak];
#    for (j in 1:ncol(t_ORN_matrix)){
#     mivalue<-mutinformation(Or_peak_value,t_ORN_matrix[,j])
#     OR_char<-paste(cluster_OR,collapse="_")
#     mi_value_subset<-data.frame(Or_gene=OR_char,mi_value=mivalue,gene=colnames(t_ORN_matrix)[j],cluster)
#     mi_value<-rbind(mi_value,mi_value_subset)
#   };
#   mi_value<-mi_value[order(mi_value$mi_value,decreasing=T),]
#   top100_mi<-mi_value[1:100,];
#   Or_mi_data<-rbind(Or_mi_data,top100_mi);
#  }
#}
#write.csv(Or_mi_data,"./ORN/remove_nopower/DEGandDEP/DEP_MI_top100_data.csv")
#pdf("./ORN/remove_nopower/DEGandDEP/DEP_MI_heatmap-top100.pdf",width=15,height=10)
#DoHeatmap(object = ORN,features=Or_mi_data$peak,label=TRUE,size = 3.5,group.colors =myUmapcolors) + scale_fill_gradientn(colors = c( "white", "#E41A1C"))+NoLegend()
#DoHeatmap(object = ORN,features=Or_mi_data$peak,label=TRUE,size = 3.5,group.colors =myUmapcolors) + scale_fill_gradientn(colors = c("#377EB8", "white", "#E41A1C"))+NoLegend()
#DoHeatmap(object = ORN,features=tau1_peak,label=TRUE,size = 3.5, group.colors =myUmapcolors) + scale_fill_gradientn(colors = ArchRPalettes$solarExtra)+NoLegend()
#dev.off();

# plot the overlap among 4 methods
library(VennDiagram)
library(RColorBrewer)
#markers <-read.csv("./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_FindAllMarkers_gene_peak_ORN.csv",row.names=1)
tau_data<-read.csv("./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_tau-0.9_cluster_specfic_data_peak_ORN.csv",row.names=1)
#cor_data<-read.csv("./ORN/remove_nopower/DEGandDEP/DEP_correlation_top100_data_peak_ORN.csv",row.names=1)
#mi_data <-read.csv("./ORN/remove_nopower/DEGandDEP/DEP_MI_top100_data.csv",row.names=1)
FeatureSelection_data<-data.frame()
pdf("./05_ORN_cluster2/07_DEG_and_DEP/01_All_ORN_cluster/DEP_2methods_venn..pdf")
for (i in levels(ORN)){
    cluster<-i;
    # 1.FindAllMarkers:
    findmarkers<-markers[markers$cluster==cluster,]$gene;
    # 2.tau:
    tau<-tau_data[tau_data$cluster==cluster,]$peak;
    # 3.correlation:
    cor<-cor_data[cor_data$cluster==cluster,]$peak;
    # 4.MI:
    #mi<-mi_data[mi_data$cluster==cluster,]$gene;
    Or_gene<-dotplot_data[dotplot_data$id==cluster,]$features.plot
    OR_char<-paste(Or_gene,collapse="_")
    title<-paste(cluster,OR_char,sep=": ")
    vennplot<-venn.diagram(
        x = list(findmarkers,tau,cor),
        category.names = c("findmarkers","tau","top100_cor"),
        filename =NULL,
        fill = brewer.pal(7, "Set2")[1:3],
        alpha = 0.50,
        main=title,
        output=TRUE
      )
    print(grid.draw(vennplot))
    grid.newpage();
}
dev.off()


# 3 methods upsetR
library(UpSetR)
listInput <- list(
        #MI = mi_data$gene, 
        tau = tau_data$peak, 
        FindAllMarkers = markers$gene,
        correlation = cor_data$peak)
pdf("./ORN/remove_nopower/DEGandDEP/DEP_3methods_upsetR.pdf")
upset(fromList(listInput), nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
dev.off()


# Motif enrichment for each OR cluster 

# For FindAllMarkers:
# test enrichment
# Get a list of motif position frequency matrices from the JASPAR database

library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'insects', all_versions = FALSE)
)
DefaultAssay(ORN) <- 'peaks_ORN_subcluster'
# add motif information
ORN <- AddMotifs(
  object = ORN,
  genome = BSgenome.Amel.antenan,
  pfm = pfm
)

# tau 

All_motif_info <- data.frame()
for (cluster in levels(ORN)){
  cluster_peak <- tau_data[tau_data$cluster==cluster,]$peak;
  enriched.motifs <- FindMotifs(ORN,features = cluster_peak);
  enriched.motifs$cluster <- cluster;
  All_motif_info <- rbind(All_motif_info,enriched.motifs)
}  
library(dplyr)
top3 <- All_motif_info %>% group_by(cluster) %>% top_n(n = 3, wt = fold.enrichment)
motif2show<-unique(top3$motif.name)

# motif enrich matrix

motif_matrix<-matrix(ncol=length(motif2show),nrow=length(levels(ORN)))
colnames(motif_matrix)<-motif2show
rownames(motif_matrix)<-levels(ORN)

last_motif_info<-All_motif_info[which(All_motif_info$motif.name%in%motif2show),]
for (i in 1:nrow(last_motif_info)){
  cluster=last_motif_info[i,]$cluster;
  motif=last_motif_info[i,]$motif.name;
  motif_matrix[cluster,motif]<-last_motif_info[i,]$fold.enrichment
}

library(pheatmap)
#count=t(scale(t(motif_matrix),scale = T,center = T))
pdf("./ORN/07_DEG_and_DEP/01_All_ORN_cluster/tau_DEP_motif_heatmap.pdf",width=20,height=20)
pheatmap(motif_matrix,cluster_cols = T,cluster_rows = F,
              color = colorRampPalette(c("white", "firebrick3"))(100),
              cellwidth = 10, cellheight = 10,
              show_rownames=T,show_colnames=T)
dev.off()


#Find markers
All_motif_info <- data.frame()
for (cluster in levels(ORN)){
  cluster_peak <- markers[markers$cluster==cluster,]$gene;
  enriched.motifs <- FindMotifs(ORN,features = cluster_peak);
  enriched.motifs$cluster <- cluster;
  All_motif_info <- rbind(All_motif_info,enriched.motifs)
}
library(dplyr)
top3 <- All_motif_info %>% group_by(cluster) %>% top_n(n = 3, wt = fold.enrichment)
motif2show<-unique(top3$motif.name)

# motif enrich matrix
motif_matrix<-matrix(ncol=length(motif2show),nrow=length(levels(ORN)))
colnames(motif_matrix)<-motif2show
rownames(motif_matrix)<-levels(ORN)

last_motif_info<-All_motif_info[which(All_motif_info$motif.name%in%motif2show),]
for (i in 1:nrow(last_motif_info)){
  cluster=last_motif_info[i,]$cluster;
  motif=last_motif_info[i,]$motif.name;
  motif_matrix[cluster,motif]<-last_motif_info[i,]$fold.enrichment
}

library(pheatmap)
#count=t(scale(t(motif_matrix),scale = T,center = T))
pdf("./ORN/remove_nopower/DEGandDEP/FindAllMarker_DEP_motif_heatmap.pdf",width=20,height=20)
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
ORN <- RunChromVAR(
  object = ORN,
  genome = BSgenome.Amel.antenan
)

p2 <- FeaturePlot(
  object = ORN,
  features = head(rownames(enriched.motifs)),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
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
  genome = BSgenome.Amel.antenan,
  pfm = pfm
)

ORN <- Footprint(
  object = ORN,
  motif.name =head(rownames(differential.activity)),
  genome = BSgenome.Amel.antenan
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
 