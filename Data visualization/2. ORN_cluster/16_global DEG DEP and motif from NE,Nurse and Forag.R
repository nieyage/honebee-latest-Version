# global DEG DEP and motif from NE,Nurse and Forager
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
set.seed(1234)
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(c("LOC412949","LOC100577496","LOC102653640","LOC727346","LOC100578352","LOC552552","LOC726019","LOC551704","LOC410623","LOC100576097","LOC409777"))
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
all_receptor_gene <- unique(c(Orco,OR_gene,IR_gene,GR_gene))
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
DefaultAssay(ORN)<-"raw_RNA"
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")

ORN_correct<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_SCT_v2.rds")
# global Stage DEG 
obj<- ORN_correct
DefaultAssay(obj)<- "SCT_v2";
obj$orig.ident<- factor(obj$orig.ident,levels=c("NE","Nurse","Forager"))
Idents(obj)<- obj$orig.ident
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
markers <- markers[which(markers$p_val<0.05),]
write.csv(markers,paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/global_specific_markers_FC0.5.csv",sep=""))
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
#obj<-ScaleData(obj,features=rownames(obj))
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
p<- Seurat::DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =c("#E69F00","#55B4E9","#009E73"),disp.min = -2,disp.max = 2,size = 6,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#F7F8F2", "#DBE9D0", "#99D3BC","#50B1D6","#03418B"))
pdf(paste0("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/global_specific_markers_heatmap_FC0.5.pdf",sep=""),width=8,height=16)
print(p)
dev.off()
DEG<- markers
#  # R MAGIC processing data and plot heatmap
#  library(reticulate)
#  use_python("/public/home/nieyg/biosoft/conda/envs/r4-base/bin/python")
#  import("magic")
#  library(Rmagic)
#  library(ggplot2)
#  magic_data<- as.matrix(ORN_correct@assays$SCT_v2@counts[markers$gene,])
#  MAGIC_data <- magic(magic_data, genes="all_genes",t=4)
#  data<- as.data.frame(MAGIC_data)
#  
#  library(ComplexHeatmap)
#  ORN_correct$orig.ident<- factor(ORN_correct$orig.ident,levels=c("NE","Nurse","Forager"))
#  mat <- magic_data[markers$gene,]
#  gene_features <- markers$gene
#  cluster_info <- sort(ORN_correct$orig.ident)
#  mat <- as.matrix(mat[,names(cluster_info)])
#  ##输入想在图上展示出来的marker基因，获得基因在热图中的位置信息
#  gene <- top10$gene
#  gene_pos <- which(rownames(mat)%in%gene)
#  row_anno<-rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene))
#  col <- c("#E69F00","#55B4E9","#009E73")
#  names(col) <- levels(cluster_info)
#  top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=col),
#    labels=levels(cluster_info),labels_gp=gpar(cex=0.5,col='white')))
#  
##给热图改个好看的颜色
#  library(circlize)
#  col_fun = colorRamp2(c(-2,-1,0,1,2),c("#F7F8F2", "#DBE9D0", "#99D3BC","#50B1D6","#03418B"))
#  pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/MAGIC_global_specific_markers_heatmap_FC0.5.pdf",width=8,height=16)
#  Heatmap(mat,
#  cluster_rows = FALSE,
#  cluster_columns = FALSE,
#  show_column_names = FALSE,
#  show_row_names = FALSE,
#  column_split = cluster_info,
#  top_annotation = top_anno,
#  column_title = NULL,
#  right_annotation = row_anno,
#  heatmap_legend_param = list(
#  title='Expression',
#  title_position='leftcenter-rot'),
#  col = col_fun)
#  dev.off()



library(AnnotationHub)
library(biomaRt)
library(dplyr)
library(goseq)
library(DOSE)
library(stringr)
library(clusterProfiler)
library(GenomicRanges)
library(AnnotationDbi)
Apis_mellifera.OrgDb <-loadDb("/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/Apis_mellifera_AH102515.OrgDb")
columns(Apis_mellifera.OrgDb)

  gene_in_cluster<-markers[markers$cluster=="NE",]$gene
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  NE_ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
  gene_in_cluster<-markers[markers$cluster=="Nurse",]$gene
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  Nurse_ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
  gene_in_cluster<-markers[markers$cluster=="Forager",]$gene
  gene.df <- bitr(gene_in_cluster, fromType = "SYMBOL",toType =  "ENTREZID",OrgDb = Apis_mellifera.OrgDb)
  Forager_ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = Apis_mellifera.OrgDb,
                ont = type,
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE);
NE_ego<- as.data.frame(NE_ego)
NE_ego$stage<-"NE"
Nurse_ego<- as.data.frame(Nurse_ego)
Nurse_ego$stage<-"Nurse"
Forager_ego<- as.data.frame(Forager_ego)
Forager_ego$stage<-"Forager"

all_ego<-rbind(NE_ego,Nurse_ego,Forager_ego)
write.csv(all_ego,"./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Global_Stage_specific_markers_GO.csv")

top5_term <- all_ego %>% group_by(stage) %>% top_n(n = 8, wt = pvalue);
top5_term$stage<-factor(top5_term$stage,levels=c("NE","Nurse","Forager"))
library(ggplot2)
pdf("./05_ORN_cluster2/06_stage_specific/UMI_corrected/cluster_DEG/Global_Stage_specific_markers_GO.pdf",width=16,height=10)
p <- ggplot(top5_term,aes(y=Count,x=Description,fill=pvalue)) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(stage~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 14),
            legend.position="right",
            legend.title = element_text(size=18),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
p
dev.off()

# recall peak 
# ORN recall peak for subcluster 
DefaultAssay(ORN)<-"ATAC"
peak<-CallPeaks(
       ORN,
       group.by = "orig.ident",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="./05_ORN_cluster2/08_call_peak_by_stage/",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(ORN),
     features = peak,
     cells = colnames(ORN)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# create a new assay using the MACS2 peak set and add it to the Seurat object
ORN[["peaks_ORN_stage"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(ORN),
  annotation = Annotation(ORN)
)
DefaultAssay(ORN) <- "raw_RNA"
saveRDS(ORN,"./05_ORN_cluster2/08_call_peak_by_stage/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_by_stage.rds")

DefaultAssay(ORN) <- 'peaks_ORN_stage'
ORN$orig.ident<- factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
Idents(ORN)<- ORN$orig.ident
markers <- FindAllMarkers(ORN,test.use = 'LR', only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.1,latent.vars = 'nCount_peaks')
markers <- markers[which(markers$p_val<0.05),]
ORN<- ScaleData(ORN)
p<- DoHeatmap(object = ORN,features=markers$gene,disp.min = -1,disp.max = 1.5,label=T,group.colors =c("#E69F00","#55B4E9","#009E73"),size = 6,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#FDF7F7", "#FDD7D5", "#F06FA2","#CA007F","#70176B"))
pdf(paste0("./05_ORN_cluster2/06_stage_specific/global_specific_DEP_heatmap.pdf",sep=""),width=8,height=16)
print(p)
dev.off()
DEP<- markers
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'insects', all_versions = FALSE)
)
pfm_honeybee <- readJASPARMatrix("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/CisBP-honeybee.jaspar", matrixClass=c("PFM", "PWM", "PWMProb"))
# Merge the two PFMatrixList objects
mergedPfmList <- c(pfm, pfm_honeybee)
DefaultAssay(ORN) <- 'peaks_ORN_stage'
ORN$orig.ident<- factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
Idents(ORN)<- ORN$orig.ident

# add motif information
ORN <- AddMotifs(
  object = ORN,
  genome = BSgenome.Amel.HAv3.1.update.chemoreceptor,
  pfm = mergedPfmList
)
ORN <- RegionStats(ORN, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)

obj <- RunChromVAR(object = ORN,genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)

DefaultAssay(obj) <- 'chromvar'
differential.activity <- FindAllMarkers(object = obj,only.pos = TRUE,mean.fxn = rowMeans,fc.name = "avg_diff")
markers <- differential.activity[which(differential.activity$p_val<0.05),]
markers <- markers[which(markers$avg_diff>0.3),]

whitePurple = c("9"='#f7fcfd',"6"='#e0ecf4',"8"='#bfd3e6',"5"='#9ebcda',"2"='#8c96c6',"4"='#8c6bb1',"7"='#88419d',"3"='#810f7c',"1"='#4d004b')
obj<- ScaleData(obj)
p<- DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =c("#E69F00","#55B4E9","#009E73"),disp.min = -2,disp.max = 2,size = 6,group.by = "orig.ident") + scale_fill_gradientn(colors =whitePurple)
pdf(paste0("./05_ORN_cluster2/06_stage_specific/global_specific_chromVAR_heatmap.pdf",sep=""),width=8,height=16)
print(p)
dev.off()

motif<- unique(gsub("-","_",markers$gene))

motif.obj <- SeuratObject::GetAssayData(ORN, slot = "motifs")
MotifPlot(ORN, motifs = motif)


pdf("./05_ORN_cluster2/06_stage_specific/MotifPlot2.pdf",width=6,height=25)
MotifPlot(
  object = obj,
  use.names = FALSE,
  motifs = motif[-19],
  ncol=3,
  assay = 'peaks_ORN_stage'
)
dev.off()

#data.use <- GetMotifData(object = ORN, assay = 'peaks_ORN_stage', slot = "pwm")
#motifs<- ConvertMotifID(object = ORN, name = motif[19])
#data.use <- data.use[motifs]
#names(x = data.use) <- GetMotifData(
#      object = object, assay = assay, slot = "motif.names"
#    )[motifs]
#p <- ggseqlogo::ggseqlogo(data = data.use)



motif<- gsub(" .*","",motif)
# enriched motif Expression heatmap
motif2TF<- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_promoter_TOBIAS/motif2TF.txt")
# get the TF info 
TF<- motif2TF[match(motif,motif2TF$V1),]$V2
fly2honeybee<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/10_fly_honeybee_ID_trans/05_fly2honeybee.csv")
TF_name<- TF

for (i in which(TF%in%fly2honeybee$fly_gene)){
  TF_name[i]<- fly2honeybee[match(TF[i],fly2honeybee$fly_gene),]$honeybee_gene_name
}
obj<- ORN_correct
DefaultAssay(obj)<- "SCT_v2";
obj<- ScaleData(obj,features=TF_name)
p<- DoHeatmap(object = obj,features=TF_name,label=T, group.colors =c("#E69F00","#55B4E9","#009E73"),disp.min = -2,disp.max = 2,size = 6,group.by = "orig.ident") + scale_fill_gradientn(colors = c("#F7F8F2", "#DBE9D0", "#99D3BC","#50B1D6","#03418B"))
pdf(paste0("./05_ORN_cluster2/06_stage_specific/TF_global_specific_markers_heatmap.pdf",sep=""),width=8,height=16)
print(p)
dev.off()

# example track (peak link to gene)
# Fig3H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
ORN_correct <- Seurat::SCTransform(ORN,
  assay="RNA",
  new.assay.name = "SCT_v2",
  ncells = 4000,
    residual.features = NULL,# all gene
    variable.features.n = 3000,
    do.scale = FALSE,
    do.center = TRUE,
    return.only.var.genes = FALSE,
    seed.use = 1448145,
    verbose = TRUE,
  vst.flavor="v2")

# filter the DEP
filter_DEP<- DEP$gene
peak_assay<- GetAssayData(ORN, assay = "peaks_ORN_stage")
peak_data<-peak_assay[filter_DEP,]
#raw_assay <- CreateAssayObject(counts = peak_data)
# add this assay to the previously created Seurat object
DefaultAssay(ORN_correct)<-"peaks_ORN_stage"

ORN_correct[["peaks_ORN_stage_filtered"]] <- CreateChromatinAssay(
  counts = peak_data,
  #fragments = fragpath[[i]],
  fragments = Fragments(ORN_correct),
  annotation = Annotation(ORN_correct)
)

DefaultAssay(ORN_correct)<-"peaks_ORN_stage_filtered"
Idents(ORN_correct)<- ORN_correct$orig.ident
obj <- RegionStats(ORN_correct, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "peaks_ORN_stage_filtered",
  expression.assay = "SCT_v2",
  genes.use = DEG$gene
)
######Visulize track and RNA exp######
idents.plot <- Idents(obj)
gene <- DEG$gene
annotations <- Annotation(object = obj[["peaks_ORN_stage_filtered"]])

gene<- gene[which(gene%in% annotations$gene_name)]
gene_data<- as.data.frame(annotations[match(gene,annotations$gene_name),])
gene_data<- gene_data[which(gene_data$seqnames!="MT"),]
gene<- gene_data$gene_name;

# how to set the track plot color ?

pdf("./05_ORN_cluster2/08_call_peak_by_stage/DEG_DEP_trackplot.pdf",width=10,height=6)
for (i in gene){
  print(i)
  p1 <- CoveragePlot(
  object = obj,
  region = i,
  features = i,
  expression.assay = "SCT_v2",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000,
  links=TRUE
)
  p2<- p1 + scale_fill_manual(values = c("#E69F00","#55B4E9","#009E73"))
print(p2)
}
dev.off()


setwd("/data/R04/liyh526/project/Honeybee/05_honeybee_bulkATAC_macs2_callpeak/2023.5.15/bulk_sc_associate/Final_figure")
dir()
library(ChIPseeker)
LOC552079.gr <- readPeakFile("LOC552079_motifs.bed", as = "GRanges")
LOC411133.gr <- readPeakFile("LOC411133_motifs.bed", as = "GRanges")
LOC725220.gr <- readPeakFile("LOC725220_motifs.bed", as = "GRanges")
LOC552079_pos_in_peak.gr <- LOC552079.gr[unique(queryHits(findOverlaps(LOC552079.gr,peak.gr))),]
LOC552079_pos_in_peak.gr <- resize(LOC552079_pos_in_peak.gr,width=25,fix="center")
LOC411133_pos_in_peak.gr <- LOC411133.gr[unique(queryHits(findOverlaps(LOC411133.gr,peak.gr))),]
LOC411133_pos_in_peak.gr <- resize(LOC411133_pos_in_peak.gr,width=20,fix="center")
LOC725220_pos_in_peak.gr <- LOC725220.gr[unique(queryHits(findOverlaps(LOC725220.gr,peak.gr))),]
LOC725220_pos_in_peak.gr <- resize(LOC725220_pos_in_peak.gr,width=20,fix="center")
# 把motif的数量按bingscore筛选一下，只保留bindscore大于均值的
LOC552079_pos_in_peak.gr <- LOC552079_pos_in_peak.gr[LOC552079_pos_in_peak.gr$V5>mean(LOC552079_pos_in_peak.gr$V5)]
LOC411133_pos_in_peak.gr <- LOC411133_pos_in_peak.gr[LOC411133_pos_in_peak.gr$V5>mean(LOC411133_pos_in_peak.gr$V5)]
LOC725220_pos_in_peak.gr <- LOC725220_pos_in_peak.gr[LOC725220_pos_in_peak.gr$V5>(mean(LOC725220_pos_in_peak.gr$V5-0.6))]
LOC552079_pos_plot <- ggplot(data=as.data.frame(LOC552079_pos_in_peak.gr)) +
  geom_segment(aes(x=start, y = 0,xend = end, yend = 0),size = 3,color="#D51B07")+
  theme_classic() +
  ylab(label = "LOC552079") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),panel.background=element_rect(fill='transparent', color='black',linetype="solid")) +
  xlim(c(934699, 948738)) +
  geom_rect(aes(xmin=935573,xmax=935782,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2) +
  geom_rect(aes(xmin=943433,xmax=943639,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2)


LOC552079_pos_plot <- ggplot(data=as.data.frame(LOC552079_pos_in_peak.gr)) +
  geom_segment(aes(x=start, y = 0,xend = end, yend = 0),size = 3,color="#D51B07")+
  theme_classic() +
  ylab(label = "LOC552079") +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),plot.margin=unit(c(0,0,0,0), "cm"),panel.background=element_rect(fill='transparent', color='black',linetype="solid")) +
  xlim(c(934699, 948738)) +
  geom_rect(aes(xmin=935573,xmax=935782,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2) +
  geom_rect(aes(xmin=943433,xmax=943639,ymin=-Inf,ymax=Inf),
            fill="gray",alpha=0.2)









library(cicero)
# convert to CellDataSet format and make the cicero object
obj.cds <- as.cell_data_set(x = obj)
obj.cicero <- make_cicero_cds(obj.cds, reduced_coordinates = reducedDims(obj.cds)$UMAP)
# get the chromosome sizes from the Seurat object
genome <- seqlengths(bone)

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(bone.cicero, genomic_coords = genome.df, sample_num = 100)

ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(bone) <- links




heatmap_output<- p
rows <- heatmap_output$data$Feature
columns <- heatmap_output$data$Cell

# 平滑处理热图数据
data_matrix <- magic_testdata
smoothed_data_matrix <- matrix(NA, nrow(data_matrix), ncol(data_matrix))
window_size <- 3  # 平滑窗口大小，可以根据需要进行调整

for (i in 1:nrow(data_matrix)) {
  for (j in 1:ncol(data_matrix)) {
    row_start <- max(1, i - window_size)
    row_end <- min(nrow(data_matrix), i + window_size)
    col_start <- max(1, j - window_size)
    col_end <- min(ncol(data_matrix), j + window_size)
    
    smoothed_data_matrix[i, j] <- mean(data_matrix[row_start:row_end, col_start:col_end])
  }
}

processed_assay <- CreateAssayObject(counts = smoothed_data_matrix)

# 将新的数据集添加到Seurat对象的assays属性中
ORN_correct$assays$processed_assay <- processed_assay








