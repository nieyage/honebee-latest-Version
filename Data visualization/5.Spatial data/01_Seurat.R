## create Seurat by gex matrix 
library(Seurat)
library(cowplot)
library(patchwork)
library(sctransform)
library(celda)
library(DropletUtils)
library(Matrix)
library(DoubletFinder)
set.seed(1234)
bin20_obj<- readRDS("/md01/nieyg/project/honeybee/data/Spatial_data/bin20.seurat.rds")
bin50_obj<- readRDS("/md01/nieyg/project/honeybee/data/Spatial_data/bin50.seurat.rds")
bin20_obj@meta.data[1:4,]

ORN<- bin20_obj

# Or2 feature plot in global umap 
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
Idents(ORN)<-ORN$seurat_clusters;
pdf("./15_spatial_data/ORN_global_tsne_umap_plot.pdf",width=7,height=5)
DimPlot(ORN, label = T, repel = TRUE, cols=myUmapcolors, reduction = "tsne",group.by = "seurat_clusters")+ ggtitle("")
DimPlot(ORN, label = F, repel = TRUE, cols=myUmapcolors, reduction = "umap",group.by = "seurat_clusters")+ ggtitle("")
dev.off()

pdf('./15_spatial_data/Or2_Or25_26_27_marker_FeaturePlot_WNN.pdf', width=18, height=4)
FeaturePlot(ORN,reduction = 'tsne',max.cutoff = 10,features = c("Or2","Or25","Or26","Or27"),order=TRUE, ncol = 4)
dev.off()


DefaultAssay(ORN)<-"SCT"
gene<-c("Or2","Or25","Or26","Or27")
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
ORN_matrix<- as.data.frame(ORN_matrix)
coexp_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or27>0&ORN_matrix$Or26>0&ORN_matrix$Or25>0),])



obj<- ORN
DefaultAssay(obj)<-"SCT"
gene<- c("Or2","LOC102655285","LOC107965761")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
colnames(ORN_matrix)<- c("Or2","Or154","Or163")
ORN_matrix<- as.data.frame(ORN_matrix)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or163,ORN_matrix$Or154,decreasing=T),]
Or154_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or154>0),])
Or163_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or163>0),])
coexp_barcode<- intersect(Or154_barcode,Or163_barcode)

list_name<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")

obj<- ORN
DefaultAssay(obj)<-"SCT"
gene<- c("Or2","LOC102655285","LOC107965761")
ORN_count<-obj@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%gene),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
ORN_matrix<- t(ORN_matrix)
colnames(ORN_matrix)<- c("Or2","Or154","Or163")
ORN_matrix<- as.data.frame(ORN_matrix)
ORN_matrix<- ORN_matrix[order(ORN_matrix$Or163,ORN_matrix$Or154,decreasing=T),]
Or154_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or154>0),])
Or163_barcode<- rownames(ORN_matrix[which(ORN_matrix$Or163>0),])
coexp_barcode<- intersect(Or154_barcode,Or163_barcode)




data<- GetAssayData(ORN)
exp_sum<- rowSums(data)[order(rowSums(data))]





/md01/nieyg/project/honeybee/data/Spatial_data/Rawdata/RawData/tentacle/tentacle.fq.gz

# fastp 
nohup fastp -i /md01/nieyg/project/honeybee/data/Spatial_data/Rawdata/RawData/tentacle/tentacle.fq.gz -o tentacle_trimmed.fq.gz &


# SE STAR mapping 

nohup STAR --runThreadN 12 \
        --genomeDir /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star \
        --readFilesIn tentacle_trimmed.fq.gz  \
        --readFilesCommand zcat \
        --outFileNamePrefix tentacle2 \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 25074000000 & 

featureCounts -T 1 -p \
-t exon -g gene_id \
-a /md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz -o tentacle.txt tentacle2Aligned.sortedByCoord.out.bam
cat tentacle.txt | cut -f1,7- > tentacle_counts.txt


# bam 2 bw files 
  bedtools  genomecov  -bg -split -ibam tentacle2Aligned.sortedByCoord.out.bam  > tentacle.bedGraph
  perl /data/R02/nieyg/pipeline/ATACseq/norm_bedGraph.pl tentacle.bedGraph tentacle.norm.bedGraph &> tentacle.norm.bedGraph.log
  sort -k1,1 -k2,2n tentacle.norm.bedGraph > tentacle.norm.sorted.bedGraph
  bedGraphToBigWig tentacle.norm.sorted.bedGraph ~/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt tentacle.norm.bw 



# Four Orco Tsne plot 
Orco<- c("Or2","LOC552552","LOC726019","LOC551704")
pdf('./15_spatial_data/Orco_marker_FeaturePlot_WNN.pdf', width=12, height=10)
FeaturePlot(ORN,reduction = 'tsne',max.cutoff = 10,features = Orco,order=TRUE, ncol = 2)
dev.off()

# upset R 
library(UpSetR)
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%Orco),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or2 = names(which(ORN_matrix[4,]>0)), 
        LOC552552 = names(which(ORN_matrix[2,]>0)), 
        LOC726019 = names(which(ORN_matrix[3,]>0)), 
        LOC551704 = names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./15_spatial_data/Orco_upsetR.pdf", width=8, height=4)
Orco<- c("Or63_b","LOC410603","LOC107963999","LOC100578045")
upset( data)
dev.off()


# Or25_27
features<- c("Or2","Or25","Or26","Or27")
pdf('./15_spatial_data/Or25_27_marker_FeaturePlot_WNN.pdf', width=12, height=10)
FeaturePlot(ORN,reduction = 'tsne',max.cutoff = 10,features = features,order=TRUE, ncol = 2)
dev.off()
# upset R 
library(UpSetR)
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or2 = names(which(ORN_matrix[1,]>0)), 
        Or25 = names(which(ORN_matrix[2,]>0)), 
        Or26 = names(which(ORN_matrix[3,]>0)), 
        Or27 = names(which(ORN_matrix[4,]>0)))
data<- fromList(listInput)
pdf("./15_spatial_data/Or25_27_upsetR.pdf", width=8, height=4)
upset( data)
dev.off()


# Or128_129
features<- c("Or2","LOC102656904","LOC102656221")
pdf('./15_spatial_data/Or128_129_marker_FeaturePlot_WNN.pdf', width=12, height=10)
FeaturePlot(ORN,reduction = 'tsne',max.cutoff = 10,features = features,order=TRUE, ncol = 2)
dev.off()
# upset R 
library(UpSetR)
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
listInput <- list(
        Or2 = names(which(ORN_matrix[3,]>0)), 
        Or128 = names(which(ORN_matrix[2,]>0)), 
        Or129= names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./15_spatial_data/Or128_129_upsetR.pdf", width=8, height=4)
upset( data)
dev.off()

features<- c("Or2","LOC100577101","LOC100577068")
library(UpSetR)
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix[,1:4]
listInput <- list(
        Or2 = names(which(ORN_matrix[3,]>0)), 
        Or39 = names(which(ORN_matrix[2,]>0)), 
        Or40= names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./15_spatial_data/Or39_40_upsetR.pdf", width=8, height=4)
upset( data)
dev.off()

features<- c("Or2","LOC102655285","LOC107965761")
library(UpSetR)
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix[,1:4]
listInput <- list(
        Or2 = names(which(ORN_matrix[3,]>0)), 
        Or154 = names(which(ORN_matrix[1,]>0)), 
        Or163= names(which(ORN_matrix[2,]>0)))
data<- fromList(listInput)
pdf("./15_spatial_data/Or154_163_upsetR.pdf", width=8, height=4)
upset( data)
dev.off()

features<-  c("Or2","Or63-b","LOC410603","LOC107963999","LOC100578045")

library(UpSetR)
ORN_count<-ORN@assays$SCT
ORN_count<-ORN_count[which(rownames(ORN_count)%in%features),]
ORN_matrix<-as.matrix(ORN_count)
ORN_matrix[,1:4]
listInput <- list(
        Or2 = names(which(ORN_matrix[4,]>0)), 
        Or63 = names(which(ORN_matrix[5,]>0)), 
        Or64 = names(which(ORN_matrix[3,]>0)), 
        Or66 = names(which(ORN_matrix[2,]>0)), 
        Or67= names(which(ORN_matrix[1,]>0)))
data<- fromList(listInput)
pdf("./15_spatial_data/Or63_67_upsetR.pdf", width=8, height=4)
upset( data)
dev.off()




Idents(ORN)<- ORN$orig.ident
pdf("./15_spatial_data/samplel_qc_plot_nopoint.pdf")
VlnPlot(ORN, features = c("nCount_Spatial", "nFeature_Spatial", "nCount_SCT","nFeature_SCT"), ncol = 4, pt.size = 0)
dev.off()

ORN_count<-ORN@assays$SCT
gene_total_UMI<- rowSums(ORN_count)
gene_total_UMI<- gene_total_UMI[order(gene_total_UMI,decreasing=TRUE)]
data<- data.frame(gene_total_UMI=gene_total_UMI,gene=names(gene_total_UMI))

data$gene<- factor(data$gene,levels=data$gene)

pdf("./15_spatial_data/seurat_gene_total_UMI.pdf",width=10,height=5)
ggplot(data[1:50,], aes(gene, gene_total_UMI))+geom_point()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
dev.off()


bulk<- read.table("/data/R02/nieyg/project/honeybee/data/Spatial_data/tentacle_counts.txt",header=TRUE)
bulk<-bulk[order(bulk$tentacle,decreasing=TRUE),]
bulk$Geneid<- factor(bulk$Geneid,levels=bulk$Geneid)
pdf("./15_spatial_data/bulk_gene_total_UMI.pdf",width=10,height=5)
ggplot(bulk[1:50,], aes(Geneid, tentacle))+geom_point()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
dev.off()



