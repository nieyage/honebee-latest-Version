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

honeybee<-readRDS("./02_All_celltype/WNN_honeybee_integrated_all_celltype.rds")

# Fig1B: 
# NE:#5CC8F2,Nurse:#009E73,F:#E69F00
pdf("./00_Figure/Fig1B-honeybee_cluster_WNN_top.pdf",width=15,height=5)
###cluster
p1 <- DimPlot(honeybee,reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(honeybee,reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(honeybee,reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off();
pdf("./00_Figure/Fig1B-honeybee_cluster_WNN_bottom.pdf",width=15,height=5)
###sample
p1 <- DimPlot(honeybee, cols=c("#5CC8F2","#009E73","#E69F00"), reduction = "umap.rna", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(honeybee, cols=c("#5CC8F2","#009E73","#E69F00"), reduction = "umap.atac",group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(honeybee, cols=c("#5CC8F2","#009E73","#E69F00"), reduction = "wnn.umap", group.by = "orig.ident", shuffle=TRUE,label = F, label.size = 5, repel = TRUE) + ggtitle("WNNUMAP")
p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))
dev.off()
# Fig1C: 
Idents(honeybee)<-honeybee$Annotation;
pdf("./00_Figure/Fig1C-honeybee_annotation_allcelltype_WNN-Annotation_first.pdf",width=7,height=5)
DimPlot(honeybee, label = T, repel = TRUE, cols=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#B3B3B3"), reduction = "wnn.umap",group.by = "Annotation")+ ggtitle("")
DimPlot(honeybee, label = F, repel = TRUE, cols=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#B3B3B3"), reduction = "wnn.umap",group.by = "Annotation")+ ggtitle("")
dev.off()

# Fig1D:
# Major celltype track and violin plot 
##Track for Marker genes promoters

DefaultAssay(honeybee) <- "peaks"
# first compute the GC content for each peak
honeybee <- RegionStats(honeybee, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(honeybee)$tx_id <-Annotation(honeybee)$gene_name
#features<-c("Or2","LOC411079","LOC410151","LOC406073","LOC409780","Obp5","Obp11","Obp4")
# link peaks to genes

honeybee <- LinkPeaks(
  object = honeybee,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Syt1","LOC411079","LOC410151","5-ht7")
)
######Visulize track and RNA exp######
idents.plot <- Idents(honeybee)
# Neuron:
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group1-5132500-5137500",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group1-5125000-5137500"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "Syt1",
  assay = "RNA"
)
p1<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Epithelial cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group5-6387000-6396000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group5-6387000-6396000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC411079",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

# glial cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group15-2484000-2489000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group15-2484000-2489000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "LOC410151",
  assay = "RNA"
)
p3<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Sheath cell 
cov_plot <- CoveragePlot(
  object = honeybee,
  region = "Group6-730000-740000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = honeybee,
  region = "Group6-730000-740000"
)
expr_plot <- ExpressionPlot(
  object = honeybee,
  features = "5-ht7",
  assay = "RNA"
)
p4<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
set<-c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F","#B3B3B3")
p1<-p1& scale_fill_manual(values=set)&labs(title="Syt1")
p2<-p2& scale_fill_manual(values=set)&labs(title="LOC411079(GRH)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)&labs(title="LOC410151(repo)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p4<-p4& scale_fill_manual(values=set)&labs(title="5-ht7") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
pdf("./00_Figure/Fig1D-Marker_gene-select-peaktrack-WNN.pdf",height=8,width=16) 
p1|p2|p3|p4
dev.off()

# Fig1E:
Neuron<-readRDS("./03_Neuron/WNN_Neuron_integrated.rds")
markers <- FindAllMarkers(Neuron, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Neuron_top<-markers[markers$cluster=="Orco-Neuron",7]
# Or2 and LOC408554
DefaultAssay(Neuron)<-"RNA"
pdf('./00_Figure/Fig1E-Neuron_marker_FeaturePlot_WNN.pdf', width=13, height=4)
p1<- FeaturePlot(Neuron, reduction = 'wnn.umap',max.cutoff = 7,features = c("LOC413063") ,order=TRUE, ncol = 1)+ggtitle("LOC413063 (pepple)")
p2<- FeaturePlot(Neuron,cols =c("lightgrey", "#183A1D"), reduction = 'wnn.umap',max.cutoff = 7,features = c("LOC551837") ,order=FALSE, ncol = 1)+ggtitle("LOC551837 (bgm)")
p3<- FeaturePlot(Neuron,cols =c("lightgrey", "#F0A04B"), reduction = 'wnn.umap',max.cutoff = 25,features = c("Or2") ,order=TRUE, ncol = 1)+ggtitle("Or2 (Orco)")
p1|p2|p3
dev.off()

# Or2 and LOC408554
DefaultAssay(Neuron)<-"RNA"
pdf('./03_Neuron/Orco-Neuron_marker_FeaturePlot_WNN.pdf', width=10, height=10)
 FeaturePlot(Neuron, reduction = 'wnn.umap',max.cutoff = 7,features = Neuron_top[1:9] ,order=FALSE, ncol = 3)
dev.off()
# Fig1F:
DefaultAssay(Neuron) <- "peaks"
# first compute the GC content for each peak
Neuron <- RegionStats(Neuron, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(Neuron)$tx_id <- Annotation(Neuron)$gene_name
Neuron <- LinkPeaks(
  object = Neuron,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("LOC413063","Or2",Neuron_top)
)
######Visulize track and RNA exp######
idents.plot <- Idents(Neuron)

# Neuron
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group5-12940000-12944000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group5-12940000-12944000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "LOC413063",
  assay = "RNA"
)
p1<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# LOC551837
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group1-22596000-22597500",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group1-22596000-22597500"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "LOC551837",
  assay = "RNA"
)
p3<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)
# Orco
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group1-5722000-5725000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group1-5722000-5725000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "Or2",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

# Orco
cov_plot <- CoveragePlot(
  object = Neuron,
  region = "Group1-5722000-5725000",
  #annotation = FALSE,
  peaks = FALSE,links = F,annotation = F
)
gene_plot <- AnnotationPlot(
  object = Neuron,
  region = "Group1-5722000-5725000"
)
expr_plot <- ExpressionPlot(
  object = Neuron,
  features = "Or2",
  assay = "RNA"
)
p2<-CombineTracks(
  plotlist = list(cov_plot,gene_plot),
  expression.plot = expr_plot,
  heights = c(10, 3),
  widths = c(10,5)
)

set<-c("#F0A04B", "#183A1D")
p1<-p1& scale_fill_manual(values=set)&labs(title="LOC413063 (pepple)")
p2<-p2& scale_fill_manual(values=set)&labs(title="Or2 (Orco)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)&labs(title="LOC551837 (bgm)") & theme(strip.text.y.left = element_blank(),strip.background = element_blank())

pdf("./00_Figure/Fig1F-Neuron-Orco-track.pdf",width=13,height=4)
p1|p3|p2
dev.off()


#Fig1G: proportion of cell: ORN vs Non ORN;
ORN<-c(2405,2032,2840)
all<-c(2773,2555,3540)
Non_ORN<- all-ORN;
data<-data.frame(group=c(rep("ORN",3),rep("Non-ORN",3)),
  Sample=rep(c("NE","Nurse","Forager"),2),
  cellnumber=c(ORN,Non_ORN))
data$group<-factor(data$group,levels=c("ORN","Non-ORN"))
data$Sample<-factor(data$Sample,levels=c("NE","Nurse","Forager"))
pdf("./00_Figure/Fig1G-ORNvsNonORN_proportion.pdf",width=4,height=4)
ggplot(data = data, aes_string(x = "group", y = "cellnumber", 
        fill = "Sample")) +  xlab("orig.ident") + ylab("Percent of cells") + 
        scale_fill_manual(values = c("#5CC8F2","#009E73","#E69F00")) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_bw();
ggplot(data = data, aes_string(x = "group", y = "cellnumber", 
        fill = "Sample")) +  xlab("orig.ident") + ylab("cell number") + 
        scale_fill_manual(values = c("#5CC8F2","#009E73","#E69F00")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()
dev.off()


#cross species
#Fig1H:
#OR vs Non-ORN porportion 
# change to 3 types
# Orco+Neuron+;Orco-Neuron+;Orco-Neuron-;
Apis_mellifera<-c(7277,55,1536);
#Fly_Joint.integrated <- readRDS("/md01/liyh526/project/Fly-cooperate/shiny_app/WNN_fly_integrated_annotation_antenna.rds")
#Fly_ORN<-readRDS("/md01/liyh526/project/Fly-cooperate/7.18run/WNN_ORN_integrated_antenna.rds")
#mosquito_all<-readRDS("/data/R02/nieyg/project/honeybee/data/publish_data/mosquito/SeuratObject1_Antenna_mergedBatches_AllCells.rds")
#mosquito_neuron<-readRDS("/data/R02/nieyg/project/honeybee/data/publish_data/mosquito/SeuratObject2_Antenna_mergedBatches_Neurons.rds")
#Drosophila_melanogaster<-c(2808,154,3756) #our data 
Drosophila_melanogaster<-c(7940,4266,25048)
Aedes_aegypti<-c(4742,433,8704);
cross_species_cellnumber<-data.frame(species=c(rep("Apis mellifera",3),rep("D.melanogaster",3),rep("Ae.aegypti",3)),
  celltype=rep(c("Orco+Neuron","Orco-Neuron","Non-neuron"),3),
  cellnumber=c(Apis_mellifera,Drosophila_melanogaster,Aedes_aegypti))
cross_species_cellnumber$celltype<-factor(cross_species_cellnumber$celltype,levels=c("Orco+Neuron","Orco-Neuron","Non-neuron"));
cross_species_cellnumber$species<-factor(cross_species_cellnumber$species,levels=c("Apis mellifera","D.melanogaster","Ae.aegypti"))
pdf("./00_Figure/Fig1H-cross_species_ORNvsNonORN_proportion_3types_publishdata.pdf",width=4,height=4)
ggplot(data = cross_species_cellnumber, aes_string(x = "species", y = "cellnumber", 
        fill = "celltype")) +  xlab(" ") + ylab("% Percent of cells") + 
        scale_fill_manual(values = c("#F0A04B","#183A1D","#D8D8D8")) + 
        geom_bar(position = "fill", stat = "identity", width = 0.6) +
        theme_bw()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
dev.off();

# Fig1I:
# OR,IR,GR gene number
chemoreceptor <- read.table("/data/R02/nieyg/ref/10X/Amel_HAv3.1/chemoreceptor_info_data.csv",header=TRUE)
OR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="OR",]$gene_name)
GR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="GR",]$gene_name)
IR_gene<- unique(chemoreceptor[chemoreceptor$gene_type=="IR",]$gene_name)
Apis_mellifera<-c(length(OR_gene),length(IR_gene),length(GR_gene));
Drosophila_melanogaster<-c(60,66,60)
Aedes_aegypti<-c(114,135,107)
cross_species_genenumber<-data.frame(species=c(rep("Apis_mellifera",3),rep("Drosophila_melanogaster",3),rep("Aedes_aegypti",3)),
  genetype=rep(c("ORs","IRs","GRs"),3),
  genenumber=c(Apis_mellifera,Drosophila_melanogaster,Aedes_aegypti))
cross_species_genenumber$genetype<-factor(cross_species_genenumber$genetype,levels=c("ORs","IRs","GRs"));
cross_species_genenumber$species<- factor(cross_species_genenumber$species,levels=c("Apis_mellifera","Drosophila_melanogaster","Aedes_aegypti"))
pdf("./00_Figure/Fig1I-cross_species_ORGRIRgene_proportion.pdf",width=4,height=4)
p<-ggplot(data = cross_species_genenumber, aes_string(x = "species", y = "genenumber", 
        fill = "genetype")) +  xlab(" ") + ylab("chemosensory receptor genes") + 
        scale_fill_manual(values = c("#9D6EB0","#EA6376","#3D63AD")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = genenumber), size = 3, hjust = 0.5, vjust = 3, position = "stack") 
dev.off();

#Fig1J:
# the OB barplot
cross_species_glomeruli<-data.frame(species=c("Apis mellifera","D. melanogaster","Ae. aegypti"),
  glomeruli=c(160,55,65));
cross_species_glomeruli$species<- factor(cross_species_glomeruli$species,levels=c("Apis mellifera","D. melanogaster","Ae. aegypti"))
pdf("./00_Figure/Fig1J-cross_species_glomeruli.pdf",width=4,height=4)
p<-ggplot(data = cross_species_glomeruli, aes_string(x = "species", y = "glomeruli", 
        fill = "species")) +  xlab(" ") + ylab("# of glomeruli") + 
        scale_fill_manual(values = c("#C8B6E2" ,"#7A86B6" ,"#495C83")) + 
        geom_bar( stat = "identity", width = 0.6) +
        theme_classic()+
        theme(axis.text.x = element_text(angle = 25,vjust = 0.5,hjust = 0.5));
p
#add gene number in plot 
p+geom_text(aes(label = glomeruli), size = 3, hjust = 0.5, vjust = 3) 
dev.off();
