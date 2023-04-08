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
ORN<- readRDS("./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_orderbytree.rds")
DefaultAssay(ORN)<-"raw_RNA"

# FigS3A

######################################################################################
# VENN DIAGRAM RAW COUNTS
######################################################################################
receptorSCT.data = ORN@assays$RNA[dotplot_feature,]
OrcoL<- c("Or2","LOC552552","LOC726019","LOC551704")
for (Orco in OrcoL) {
  print(Orco)
  positiveL <- names( receptorSCT.data[Orco,][receptorSCT.data[Orco,] >= 1] )
  positiveL2 <- names( receptorSCT.data[Orco,][receptorSCT.data[Orco,] >= 2] )
  ORN@meta.data[paste0(Orco, '_norExp1')] <- 
    as.numeric(
      lapply(rownames(ORN@meta.data), function(x){
        # print(x)
        if (as.character(x) %in% positiveL) {
          return(1)
        } else {return(0)}
      })
    )
  
  if (Orco == 'Or2') {
    ORN@meta.data[paste0(Orco, '_norExp2')] <- 
      as.numeric(
        lapply(rownames(ORN@meta.data), function(x){
          # print(x)
          if (as.character(x) %in% positiveL2) {
            return(1)
          } else {return(0)}
        })
      )
  }
}
library(VennDiagram)
dataset1 <-  row.names(ORN@meta.data[ORN@meta.data$Or2_norExp1 == 1,])
dataset2 <-  row.names(ORN@meta.data[ORN@meta.data$LOC552552_norExp1 == 1,])
dataset3 <-  row.names(ORN@meta.data[ORN@meta.data$LOC726019_norExp1 == 1,])
dataset4 <-  row.names(ORN@meta.data[ORN@meta.data$LOC551704_norExp1 == 1,])
name1 <- 'Or2 , norm.exp > 1'
name2 <- 'LOC552552 , norm.exp > 1'
name3 <- 'LOC726019 , norm.exp > 1'
name4 <- 'LOC551704 , norm.exp > 1'
vennplot<- venn.diagram(
  x = list(dataset1, dataset2, 
           dataset3, dataset4),
  category.names = c(name1, name2, name3, name4),
  fill = c("#7876B1","#66E1E6","#ffcbcb","#0F3057"),
  filename = NULL,
  output=TRUE
)
pdf("./00_Figure/FigS3A-Orco_Venn.pdf")
grid.draw(vennplot)
dev.off()

# FigS3B Orco Violin plot 
DefaultAssay(ORN) <- "RNA"
Idents(ORN)<-ORN$subcluster
pdf('./00_Figure/FigS3B-Orcocoreceptor_VlnPlot_RNA.pdf',width=25, height=10)
print( VlnPlot(ORN, features = OrcoL, ncol = 1, pt.size = 0.1) )
dev.off()

# FigS3C  SCATTER PLOTS (Or2 Vs potential Orco)

DefaultAssay(ORN) <- "RNA"
scatterColors <- c('#A06CB4', '#DF6C78', '#911A2E', '#CD9139', '#B4B4B6', '#21918c', '#3b528b', '#440154')

### Plotting scatter plot: 1. single-cell level, 2. cluster level

ORN$Or2_UMIs <- ORN@assays$SCT@counts['Or2',]
ORN$LOC552552_UMIs <- ORN@assays$SCT@counts['LOC552552',]
ORN$LOC726019_UMIs <- ORN@assays$SCT@counts['LOC726019',]
ORN$LOC551704_UMIs <- ORN@assays$SCT@counts['LOC551704',]

ORN$Or2_Exp <- ORN@assays$RNA@data['Or2',]
ORN$LOC552552_Exp <- ORN@assays$RNA@data['LOC552552',]
ORN$LOC726019_Exp <- ORN@assays$RNA@data['LOC726019',]
ORN$LOC551704_Exp <- ORN@assays$RNA@data['LOC551704',]

#.....................................................................................
#  Or2 vs. LOC552552 UMI
# ....................................................................................
# 1. single-cell level
pmain <- ORN@meta.data %>%
  ggplot( aes(Or2_UMIs, LOC552552_UMIs) ) + 
  geom_point(size=0.3) + 
  scale_x_log10() + scale_y_log10()
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = ORN@meta.data, aes(x = Or2_UMIs), fill=scatterColors[1]) +scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ORN@meta.data, aes(x = LOC552552_UMIs), fill=scatterColors[2]) + scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.3, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.3, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS3C-1-neuronFigures_Or2-v-LOC552552_UMI.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height =6000,p3)

pmain <- ORN@meta.data %>%
  ggplot( aes(Or2_Exp, LOC552552_Exp) ) + 
  geom_point(size=0.3) +
  geom_hline(yintercept=1, color=scatterColors[3])+
  geom_vline(xintercept = 2, color=scatterColors[3])
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = ORN@meta.data, aes(x = Or2_Exp), fill=scatterColors[1]) +
  geom_vline(xintercept = 2, color=scatterColors[3])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = ORN@meta.data, aes(x = LOC552552_Exp), fill=scatterColors[2]) + 
  geom_vline(xintercept = 1, color=scatterColors[3])+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS3C-1-neuronFigures_Or2-v-LOC552552_Norm.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 6000, height =6000,p3)

library(tidyr)
library(ggrepel)
# 2. in cluster level
dotplot_data<-p$data;

plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC552552')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC552552, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC552552')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC552552), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS3C-1neuronFigures_Or2-v-LOC552552_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC726019')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 
pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC726019, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC726019')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC726019), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS3C-2-neuronFigures_Or2-v-LOC726019_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)


plot_avgExp.clusters.df <- 
  dotplot_data %>%
  filter(features.plot %in% c('Or2', 'LOC551704')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  plot_avgExp.clusters.df %>%
  ggplot(aes(Or2, LOC551704, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Or2') + ylab('Average expression of LOC551704')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = plot_avgExp.clusters.df, aes(x = Or2), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = plot_avgExp.clusters.df, aes(x = LOC551704), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
p3 <- ggdraw(p2)
p3

filename <- paste0("./00_Figure/FigS3C-3neuronFigures_Or2-v-LOC551704_clusters.pdf")
ggsave(filename, limitsize = FALSE, units = "px", width = 2000, height = 2000,p3)

# FigS3D: Orco track plot 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN) <- "peaks"
# first compute the GC content for each peak
ORN <- RegionStats(ORN, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
Annotation(ORN)$tx_id <-Annotation(ORN)$gene_name
#features<-c("Or2","LOC411079","LOC410151","LOC406073","LOC409780","Obp5","Obp11","Obp4")
# link peaks to genes

ORN <- LinkPeaks(
  object = ORN,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = OrcoL
)
######Visulize track and RNA exp######
idents.plot <- Idents(ORN)

pdf("./00_Figure/FigS3D-Orco_gene-peaktrack-RNAexp-WNN.pdf",height=12,width=6)
# Or2
p1 <- CoveragePlot(
  object = ORN,
  region = "Group1-5723000-5724000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)


# LOC552552
p2 <- CoveragePlot(
  object = ORN,
  region = "Group8-7228000-7230000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)


# LOC726019
p3 <- CoveragePlot(
  object = ORN,
  region = "Group11-16143000-16145000",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)
#LOC726019

p4 <- CoveragePlot(
  object = ORN,
  region = "Group14-4996500-4997500",
  #extend.upstream = 0000,
  annotation=TRUE,
  peaks = F,
  #extend.downstream = -30000,
  links=F
)
set<-c(myUmapcolors,myUmapcolors)
p1<-p1& scale_fill_manual(values=set)
p2<-p2& scale_fill_manual(values=set) & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set) & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p4<-p4& scale_fill_manual(values=set) & theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p1|p2|p3|p4
dev.off()

