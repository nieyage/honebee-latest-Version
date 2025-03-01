# enhance the C1234 DEG and heatmap
library(Signac)
library(Seurat)
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
library(ggplot2)
library(dplyr)
obj<- readRDS("./00_Figure/Fig4/Fig4-last-data-obj.rds")
Idents(obj)<-obj$group_manully
obj_features<- c("Or63-b","LOC410603","LOC107963999","LOC100578045")
DefaultAssay(obj)<- "SCT"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
markers<- markers[!duplicated(rownames(markers)),]
write.csv(markers,"./00_Figure/Fig4/Fig4-DEG_markers-DEG_heatmap_need2select_showmarkers.csv")
# downsample C1 to 50 
C1_barcode<-rownames(obj@meta.data[obj$group_manully=="C1",])
C2_barcode<-rownames(obj@meta.data[obj$group_manully=="C2",])
C3_barcode<-rownames(obj@meta.data[obj$group_manully=="C3",])
C4_barcode<-rownames(obj@meta.data[obj$group_manully=="C4",])
down_C1<- sample(C1_barcode,50)
DefaultAssay(obj)<- "SCT"
obj<- ScaleData(obj)
down_obj<-subset(obj,cell=c(down_C1,C2_barcode,C3_barcode,C4_barcode))
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
markers$cluster<- factor(markers$cluster,levels=c("C1","C2","C3","C4"))
#top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
markers<- markers[order(markers$cluster),]
duplicated_genes <- markers$gene[duplicated(markers$gene) | duplicated(markers$gene, fromLast = TRUE)]

LOC413399
LOC409565
Foxp
# 删除包含重复基因的行
new_data <- markers[!markers$gene %in% duplicated_genes, ]

markers<- new_data
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC);

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC);

C1_C2 <- FindMarkers(down_obj,ident.1="C1",ident.2="C2",  min.pct = 0.05, logfc.threshold = 0.3)
C2_gene<- rownames(C1_C2[C1_C2$avg_log2FC<0,])
order<- C1_C2[C1_C2$avg_log2FC<0,]
order<- order[order(order$avg_log2FC,decreasing=F),]
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[7:9])

last<- c("Or63-b","LOC413153","LOC408914","LOC551098","LOC413399",#C1
  "LOC410603","LOC410307","LOC411277","LOC413968","LOC102653594",#C2
  "LOC551123","LOC113218600","LOC102653599","Foxp","LOC410246",#C3
  "LOC409565","LOC724195","LOC408935","LOC413159","LOC724962" #C4
  )
pdf("./00_Figure/MBE/SCT-DEG_markers-DEG_heatmap.pdf",width=5,height=4)
# DoHeatmap(object = obj,features=markers$gene,label=T, group.colors =color_for_cluster,
#   disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
# DoHeatmap(object = obj,features=top10$gene,label=T, group.colors =color_for_cluster,
#   disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
# DoHeatmap(object = obj,features=top5$gene,label=T, group.colors =color_for_cluster,
#   disp.min = -1,disp.max = 1,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
# DoHeatmap(object = down_obj,features=top15$gene,label=T, group.colors =color_for_cluster,
#   disp.min = -1.5,disp.max = 1.5,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
# DoHeatmap(object = down_obj,features=markers$gene,label=T, group.colors =color_for_cluster,
#   disp.min = -1.5,disp.max = 1.5,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
# DoHeatmap(object = down_obj,features=top10$gene,label=T, group.colors =color_for_cluster,
#   disp.min = -1.5,disp.max = 1.5,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = down_obj,features=last,label=T, group.colors =color_for_cluster,
  disp.min = -1.5,disp.max = 1.5,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
# DoHeatmap(object = down_obj,features=unique(c(C2_gene)),label=T, group.colors =color_for_cluster,
#   disp.min = -1.5,disp.max = 1.5,size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

pdf("./00_Figure/MBE/SCT-DEG_markers-DEG_heatmap.pdf",width=5,height=8)
DoHeatmap(object = down_obj,features=unique(c(C2_gene)),label=T, group.colors =color_for_cluster,
  size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();
pdf("./00_Figure/MBE/dittoSeq—SCT-DEG_markers-DEG_heatmap.pdf",width=5,height=4)

library(dittoSeq)
dittoHeatmap(down_obj,heatmap.colors = colorRampPalette(solarExtra[3:8])(50), top5$gene,annot.by = c("group_manully"))
dittoHeatmap(down_obj,heatmap.colors = colorRampPalette(solarExtra[3:8])(50), top15$gene,annot.by = c("group_manully"))
dittoHeatmap(down_obj,heatmap.colors = colorRampPalette(solarExtra[3:8])(50), top10$gene,annot.by = c("group_manully"))
dev.off()

library(Scillus)

pdf("./00_Figure/MBE/scillus—SCT-DEG_markers-DEG_heatmap.pdf",width=5,height=4)
plot_heatmap(dataset = down_obj, 
             markers = top5$gene,
             sort_var = c("group_manully"),
             anno_var = c("group_manully"),
             anno_colors = list("Set2"))
plot_heatmap(dataset = down_obj, 
             markers = top15$gene,
             sort_var = c("group_manully"),
             anno_var = c("group_manully"),
             anno_colors = list("Set2"))
plot_heatmap(dataset = down_obj, 
             markers = top10$gene,
             sort_var = c("group_manully"),
             anno_var = c("group_manully"),
             anno_colors = list("Set2"))
dev.off()


# Find the DEP among the 4 cluster:
DefaultAssay(obj)<- "peaks_ORN_subcluster"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
write.csv(markers,"./00_Figure/Fig4/Fig4-DEG_markers-DEP_heatmap.csv")

sambaNight<- c("#1873CC", "#1798E5" ,"#00BFFF", "#4AC596" ,"#00CC00" ,"#A2E700", "#FFFF00" ,"#FFD200","#FFA500")
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

obj<- ScaleData(obj)
DefaultAssay(down_obj)<- "peaks_ORN_subcluster"
down_obj<- ScaleData(down_obj)
markers <- FindAllMarkers(down_obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
markers$cluster<- factor(markers$cluster,levels=c("C1","C2","C3","C4"))
markers<- markers[order(markers$cluster),]


pdf("./00_Figure/MBE/Fig4-DEP_markers-DEP_heatmap.pdf",width=7,height=5)
DoHeatmap(object = down_obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -1.5,disp.max = 1.5,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
DoHeatmap(object = down_obj,features=markers$gene,label=T, group.colors =color_for_cluster,
  disp.min = -2,disp.max = 2,size = 2,group.by = "Annotation") + scale_fill_gradientn(colors = blueYellow[3:8])
dev.off();

C1_C2 <- FindMarkers(down_obj,ident.1="C1",ident.2="C2",  min.pct = 0.05, logfc.threshold = 0.3)
C2_gene<- rownames(C1_C2[C1_C2$avg_log2FC<0,])
order<- C1_C2[C1_C2$avg_log2FC<0,]
order<- order[order(order$avg_log2FC,decreasing=F),]


# MP SP coexp OR exp ranking 
# Fig3: Cis-Elements Orchestrating Olfactory Receptor Co-Expression 
library(Signac)
library(Seurat)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Amel.antenan)
library(patchwork)

input_data<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/MBE/MP_SP_input.csv")
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak_latest.rds")
UMI_data<- ORN@assays$raw_RNA;
ORN_count<- UMI_data[input_data$gene,]
#gene_total_UMI<- rowSums(ORN_count)
#gene_total_UMI<- gene_total_UMI[order(gene_total_UMI,decreasing=TRUE)]
#data<- data.frame(gene_total_UMI=gene_total_UMI,gene=names(gene_total_UMI))
input_data$ranking<- NA;
input_data$gene_total_UMI<- NA;

for (i in 1:nrow(input_data)){
	cell_group<- input_data[i,3]
	cell<- rownames(ORN@meta.data[ORN@meta.data$cell_group==cell_group,])
	gene<- input_data[i,1]
	tmp_data<- UMI_data[,cell]
	gene_total_UMI<- rowSums(tmp_data)
	gene_total_UMI<- gene_total_UMI[order(gene_total_UMI,decreasing=TRUE)]
	data<- data.frame(gene_total_UMI=gene_total_UMI,gene=names(gene_total_UMI),ranking=c(1:length(gene_total_UMI)))
	input_data$ranking[i]<- data[gene,3]
	input_data$gene_total_UMI[i]<- data[gene,1]
}
library(ggpubr)
input_data$type<-factor(input_data$type,levels=c("MP","SP"))
colors_for_exp_pattern<- c( "#DE7C5B","#FBD277")

my_comparisons <- list(c("MP", "SP"))
pdf("./00_Figure/MBE/MPvsSP_UMI.pdf",width=6,height=4)
p1<- ggboxplot(input_data, x="type", y="ranking", color = "type",width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons,method='t.test')+theme(legend.position="none")+ylab("The ranking of total OR UMI")+
scale_color_manual(values = colors_for_exp_pattern)+
  geom_jitter(aes(color = factor(type)), width = 0.1, alpha = 1)

p2<- ggboxplot(input_data, x="type", y="gene_total_UMI", color = "type",width=0.6, notch = F)+
stat_compare_means(comparisons = my_comparisons,method='t.test')+theme(legend.position="none")+ylab("Total UMI of OR")+
scale_color_manual(values = colors_for_exp_pattern)+
  geom_jitter(aes(color = factor(type)), width = 0.1, alpha = 1)

p1+ p2
dev.off()


# step1: plot the OR ggtree 
# sequence tree 
library(Biostrings)
library(muscle)
library(ape);
library(ggtree)
library(tidytree)
library(ggtreeExtra)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(reshape2)
library(seqinr)
library(phangorn)
library(msa)

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', 
                    '#E0D4CA', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', 
                    '#6778AE', '#B53E2B', '#DCC1DD', '#CCE0F5', '#625D9E', 
                    '#68A180', '#968175', '#FF9999', '#344CB7', '#FFCC1D', 
                    '#24A19C', '#FF9A76',"#BC80BD", "#CCEBC5", "#FFED6F", 
                    "#E95C59", "#476D87",'#9FA3A8')
OR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
GR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/GR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
IR_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/IR_transcript_pep.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
supply_fasta<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/supply.aa", format="fasta",nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
pub<-readAAStringSet("/md01/nieyg/ref/10X/Amel_HAv3.1/OR_transcript_pep_tree.aa")
tree_anno<- read.csv("./00_Figure/FigS5/FigS5_OR_gene_tree_annotation_info.csv")

#OR2 is placed in the last column;
all_OR_gene_fasta<- c(pub,supply_fasta[c(1,4,5)])
all_OR_gene_fasta<- all_OR_gene_fasta[which(names(all_OR_gene_fasta)%in% tree_anno$OR_gene),]
names(all_OR_gene_fasta) <- tree_anno$last_name[match(names(all_OR_gene_fasta),tree_anno$OR_gene)]

tree_anno<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/FigS5/OR_gene_tree_anno2.csv")
rownames(tree_anno)<- tree_anno$last_name
data<-  as.data.frame(tree_anno[,c(5,6,8)])
data[which(data$exp_pattern2=="undefine"),]$exp_pattern2="others"
data[which(data$exp_pattern2=="single_OR"),]$exp_pattern2="others"

data_long<- melt(data, id.vars = c("last_name"), #需保留的不参与聚合的变量列名
                  measure.vars = c('exp_pattern2'),
                  variable.name = c('POS'),#聚合变量的新列名
                  value.name = 'value')
order<- c("co-exp","others")
data_long$value<- factor(data_long$value,levels=order)


alignment <- msa(all_OR_gene_fasta, method = "Muscle")
phangAlign <- as.phyDat(alignment, type = "AA")
mt = modelTest(phangAlign)
# extract best model
best_model <- as.pml(mt)
# 计算距离矩阵
dm <- dist.ml(phangAlign, model=best_model$model, 
              bf=best_model$bf, Q=best_model$Q,
              k=best_model$k, shape=best_model$shape)
treeNJ <- NJ(dm)
write.tree(treeNJ, file = "cyp_nj.nwk")
fit <- pml(treeNJ, phangAlign)
fitML <- optim.pml(fit, model=best_model$model, 
                   optInv = TRUE, optGamma = TRUE)
write.tree(fitML$tree, file = "cyp_ml.nwk")
pdf("./00_Figure/MBE/NJ_OR_sequence_protein_similarity-tree_add_anno.pdf",width=10,height=15)
p1<- ggtree(treeNJ,color="black",linetype=1,size=1.5,ladderize = T) + geom_tiplab(hjust = -0.05,size=5,fontface="italic")
# p1
p2<- p1+ new_scale_fill() + 
      geom_fruit(
          data=data_long,
          geom=geom_tile,
          mapping=aes(y=last_name, x=POS,fill=value),
          offset=0.04,   # The distance between external layers, default is 0.03 times of x range of tree.
          pwidth=0.1 # width of the external layer, default is 0.2 times of x range of tree.
      ) +
      scale_fill_manual(
          values=c("#DE7C5B", "grey"),
          guide=guide_legend(keywidth=1, keyheight=1, order=3)
      ) 
p2;
# create the chemoreceptor info 
value_colors <- c("#DE7C5B", "grey")
names(value_colors) <- order
groupInfo <- split(data_long$last_name, data_long$value)
tree <- groupOTU(treeNJ, groupInfo)
p1<- ggtree(tree,aes(color=group),linetype=1,size=1.5,ladderize = T) + 
scale_color_manual(values = value_colors) +
geom_tiplab(hjust = -0.05,size=4,fontface="italic")
p1
dev.off()



pdf("./00_Figure/MBE/fitML_OR_sequence_protein_similarity-tree_add_anno.pdf",width=10,height=15)
p1<- ggtree(fitML$tree,color="black",linetype=1,size=1.5,ladderize = T) + geom_tiplab(hjust = -0.05,size=5,fontface="italic")
# p1
p2<- p1+ new_scale_fill() + 
      geom_fruit(
          data=data_long,
          geom=geom_tile,
          mapping=aes(y=last_name, x=POS,fill=value),
          offset=0.04,   # The distance between external layers, default is 0.03 times of x range of tree.
          pwidth=0.1 # width of the external layer, default is 0.2 times of x range of tree.
      ) +
      scale_fill_manual(
          values=c("#129FAF", "#DE7C5B","#FBD277","grey"),
          guide=guide_legend(keywidth=1, keyheight=1, order=3)
      ) 
p2;
# create the chemoreceptor info 
value_colors <- c("#129FAF", "#DE7C5B", "#FBD277", "grey")
names(value_colors) <- order
groupInfo <- split(data_long$last_name, data_long$value)
tree <- groupOTU(fitML$tree, groupInfo)
p1<- ggtree(tree,aes(color=group),linetype=1,size=1.5,ladderize = T) + 
scale_color_manual(values = value_colors) +
geom_tiplab(hjust = -0.05,size=4,fontface="italic")
p1
dev.off()
