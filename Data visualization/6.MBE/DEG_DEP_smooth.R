Idents(down_obj)<-down_obj$group_manully
last<- c("Or63-b","LOC413153","LOC408914","LOC551098","LOC413399",#C1
  "LOC410603","LOC410307","LOC411277","LOC413968","LOC102653594",#C2
  "LOC113219162","LOC411586","LOC413416","Foxp","LOC410246",#C3
  "LOC409565","LOC100577081","LOC551673","LOC413014","LOC102654426" #C4
  )

barcode_label<-data.frame(barcode=colnames(down_obj),label=Idents(down_obj))
DefaultAssay(down_obj)<-"raw_RNA"

DefaultAssay(obj)<- "SCT"
obj<- ScaleData(obj)
down_obj<-subset(obj,cell=c(down_C1,C2_barcode,C3_barcode,C4_barcode))
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)
markers$cluster<- factor(markers$cluster,levels=c("C1","C2","C3","C4"))
#top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
markers<- markers[order(markers$cluster),]
duplicated_genes <- markers$gene[duplicated(markers$gene) | duplicated(markers$gene, fromLast = TRUE)]


new_data <- markers[!markers$gene %in% duplicated_genes, ]
markers<- new_data

all_DEG<- c(markers$gene[1:22],"LOC410603","LOC410307","LOC411277","LOC413968","LOC102653594",markers$gene[24:53])

write.csv(markers[all_DEG,],"C1234_DEG_cluster.csv")
obj_data<-as.data.frame(t(as.matrix(down_obj@assays$SCT[all_DEG,])))
barcode_label<- barcode_label[order(barcode_label$label),]
obj_data<- obj_data[rownames(barcode_label),]


C39_data<-obj_data[rownames(barcode_label[barcode_label$label=="C1",]),]
C40_data<-obj_data[rownames(barcode_label[barcode_label$label=="C2",]),]
C41_data<-obj_data[rownames(barcode_label[barcode_label$label=="C3",]),]
C42_data<-obj_data[rownames(barcode_label[barcode_label$label=="C4",]),]

C39_data<- C39_data[which(rowSums(C39_data)>0),]
C40_data<- C40_data[which(rowSums(C40_data)>0),]
C41_data<- C41_data[which(rowSums(C41_data)>0),]
C42_data<- C42_data[which(rowSums(C42_data)>0),]

last_data_heatmap<- rbind(C39_data,C40_data,C41_data,C42_data)



#order_heatmap<- read.csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/00_Figure/Fig4/Fig4E_4gene_heatmap_data.csv",row.names=1)
#cell_order<- rownames(order_heatmap)
#cell_order<- cell_order[which(cell_order%in% rownames(last_data_heatmap))]
#DEG_heatmap_data<- last_data_heatmap[cell_order,last]

smooth_column <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 2:(length(col) - 1)) {
    smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
  }
  smoothed[1] <- (col[1] + col[2]) / 2
  smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
  return(smoothed)
}

smooth_column_5cell <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 3:(length(col) - 2)) {
    smoothed[i] <- (col[i - 2] + col[i - 1] + col[i] + col[i + 1] + col[i + 2]) / 5
  }
  smoothed[1] <- (col[1] + col[2] + col[3]) / 3
  smoothed[2] <- (col[1] + col[2] + col[3] + col[4]) / 4
  smoothed[length(col) - 1] <- (col[length(col) - 3] + col[length(col) - 2] + col[length(col) - 1] + col[length(col)]) / 4
  smoothed[length(col)] <- (col[length(col) - 2] + col[length(col) - 1] + col[length(col)]) / 3
  return(smoothed)
}
smooth_column_10cell <- function(col) {
  smoothed <- numeric(length(col))
  
  # 中间部分：10个单元的滑动平均
  for (i in 6:(length(col) - 5)) {
    smoothed[i] <- mean(col[(i - 5):(i + 4)])
  }
  
  # 边界部分的特殊处理
  smoothed[1] <- mean(col[1:3])  # 前3个单元
  smoothed[2] <- mean(col[1:4])  # 前4个单元
  smoothed[3] <- mean(col[1:5])  # 前5个单元
  smoothed[4] <- mean(col[1:6])  # 前6个单元
  smoothed[5] <- mean(col[1:7])  # 前7个单元
  
  smoothed[length(col) - 4] <- mean(col[(length(col) - 6):length(col)])  # 后7个单元
  smoothed[length(col) - 3] <- mean(col[(length(col) - 5):length(col)])  # 后6个单元
  smoothed[length(col) - 2] <- mean(col[(length(col) - 4):length(col)])  # 后5个单元
  smoothed[length(col) - 1] <- mean(col[(length(col) - 3):length(col)])  # 后4个单元
  smoothed[length(col)] <- mean(col[(length(col) - 2):length(col)])  # 后3个单元
  
  return(smoothed)
}

library(pheatmap)
# 对每一列进行平滑
my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[6:8])
barcode_label_pheatmap<-data.frame(label=c(rep("C1",nrow(C39_data)),rep("C2",nrow(C40_data)),rep("C3",nrow(C41_data)),rep("C4",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-c(rownames(C39_data),rownames(C40_data),rownames(C41_data),rownames(C42_data))
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)

DEG_heatmap_data<- last_data_heatmap[,all_DEG]

last_data_heatmap_smoothed_data <- as.data.frame(lapply(DEG_heatmap_data, smooth_column))
rownames(last_data_heatmap_smoothed_data)<- rownames(DEG_heatmap_data)

last_data_heatmap_smoothed_data_10cell <- as.data.frame(lapply(DEG_heatmap_data, smooth_column_10cell))
rownames(last_data_heatmap_smoothed_data_10cell)<- rownames(DEG_heatmap_data)

last_data_heatmap_smoothed_data_5cell <- as.data.frame(lapply(DEG_heatmap_data, smooth_column_5cell))
rownames(last_data_heatmap_smoothed_data_5cell)<- rownames(DEG_heatmap_data)

write.csv(last_data_heatmap,"All_DEG_heatmap_matrix.csv")

barcode_label_pheatmap$barcode<- rownames(barcode_label_pheatmap)
write.csv(barcode_label,"barcode_label.csv")

pdf("./00_Figure/MBE/All_DEG_heatmap.pdf",height=8,width=8)
pheatmap(t(DEG_heatmap_data),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(solarExtra)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(solarExtra)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data_5cell),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(solarExtra)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data_10cell),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(solarExtra)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()

# RNA DEG to peak and peak heatmap 
DefaultAssay(down_obj) <- "peaks_ORN_subcluster"
# first compute the GC content for each peak
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
down_obj <- RegionStats(down_obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
down_obj <- LinkPeaks(
  object = down_obj,
  peak.assay = "peaks_ORN_subcluster",
  expression.assay = "raw_RNA",
  pvalue_cutoff = 1,
  min.cells = 5,
  genes.use = last,score_cutoff = 0
)
lnk <- Links(object = down_obj[["peaks_ORN_subcluster"]])
highest_score_link <- as.data.frame(lnk) %>%
  group_by(gene) %>%
  filter(score == max(score))
highest_score_link<- as.data.frame(highest_score_link)
rownames(highest_score_link)<- highest_score_link$gene;

RNA_DEG_peak <- na.omit(highest_score_link[last,]$peak)


RNA_DEG_peak<- rownames(markers)

DefaultAssay(down_obj)<- "peaks_ORN_subcluster"
down_obj<- ScaleData(down_obj)
markers <- FindAllMarkers(down_obj, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
markers$cluster<- factor(markers$cluster,levels=c("C1","C2","C3","C4"))
markers<- markers[order(markers$cluster),]

top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC);

RNA_DEG_peak<- top5$gene
barcode_label<-data.frame(barcode=colnames(down_obj),label=Idents(down_obj))
DefaultAssay(down_obj)<-"peaks_ORN_subcluster"
obj_data<-as.data.frame(t(as.matrix(down_obj@assays$peaks_ORN_subcluster[RNA_DEG_peak,])))
barcode_label<- barcode_label[order(barcode_label$label),]
obj_data<- obj_data[rownames(barcode_label),]

C39_data<-obj_data[rownames(barcode_label[barcode_label$label=="C1",]),]
C40_data<-obj_data[rownames(barcode_label[barcode_label$label=="C2",]),]
C41_data<-obj_data[rownames(barcode_label[barcode_label$label=="C3",]),]
C42_data<-obj_data[rownames(barcode_label[barcode_label$label=="C4",]),]

C39_data<- C39_data[which(rowSums(C39_data)>0),]
C40_data<- C40_data[which(rowSums(C40_data)>0),]
C41_data<- C41_data[which(rowSums(C41_data)>0),]
C42_data<- C42_data[which(rowSums(C42_data)>0),]

last_data_heatmap<- rbind(C39_data,C40_data,C41_data,C42_data)
write.csv(last_data_heatmap,"All_DEP_heatmap_matrix.csv")


color_for_cluster<- c("#4BA9D1",my47colors[6:8])
barcode_label_pheatmap<-data.frame(label=c(rep("C1",nrow(C39_data)),rep("C2",nrow(C40_data)),rep("C3",nrow(C41_data)),rep("C4",nrow(C42_data))))
rownames(barcode_label_pheatmap)<-c(rownames(C39_data),rownames(C40_data),rownames(C41_data),rownames(C42_data))
col <- color_for_cluster[1:length(unique(barcode_label_pheatmap$label))]
names(col)<-unique(barcode_label_pheatmap$label)
ann_colors= list(label = col)

DEG_heatmap_data<- last_data_heatmap[,RNA_DEG_peak]

last_data_heatmap_smoothed_data <- as.data.frame(lapply(DEG_heatmap_data, smooth_column))
rownames(last_data_heatmap_smoothed_data)<- rownames(DEG_heatmap_data)

last_data_heatmap_smoothed_data_10cell <- as.data.frame(lapply(DEG_heatmap_data, smooth_column_10cell))
rownames(last_data_heatmap_smoothed_data_10cell)<- rownames(DEG_heatmap_data)

last_data_heatmap_smoothed_data_5cell <- as.data.frame(lapply(DEG_heatmap_data, smooth_column_5cell))
rownames(last_data_heatmap_smoothed_data_5cell)<- rownames(DEG_heatmap_data)
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")

pdf("./00_Figure/MBE/top5_DEP_heatmap.pdf",height=8,width=8)
pheatmap(t(DEG_heatmap_data),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(blueYellow)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(blueYellow)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data_5cell),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(blueYellow)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
pheatmap(t(last_data_heatmap_smoothed_data_10cell),
             cluster_cols = F,
             cluster_rows = F,
             scale = "row",
             color = colorRampPalette(blueYellow)(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             #annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()










C3_C4 <- FindMarkers(down_obj,ident.1="C3",ident.2="C4",  min.pct = 0.05, logfc.threshold = 0.3)

order<- C3_C4[C3_C4$avg_log2FC<0,]
order<- order[order(order$avg_log2FC,decreasing=F),]

C4_gene<- rownames(order[1:30,])

order<- C3_C4[C3_C4$avg_log2FC>0,]
order<- order[order(order$avg_log2FC,decreasing=T),]

C3_gene<- rownames(order[1:30,])

pdf("./00_Figure/MBE/C34SCT-DEG_markers-DEG_heatmap.pdf",width=5,height=8)
DoHeatmap(object = down_obj,features=unique(c(C3_gene)),label=T, group.colors =color_for_cluster,
  size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = down_obj,features=unique(c(C4_gene)),label=T, group.colors =color_for_cluster,
  size = 2,group.by = "group_manully") + scale_fill_gradientn(colors = solarExtra[3:8])

dev.off();


my47colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
         '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
         '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
         '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', '#116530', '#678983',
         '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76')
color_for_cluster<- c("#4BA9D1",my47colors[7:9])












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


obj<-subset(ORN,idents="p1:14");
obj_features<- unique(dotplot_data[dotplot_data$id%in%"p1:14",]$features.plot)
DefaultAssay(obj)<-"SCT"
obj_data<-as.data.frame(t(as.matrix(obj@assays$SCT[obj_features,])))

obj_data<-obj_data[order(obj_data$Or25,obj_data$Or26,obj_data$Or27,decreasing=T),]
# 定义平滑函数
smooth_column <- function(col) {
  smoothed <- numeric(length(col))
  for (i in 2:(length(col) - 1)) {
    smoothed[i] <- (col[i - 1] + col[i] + col[i + 1]) / 3
  }
  smoothed[1] <- (col[1] + col[2]) / 2
  smoothed[length(col)] <- (col[length(col) - 1] + col[length(col)]) / 2
  return(smoothed)
}

# 对每一列进行平滑
smoothed_data <- as.data.frame(lapply(obj_data, smooth_column))

pdf("./00_Figure/Fig4/Fig4A-OR25_27_heatmap-smooth.pdf",width=8,height=2)
    pheatmap(t(smoothed_data),
             cluster_cols = F,
             cluster_rows = F,
             border_color=NA,
             color = colorRampPalette(c("white", "#CA0002"))(100),
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
dev.off()


