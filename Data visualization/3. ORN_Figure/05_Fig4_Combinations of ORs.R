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
ORN<- readRDS("./05_ORN_cluster/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_order_by_tree_recall_peak.rds")
DefaultAssay(ORN)<-"raw_RNA"
dotplot_data<-read.csv("./05_ORN_cluster/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
OR_Freq<- as.data.frame(table(dotplot_data$features.plot))
feature<-OR_Freq[OR_Freq$Freq>1,]$Var1
combination_data<- dotplot_data[dotplot_data$features.plot%in% feature,]
#write.csv(combination_data,"./05_ORN_cluster/02_second_cluster/06_rm_without_power/combination_data.csv")
combination_data<- read.csv("./05_ORN_cluster/02_second_cluster/06_rm_without_power/combination_data.csv")


# show a typical combination;
# select a beautiful track to show :
library(ggplot2)
library(gggenes)
library(pheatmap)
library(dittoSeq)
library(cowplot)
library(lsa)
library(ggpubr)
library(RColorBrewer)
log2FCdata<-data.frame();
DefaultAssay(ORN)<- "integratedRNA_onecluster"
pdf("./05_ORN_cluster/03_coexp_cluster/combination_data_multi_OR_with_powerful.pdf",width=14,height=6)
for (cluster in unique(combination_data$combination)){
  print(cluster)
  combination_group<- combination_data[which(combination_data$combination==cluster),];
  cluster_info<- unique(combination_group$id)
  obj<-subset(ORN,idents=cluster_info);
  obj_features<- unique(combination_group$features.plot)
  if(length(obj_features)> 1){
# add max_exp OR label for each cell
  ORN_count<-obj@assays$raw_RNA
  barcode_label<-data.frame(barcode=colnames(obj),label=obj$subcluster)
  ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
  ORN_matrix<-as.matrix(ORN_count)
  ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
  ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
  barcode_label<-barcode_label[colnames(ORN_matrix),]
  barcode_label<-barcode_label[order(barcode_label$label),]
# p1 cell cosine simility heatmap 
    embeddings <- Embeddings(object = ORN, reduction = "pca")[,1:50]
    embeddings <- embeddings[barcode_label$barcode,]
    trans_dist <- 1-cosine(t(embeddings))
    barcode_label_pheatmap<-data.frame(label=barcode_label$label)
    rownames(barcode_label_pheatmap)<-barcode_label$barcode
    col<-brewer.pal(12,"Set3")[1:length(unique(barcode_label_pheatmap$label))]
    names(col)<-unique(barcode_label_pheatmap$label)
    ann_colors= list(label = col)
    #trans_dist<-trans_dist[rownames(barcode_label_pheatmap),rownames(barcode_label_pheatmap)]
    p1<-pheatmap(trans_dist,
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("#4B0C1A", "#C25B3F","#F1ECEC", "#3082BD","#1C214F"))(100),
             annotation_col = barcode_label_pheatmap,
             annotation_colors = ann_colors,
             annotation_row = barcode_label_pheatmap,
             annotation_legend = TRUE,
             show_rownames=F,
             show_colnames=F)
# calculate the cosine simility within group and between groups;
    rownames(barcode_label)<-barcode_label$barcode
    within_group<-c()
    between_group<-c()
    for (i in 1:nrow(trans_dist)){
    for (j in 1:ncol(trans_dist)){
        if(i!=j){
            if(barcode_label[rownames(trans_dist)[i],2]==barcode_label[colnames(trans_dist)[j],2]){
                within_group<-c(within_group,trans_dist[i,j]);
            }
            else{between_group<-c(between_group,trans_dist[i,j])}
        }
    }
    }
# calculate the FC 
    if(!is.null(median(between_group))){
    log2FC<-log2(median(between_group))-log2(median(within_group));
    test<-wilcox.test(within_group,between_group);
    pvalue<-test$p.value;
    data_subset<-data.frame(cluster,log2FC,pvalue)
    log2FCdata<-rbind(log2FCdata,data_subset)
    }
# plot density line 
# manage data
    type<-c(rep("within-OR",length(within_group)),rep("between-OR",length(between_group)))
    var<-c(within_group,between_group)
    data<-data.frame(type,var)
    data$type<-factor(data$type,levels=c("within-OR","between-OR"))
    p2<-ggplot(data, aes(x=var, fill=type)) + xlab("Transcriptome distance")+
                  geom_density(alpha=.25) + theme_classic() 
    # t-test
    p3 <- ggboxplot(data, x="type", y="var", color = "type")+stat_compare_means()+guides(fill = "none")
#raw counts heatmap 
# Heat map of expression  value for Or25-27 for cluster 5_1,5_2,14-1 and others( random select a few as control) 
    obj_barcode<-colnames(obj)
    all_barcode<-colnames(ORN)
    random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
    obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
    subcluster<- obj$subcluster
      for (i in 1:length(subcluster)){
        if(subcluster[i] %in% cluster_info){subcluster[i]=subcluster[i]}
        else{subcluster[i]="other"}
      }
    obj$subcluster<- subcluster
    Idents(obj)<-obj$subcluster
    DefaultAssay(obj)<-"raw_RNA"
    obj_data<-as.data.frame(t(as.matrix(obj@assays$raw_RNA[obj_features,])))
    barcode_info<-data.frame(subcluster=obj$subcluster,barcode=colnames(obj))
    rownames(barcode_info)<-colnames(obj)
    barcode_info<- barcode_info[order(barcode_info$subcluster),]
    obj_data<-obj_data[rownames(barcode_info),]
    #barcode_info<- as.data.frame(barcode_info[,-2])
    col<-c(col,"grey")
    names(col)<-unique(subcluster)
    ann_colors2= list(label = col)
    #barcode_info<-barcode_info[,-2]
    barcode_info_h<-data.frame(subcluster=obj$subcluster)
    rownames(barcode_info_h)<-colnames(obj)
    p4<-pheatmap(t(obj_data),
             cluster_cols = F,
             cluster_rows = F,
             color = colorRampPalette(c("white", "red"))(100),
             annotation_col = barcode_info_h,
             annotation_colors = ann_colors2,
             #annotation_row = barcode_info,
             annotation_legend = TRUE,
             show_rownames=T,
             show_colnames=F
      )
    top_right<-plot_grid(p2,p3,labels = c(" "," "),rel_widths = c(2, 1))
    right<-plot_grid(top_right,p4$gtable,ncol = 1,labels = c(" "," "))
    last<-plot_grid(p1$gtable, right, labels = c(' ', ''), label_size = 12, ncol = 2)
    title <- ggdraw() + 
      draw_label(
        paste("combination",cluster,":","log2FC=",log2FC),
        fontface = 'bold',
        x = 0,
        hjust = 0
      )
    add_title<-plot_grid(title, last,ncol = 1 , rel_heights = c(0.1, 1) )
    print(add_title)
    }}
dev.off()


# # show a typical combination track plot 
# Extract result pages
library(pdftools)
# inputs
infile <- "./05_ORN_cluster/03_coexp_cluster/combination_data_multi_OR_with_powerful.pdf"  # input pdf
outfile <-  "./05_ORN_cluster/03_coexp_cluster/combination_data_multi_OR_with_powerful_extract.pdf"
num <- pdf_length(infile)/3
pdf_subset(infile, pages = (1:num)*3, output = outfile)

# Fig4H coexp track plot
## Track from scATAC-seq for all multiple cluster 
library(BSgenome.Amel.HAv3.1.update.chemoreceptor)
DefaultAssay(ORN)<-"peaks_ORN_subcluster"
pdf("./05_ORN_cluster/03_coexp_cluster/combination_data_multi_OR_with_powerful_trackplot.pdf",width=10,height=10)
for (cluster in unique(combination_data$combination)){
  print(cluster)
  combination_group<- combination_data[which(combination_data$combination==cluster),];
  cluster_info<- unique(combination_group$id)
  obj<-subset(ORN,idents=cluster_info);
  obj_features<-obj_features<- unique(combination_group$features.plot)
  obj_barcode<-colnames(obj)
  all_barcode<-colnames(ORN)
  random_barcode<-sample(setdiff(all_barcode,obj_barcode),length(obj_barcode))
  obj<-subset(ORN,cells=c(obj_barcode,random_barcode))
  subcluster<- obj$subcluster
    for (i in 1:length(subcluster)){
      if(subcluster[i] %in% cluster_info){subcluster[i]=subcluster[i]}
      else{subcluster[i]="other"}
    }
  Idents(obj)<-subcluster
  # first compute the GC content for each peak
  obj <- RegionStats(obj, genome = BSgenome.Amel.HAv3.1.update.chemoreceptor)
  Annotation(obj)$tx_id <-gsub("_g","-g",Annotation(obj)$gene_name)
  Annotation(obj)$tx_id <-Annotation(obj)$transcript_id
  # link peaks to genes
  obj <- LinkPeaks(
    object = obj,
    peak.assay = "peaks_ORN_subcluster",
    expression.assay = "raw_RNA",
    genes.use = obj_features
  )
  ######Visulize track and RNA exp######
  idents.plot <- Idents(obj)
  # plot region 
  start<- min(start(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  end  <- max(end(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]))
  seq  <- as.character(Annotation(obj)[which(Annotation(obj)$gene_name%in%obj_features),]@seqnames@values)
  ranges.show <- paste(seq,start,end,sep="-")
  col<-brewer.pal(12,"Set3")[1:length(cluster_info)]
  col<- c(col,"grey")
  p1<-CoveragePlot(
    object = obj,
    region = ranges.show,
    window = 150,
    extend.upstream = 100,
    annotation = TRUE,
    extend.downstream = 100,
    tile = TRUE,
    tile.size = 100,
    tile.cells = 30,
    links=F
  )
  #p1<- p1&scale_fill_manual(values = col )
  print(p1)
}
dev.off()


# UpsetR plot of Observation and Expectation
library(UpSetR)
library(reshape2)  
pdf("./05_ORN_cluster/03_coexp_cluster/combination_data_multi_OR_with_powerful_upsetRplot_Observation.pdf",width=10,height=10)
for (cluster in unique(combination_data$combination)){
  print(cluster)
  combination_group<- combination_data[which(combination_data$combination==cluster),];
  cluster_info<- unique(combination_group$id)
  obj<-subset(ORN,idents=cluster_info);
  obj_features<-obj_features<- unique(combination_group$features.plot)
  obj_barcode<-colnames(obj)
  if(length(obj_features)>1){
  # 1: Observation
  ORN_count<-obj@assays$RNA
  ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
  ORN_matrix<-as.matrix(ORN_count)
  ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
  ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
  data<- as.data.frame(ORN_matrix)
  data$features<- rownames(data)
  data_long<-melt(data, id.vars = c("features"),
                 measure.vars = c(colnames(data)[-length(colnames(data))]),
                 variable.name = c('barcode'),
                 value.name = 'value')
  data_long<- data_long[data_long$value>0,]
  listInput <- split(data_long$barcode,data_long$features)
  p_observation<- upset(fromList(listInput), nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))

print(p_observation)
}}
dev.off()

  # 2: Expectation by random dropout
pdf("./05_ORN_cluster/03_coexp_cluster/combination_data_multi_OR_with_powerful_upsetRplot_Expectation.pdf",width=10,height=10)
for (cluster in unique(combination_data$combination)){
  print(cluster)
  combination_group<- combination_data[which(combination_data$combination==cluster),];
  cluster_info<- unique(combination_group$id)
  obj<-subset(ORN,idents=cluster_info);
  obj_features<-unique(combination_group$features.plot)
  obj_barcode<-colnames(obj)
  ORN_count<-obj@assays$RNA
  ORN_count<-ORN_count[which(rownames(ORN_count)%in%obj_features),]
  ORN_matrix<-as.matrix(ORN_count)
  ORN_matrix<-ORN_matrix[,colSums(ORN_matrix)>0]
  ORN_matrix<-ORN_matrix[rowSums(ORN_matrix)>0,]
  data<- as.data.frame(ORN_matrix)
  data$features<- rownames(data)
  data_long<-melt(data, id.vars = c("features"),
                 measure.vars = c(colnames(data)[-length(colnames(data))]),
                 variable.name = c('barcode'),
                 value.name = 'value')
  data_long<- data_long[data_long$value>0,]
  listInput <- split(data_long$barcode,data_long$features)
  total_cell <- length(colnames(ORN_matrix))
  if(length(obj_features)==4){
    gene1<- obj_features[1]
    gene2<- obj_features[2]
    gene3<- obj_features[3]
    gene4<- obj_features[4]
    gene1_dropout <- 1-length(listInput[[1]])/length(colnames(ORN_matrix))
    gene2_dropout <- 1-length(listInput[[2]])/length(colnames(ORN_matrix))
    gene3_dropout <- 1-length(listInput[[3]])/length(colnames(ORN_matrix))
    gene4_dropout <- 1-length(listInput[[4]])/length(colnames(ORN_matrix))
    input <- c(
      "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout*gene4_dropout,
      "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout*gene4_dropout,
      "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout*gene4_dropout,
      "gene4" =total_cell*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene3_dropout,
      "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout*gene4_dropout, 
      "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout*gene4_dropout, 
      "gene1&gene4" =total_cell*(1-gene1_dropout)*(1-gene4_dropout)*gene2_dropout*gene3_dropout, 
      "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout*gene4_dropout, 
      "gene2&gene4" =total_cell*(1-gene2_dropout)*(1-gene4_dropout)*gene1_dropout*gene3_dropout, 
      "gene3&gene4" =total_cell*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene2_dropout, 
      "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*gene4_dropout,
      "gene1&gene2&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene4_dropout)*gene3_dropout,
      "gene1&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout, 
      "gene2&gene3&gene4"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout,
      "gene1&gene2&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout))}
  if(length(obj_features)==5){
    gene1_dropout <- 1-length(listInput[[1]])/length(colnames(ORN_matrix))
    gene2_dropout <- 1-length(listInput[[2]])/length(colnames(ORN_matrix))
    gene3_dropout <- 1-length(listInput[[3]])/length(colnames(ORN_matrix))
    gene4_dropout <- 1-length(listInput[[4]])/length(colnames(ORN_matrix))
    gene5_dropout <- 1-length(listInput[[5]])/length(colnames(ORN_matrix))
    input <- c(
      "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout*gene4_dropout*gene5_dropout,
      "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout*gene4_dropout*gene5_dropout,
      "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout*gene4_dropout*gene5_dropout,
      "gene4" =total_cell*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene3_dropout*gene5_dropout,
      "gene5" =total_cell*(1-gene5_dropout)*gene1_dropout*gene2_dropout*gene3_dropout*gene4_dropout,
      "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout*gene4_dropout*gene5_dropout, 
      "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout*gene4_dropout*gene5_dropout, 
      "gene1&gene4" =total_cell*(1-gene1_dropout)*(1-gene4_dropout)*gene2_dropout*gene3_dropout*gene5_dropout, 
      "gene1&gene5" =total_cell*(1-gene1_dropout)*(1-gene5_dropout)*gene2_dropout*gene3_dropout*gene4_dropout, 
      "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout*gene4_dropout*gene5_dropout, 
      "gene2&gene4" =total_cell*(1-gene2_dropout)*(1-gene4_dropout)*gene1_dropout*gene3_dropout*gene5_dropout, 
      "gene2&gene5" =total_cell*(1-gene2_dropout)*(1-gene5_dropout)*gene1_dropout*gene3_dropout*gene4_dropout, 
      "gene3&gene4" =total_cell*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene2_dropout*gene5_dropout,
      "gene3&gene5" =total_cell*(1-gene3_dropout)*(1-gene5_dropout)*gene1_dropout*gene2_dropout*gene4_dropout,
      "gene4&gene5" =total_cell*(1-gene4_dropout)*(1-gene5_dropout)*gene1_dropout*gene2_dropout*gene3_dropout,
      "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*gene4_dropout*gene5_dropout,
      "gene1&gene2&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene4_dropout)*gene3_dropout*gene5_dropout,
      "gene1&gene2&gene5"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene5_dropout)*gene3_dropout*gene4_dropout,
      "gene1&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout*gene5_dropout,
      "gene1&gene3&gene5"=total_cell*(1-gene1_dropout)*(1-gene3_dropout)*(1-gene5_dropout)*gene2_dropout*gene4_dropout,
      "gene2&gene3&gene4"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout*gene5_dropout,
      "gene2&gene3&gene5"=total_cell*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene5_dropout)*gene1_dropout*gene4_dropout,
      "gene3&gene4&gene5"=total_cell*(1-gene3_dropout)*(1-gene4_dropout)*(1-gene5_dropout)*gene1_dropout*gene2_dropout,
      "gene1&gene2&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene5_dropout,
      "gene1&gene2&gene3&gene5"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene5_dropout)*gene4_dropout,
      "gene1&gene2&gene5&gene4"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene5_dropout)*(1-gene4_dropout)*gene3_dropout,
      "gene1&gene5&gene3&gene4"=total_cell*(1-gene1_dropout)*(1-gene5_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene2_dropout,
      "gene5&gene2&gene3&gene4"=total_cell*(1-gene5_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*gene1_dropout,
      "gene1&gene2&gene3&gene4&gene5"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout)*(1-gene4_dropout)*(1-gene5_dropout))}
  if(length(obj_features)==3){
        gene1<- obj_features[1]
        gene2<- obj_features[2]
        gene3<- obj_features[3]
        gene1_dropout <- 1-length(listInput[[1]])/length(colnames(ORN_matrix))
        gene2_dropout <- 1-length(listInput[[2]])/length(colnames(ORN_matrix))
        gene3_dropout <- 1-length(listInput[[3]])/length(colnames(ORN_matrix))
         input <- c("gene1" =total_cell*(1-gene1_dropout)*gene2_dropout*gene3_dropout,
          "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout*gene3_dropout,
          "gene3" =total_cell*(1-gene3_dropout)*gene1_dropout*gene2_dropout,
          "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout)*gene3_dropout, 
          "gene1&gene3" =total_cell*(1-gene1_dropout)*(1-gene3_dropout)*gene2_dropout, 
          "gene2&gene3" =total_cell*(1-gene2_dropout)*(1-gene3_dropout)*gene1_dropout, 
          "gene1&gene2&gene3"=total_cell*(1-gene1_dropout)*(1-gene2_dropout)*(1-gene3_dropout))}
  if(legend(obj_features)==2){
        gene1<- obj_features[1]
    gene2<- obj_features[2]
    gene1_dropout <- 1-length(listInput[[1]])/length(colnames(ORN_matrix))
    gene2_dropout <- 1-length(listInput[[2]])/length(colnames(ORN_matrix))
     input <- c(
      "gene1" =total_cell*(1-gene1_dropout)*gene2_dropout,
      "gene2" =total_cell*(1-gene2_dropout)*gene1_dropout,
      "gene1&gene2" =total_cell*(1-gene1_dropout)*(1-gene2_dropout))}
  data <- UpSetR::fromExpression(input)
  colnames(data)<- names(listInput)
  p_expectation<- upset(data, nintersects = 30, mb.ratio = c(0.5, 0.5), order.by = c("freq", "degree"), decreasing = c(FALSE,FALSE))
  print(p_expectation)
}
dev.off()


