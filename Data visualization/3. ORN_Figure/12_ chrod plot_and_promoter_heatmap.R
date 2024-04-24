###  The Coexp circle plot in honeybee ###
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
dotplot_data<-read.csv("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/dotplot_data_remove_nopower_latest.csv")
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
dotplot_feature <- ORgene_name_trans[match(dotplot_feature,ORgene_name_trans$OR_gene),]$last_name


ORN <- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_latest.rds")
cluster_info<-as.data.frame(table(dotplot_data$id))
multiOR_cluster<-as.character(cluster_info[cluster_info$Freq>1,1])

multiOR_cluster_OR<- unique(dotplot_data[which(dotplot_data$id%in%multiOR_cluster),]$features.plot)
##########################################################
### chord plot
library(chorddiag)

highReceptorRNA = ORN@assays$raw_RNA@counts[multiOR_cluster_OR,]
highReceptorRNA.data = ORN@assays$raw_RNA@data[multiOR_cluster_OR,]

write.csv(highReceptorRNA,"./highReceptorRNA.csv")
write.csv(highReceptorRNA.data,"./highReceptorRNA.data.csv")

# change the OR gene name 

ORgene_name_trans<- read.csv("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/OR_nameing/OR_gene_naming_result.csv")

highReceptorRNA.expCol=highReceptorRNA[,colSums(highReceptorRNA) != 0]
highReceptorRNA.data.expCol=highReceptorRNA.data[,colSums(highReceptorRNA.data) != 0]

rownames(highReceptorRNA.expCol) <- ORgene_name_trans[match(rownames(highReceptorRNA.expCol),ORgene_name_trans$OR_gene),]$last_name
rownames(highReceptorRNA.data.expCol) <- ORgene_name_trans[match(rownames(highReceptorRNA.data.expCol),ORgene_name_trans$OR_gene),]$last_name

# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.expCol), )
coexp.df <- data.frame(matrix(ncol = 57, nrow = 57))
colnames(coexp.df) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df) <- rownames(highReceptorRNA.expCol)

coexp.df2 <- data.frame(matrix(ncol = 57, nrow = 57))
colnames(coexp.df2) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df2) <- rownames(highReceptorRNA.expCol)

coexp.df3 <- data.frame(matrix(ncol = 57, nrow = 57))
colnames(coexp.df3) <- rownames(highReceptorRNA.expCol)
rownames(coexp.df3) <- rownames(highReceptorRNA.expCol)


# normalized expression
for (i in 1:57) {
  for (j in 1:57) {
     if (i==j) { coexp.df[i,j] <- 0 }
    else {
      lapply(1:length(ncol(highReceptorRNA.data.expCol)), function(x){
        print(c('i = ', i))
        print(c('j = ', j))
        iExpressedCells.mx = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1]
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1])
        # print(jExpressedN)
        gene1<- multiOR_cluster_OR[i]
        gene2<- multiOR_cluster_OR[j]
        cluster_id<-dotplot_data[which(dotplot_data$features.plot%in%c(gene1,gene2)),]$id;
        if(length(which(duplicated(cluster_id)))){
        coexp.df[i,j]<<-jExpressedN
        iExpressedCells.mx2 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 2]
        jExpressedN2=NCOL(iExpressedCells.mx2[, iExpressedCells.mx2[j,] >= 2])
        # print(jExpressedN)
        coexp.df2[i,j]<<-jExpressedN2
        
        iExpressedCells.mx3 = highReceptorRNA.data.expCol[,highReceptorRNA.data.expCol[i,] >= 1.5]
        jExpressedN3=NCOL(iExpressedCells.mx3[, iExpressedCells.mx3[j,] >= 1.5])
        # print(jExpressedN)
	    coexp.df3[i,j]<<-jExpressedN3;}
        else{coexp.df[i,j]<<-0;
        	coexp.df2[i,j]<<-0;
        	coexp.df3[i,j]<<-0;
        }
      })
    }
  }
}

# just keep the coexp_pair
write.csv(coexp.df,"Fig2-coexp_OR_chrod_plot_coexp.df.csv")
write.csv(coexp.df2,"Fig2-coexp_OR_chrod_plot_coexp.df2.csv")
write.csv(coexp.df3,"Fig2-coexp_OR_chrod_plot_coexp.df3.csv")

coexp.df<- read.csv("/Users/fraya/Documents/project/honeybee/Fig2-coexp_OR_chrod_plot_coexp.df.csv",row.names = 1)
coexp.df2<- read.csv("/Users/fraya/Documents/project/honeybee/Fig2-coexp_OR_chrod_plot_coexp.df2.csv",row.names = 1)
coexp.df3<- read.csv("/Users/fraya/Documents/project/honeybee/Fig2-coexp_OR_chrod_plot_coexp.df3.csv",row.names = 1)

# names(hex_codes1) <- colnames(coexp.df)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                    '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
                    '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
                    '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
                    '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
                    "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
                    "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C"  , "#377EB8", "#4DAF4A", 
                    "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
hex_codes1<- myUmapcolors
#
OR<- c("Or1"   ,"Or3"  ,"Or4"  ,"Or5"  , "Or12" , "Or15"  ,"Or16" ,
       "Or25"  ,"Or26" ,"Or27" ,"Or29" , "Or30"  ,"Or31" ,"Or35" , "Or39" ,
       "Or40" ,"Or41" , "Or54"  ,"Or55" ,"Or56" ,"Or57" ,"Or58" ,"Or59" ,         
       "Or60"  ,"Or61" ,"Or62","Or64","Or65","Or67" ,"Or91","Or94","Or96","Or98",
       "Or100P","Or101","Or102" ,"Or105","Or110F","Or111",
       "Or124C","Or128","Or129F","Or132","Or133N","Or144","Or145", "Or146","LOC102656567.b",    
       "Or150", "Or151", "Or154" , "Or163","Or165C","Or166C","Or167PAR","Or169","Or170")
coexp.df<- coexp.df[OR,OR]
pdf("./Fig2E-coexp_OR_chrod_plot.pdf",width=6,height=6)
chorddiag(as.matrix(coexp.df), groupColors = myUmapcolors, groupnamePadding = 10,groupnameFontsize=10, showTicks = FALSE)
dev.off()

# the coexp cluster promoter region 
# the promoter region of coexp_OR 
promoter_region<- read.table("/md01/nieyg/project/honeybee/honebee-latest-Version/04_chemoreceptor/dotplot_OR/promoter_location/OR_promoter_merged_last.bed",sep="\t")
promoter_region$string<- paste(promoter_region$V1,promoter_region$V2,promoter_region$V3,sep="-")
coexp_OR_promoter_region<- promoter_region
donhavepromoter<- setdiff(OR_gene,promoter_region$V4)
for(i in donhavepromoter){
	print(i);
	tmp<- chemoreceptor[chemoreceptor$gene_name==i,];
	tmp<- tmp[1,]
	if(tmp$strand=="+"){
		start<- tmp$start-500;
		end<- tmp$start+500;
		seqnames<- tmp$seqnames;
		tmp_data<- data.frame(V1=seqnames,V2=start,V3=end,V4=i,V5="*",string=paste(seqnames,start,end,sep="-"))
		coexp_OR_promoter_region<- rbind(coexp_OR_promoter_region,tmp_data)
	}
	if(tmp$strand=="-"){
	start<- tmp$end-500;
	end<- tmp$end+500;
	seqnames<- tmp$seqnames;
	tmp_data<- data.frame(V1=seqnames,V2=start,V3=end,V4=i,V5="*",string=paste(seqnames,start,end,sep="-"))
	coexp_OR_promoter_region<- rbind(coexp_OR_promoter_region,tmp_data)
}
}
coexp_OR_promoter_region$V4 <- ORgene_name_trans[match(coexp_OR_promoter_region$V4,ORgene_name_trans$OR_gene),]$last_name
coexp_OR_promoter_region<- na.omit(coexp_OR_promoter_region)
library(GenomicRanges)
# 将字符串转换成 GRanges 对象
granges_obj <- GRanges(
  seqnames = coexp_OR_promoter_region$V1,  # 假设所有序列都来自 "Group15"
  ranges = IRanges(
    start = coexp_OR_promoter_region$V2,  # 提取起始位置
    end = coexp_OR_promoter_region$V3     # 提取结束位置
  )
)
names(granges_obj)<- coexp_OR_promoter_region$V4
DefaultAssay(ORN)<- "ATAC"
macs2_counts <- FeatureMatrix(
     fragments = Fragments(ORN),
     features = granges_obj,
     cells = colnames(ORN)
     )
rownames(macs2_counts)<- coexp_OR_promoter_region$V4


# make the trans dist tree 
object<- ORN
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- rowMeans(x = macs2_counts[, cells])
})
data.dims <- do.call(what = "rbind", args = data.dims)
rownames(data.dims) <- levels(x = object)
dotplot_feature<-unique(rev(as.character(dotplot_data$features.plot)));
dotplot_feature <- ORgene_name_trans[match(dotplot_feature,ORgene_name_trans$OR_gene),]$last_name
order<- dotplot_feature[which(dotplot_feature%in%rownames(macs2_counts))]
data.dims<- data.dims[,order]

singleOR_cluster<-as.character(cluster_info[cluster_info$Freq==1,1])
singleOR_cluster<-singleOR_cluster[-c(1,3,5,14,15,22)]
singleOR_cluster<- c(singleOR_cluster,"57")
singleOR_cluster_OR<- unique(dotplot_data[which(dotplot_data$id%in%singleOR_cluster),]$features.plot)
singleOR_cluster_OR <- ORgene_name_trans[match(singleOR_cluster_OR,ORgene_name_trans$OR_gene),]$last_name
singleOR_cluster_OR <- na.omit(singleOR_cluster_OR);
singleOR_cluster_OR <- singleOR_cluster_OR[-c(3,14)]


data.dims<-data.dims[,which(colnames(data.dims)%in%gsub("\\.","-",colnames(data.dims_last)))]
data.dims<-data.dims[,gsub("\\.","-",colnames(data.dims_last))]

write.csv(data.dims,"OR_promoter_data_tmp.csv")


data.dims_last<- read.csv("OR_promoter_data.csv",row.names=1)
data.dims_last[is.na(data.dims_last)]<- 0
data.dims_last<-data.dims_last[rev(rownames(data.dims_last)),]


##########################################################
pdf("./00_Figure/Fig2/Fig2G-promoter-region-coexp.pdf",width=6,height=4.5)
pheatmap(data.dims_last,
	cluster_cols=FALSE,           
	cluster_rows=FALSE,scale="row",
	color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
    show_colnames = T
    #annotation_col=anno,
    #annotation_colors=an_cols
	)
dev.off()

# Fig2F: the intergenic region length :single vs coexp 
data<- read.csv("./OR_info.csv")
data<- na.omit(data)
data$dist[data$dist>15000]<- 15000;
data$exp_type[data$exp_type=="co-exp_MP"]<- "co-exp"
data$exp_type[data$exp_type=="co-exp_RT"]<- "co-exp"

library(ggpubr)
library(cowplot)
#my_comparisons <- list( c("single_OR", "co-exp_MP"), c("single_OR", "co-exp_RT"), c("co-exp_MP", "co-exp_RT") )
data$exp_type<- factor(data$exp_type,levels=c("single_OR","co-exp"))
pdf("./00_Figure/Fig2/Fig2F_singlevscoexp_Intergenic_region_length_distribution.pdf",width=4,height=4)

ggboxplot(data, x="exp_type", y="dist", 
  bxp.errorbar = T,width=0.6, notch = F)+
stat_compare_means( method='t.test')+theme(legend.position="none")+ylab("Intergenic region length")
dev.off()



















