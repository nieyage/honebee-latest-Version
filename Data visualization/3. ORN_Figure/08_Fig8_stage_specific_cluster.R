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

## stage specific cluster
ORN$orig.ident<-factor(ORN$orig.ident,levels=c("NE","Nurse","Forager"))
data<-as.data.frame(table(ORN$orig.ident,ORN$subcluster))
ORN_sample_number<-as.data.frame(table(ORN$orig.ident))
NE_all_number<-ORN_sample_number[which(ORN_sample_number$Var1=="NE"),2]
Nurse_all_number<-ORN_sample_number[which(ORN_sample_number$Var1=="Nurse"),2]
Forager_all_number<-ORN_sample_number[which(ORN_sample_number$Var1=="Forager"),2]
total<-NE_all_number+Nurse_all_number+Forager_all_number
#  chisq.test;
#  fisher.test more sensitive 
fisher<-data.frame()
for(i in unique(data$Var2)){
    Clu<-data[data$Var2==i,];
    Clu<-Clu$Freq
    all<-c(NE_all_number,Nurse_all_number,Forager_all_number)
    df <- as.data.frame(rbind(Clu,all))
    colnames(df) <- c("NE","Nurse","Forager")
    print(i)
    #fisher.test(df)
    #chisq.test(df)
    p_value<-fisher.test(df)$p.value;
    #p_value<-chisq.test(df)$p.value;
    #X_squared<-fisher.test(df)$statistic;
    fisher[1,i]<-p_value;
    #chisq[2,i]<-X_squared;
}
rownames(fisher)<-c("p_value")
Idents(ORN)<-ORN$subcluster
#  sample propotion in each cluster;
library(dittoSeq)
p<-dittoBarPlot(ORN, "orig.ident", group.by = "subcluster",data.out = TRUE)
# change levels
p$data$grouping<-factor(p$data$grouping,levels=levels(factor(ORN$subcluster)))
p$data$label<-factor(p$data$label,levels=c("NE","Nurse","Forager"))
por_data<-p$data;
por_data$grouping<-as.character(por_data$grouping)
por_data<-por_data[por_data$grouping%in%colnames(fisher),]
por_data$p_value<-100
for(i in 1:length(por_data$grouping)){
  df<-fisher[,por_data$grouping[i]];
  por_data$p_value[i]=df[1]
  }
por_data$signif<-" "
for(i in 1:length(por_data$grouping)){
  if(por_data$p_value[i]<0.05){
    por_data$signif[i]="*"
  }
  };
por_data$grouping<-factor(por_data$grouping,levels=levels(Idents(ORN)))
one_class_por_data<-por_data
pdf("./05_ORN_cluster2/06_stage_specific/stage_specific_cluster_proportion.pdf",height=2,width=12)
ggplot(data = por_data, aes_string(x = "grouping", y = "percent", 
        fill = "label")) +  xlab("orig.ident") + ylab("Percent of cells") + 
        scale_fill_manual(values = c("#5CC8F2","#009E73","#E69F00")) + 
        geom_col()+ scale_y_continuous(breaks = c(0, Forager_all_number/total,1-NE_all_number/total, 1))+ coord_cartesian(ylim = c(0,1))+
geom_hline(aes(yintercept=Forager_all_number/total),linetype="dashed",col="black")+
geom_hline(aes(yintercept=1-NE_all_number/total),linetype="dashed",col="black")+
geom_text(aes(label=signif,y =0))+theme(axis.text.x = element_text(angle = -90,vjust = 0.5,hjust = 0.5))
dev.off()

## Stage specific 
signif<-por_data[por_data$p_value<0.05,]
global_ratio<-c(NE_all_number/total,Nurse_all_number/total,Forager_all_number/total)
signif$global_ratio<-c(rep(global_ratio,nrow(signif)/3))
signif$FC<-signif$percent/signif$global_ratio;

#all_class_Data<-rbind(one_class_data,multiple_classes_data)
cluster<-unique(signif$grouping)
#cluster<-cluster[-12]
all_class_Data<-dotplot_data
cluster_specific<-data.frame()
for (i in cluster){
  signif_subset<-signif[which(signif$grouping==i),];
  sample_specific<-signif_subset[which(signif_subset$FC==max(signif_subset$FC)),]$label
  OR_specfic<-all_class_Data[which(all_class_Data$id==i),]$features.plot;
  if(length(OR_specfic)!=1){
    OR_specfic<-paste(OR_specfic,collapse=",")
  }
  cluster_specific_subset<-data.frame(cluster=i,OR_specfic,sample_specific)
  cluster_specific<-rbind(cluster_specific,cluster_specific_subset)
} 
write.csv(cluster_specific,"./05_ORN_cluster2/06_stage_specific/cluster_specific.csv")

