# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(SCopeLoomR)
  library(KernSmooth)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(ComplexHeatmap)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)
  library(patchwork)
  set.seed(1234)
})


ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
DefaultAssay(ORN)<-"raw_RNA"

## Get data from sce object:
exprMat <- ORN@assays$raw_RNA
cellInfo <- data.frame(CellType=Idents(ORN))

## Color to assign to the variables (same format as for NMF::aheatmap)
#colVars <- list(CellType=c("microglia"="forestgreen", 
#                           "endothelial-mural"="darkorange", 
#                           "astrocytes_ependymal"="magenta4", 
#                           "oligodendrocytes"="hotpink", 
#                           "interneurons"="red3", 
#                           "pyramidal CA1"="skyblue", 
#                           "pyramidal SS"="darkblue"))
#colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
#saveRDS(colVars, file="int/colVars.Rds")
#plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


### Initialize settings
library(RcisTarget)
library(feather)
dfile <- "/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/up500_trans_cisbp.genes_vs_motifs.rankings.feather"
db <- importRankings(dfile, indexCol = "motifs")       
#(setting 'indexCol = "motifs"' makes sure the 'motifs' column is made the first column in the tibble)
db@rankings[1:5,1:5]
names(db@rankings)[1] <- "features"
write_feather(db@rankings, "/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/up500_trans_cisbp.genes_vs_motifs.rankings2.feather")

scenicOptions <- initializeScenic(org="hgnc", 
                                  nCores=2,
                                  dbDir=mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "GOAT")


scenicOptions<- initializeScenic(org = "dmel", 
	dbDir = "/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/", 
	dbs ="up500_trans_cisbp.genes_vs_motifs.rankings2.feather", nCores = 12)

# Step1: Co-expression network

# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)
# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)

# Build and score the GRN (runSCENIC_…)
library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save stat







