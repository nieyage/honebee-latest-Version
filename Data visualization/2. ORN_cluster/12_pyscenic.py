# make exp matrix 
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
exprMat <- as.matrix(ORN@assays$raw_RNA@counts)
cellInfo <- data.frame(CellType=Idents(ORN))
write.csv(exprMat,"./09_GRN/02_SCENIC/ORN_exprMat.csv")
write.csv(cellInfo,"./09_GRN/02_SCENIC/ORN_cellInfo.csv")

# make the motif2TF.tbl file of honeybee motif 
grep ">" CisBP-honeybee.jaspar|cut -d ">" -f 2 > honeybee_motif_TF.list

data<- read.csv("honeybee_motif_TF.list",sep=" ")
data$motif_name<- data$motif_id
data$source_name <- "CisBP"
data$source_version <-"None"
data$gene_name <- data$motif_description
data$motif_similarity_qvalue <- 0.000000
data$similar_motif_id <- "None"
data$similar_motif_description<- "None"
data$orthologous_identity  <- 1.000000 
data$orthologous_gene_name <- "None"
data$orthologous_species <-"None"
data$description<- "gene is directly annotated"
data<- data[,c(1,3,2,4:13)]
write.table(data,"honeybee_cisbp_motif2TF.tbl",row.names=F,sep="\t")

sed -i 's/"//g' honebee_cisbp_motif2TF.tbl

#motif_id motif_name  motif_description source_name source_version  gene_name motif_similarity_qvalue similar_motif_id  similar_motif_description orthologous_identity  orthologous_gene_name orthologous_species description
jaspar__MA0002.2  MA0002.2  RUNX1 jaspar  2016  RUNX1 0.000000  None  None  1.000000  None  None  gene is directly annotated
#/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/honeybee_cisbp_motif2TF.tbl

# Step1: co-expression modules



import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

DATA_FOLDER="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/03_pyscenic/"

DATABASES_GLOB = os.path.join("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/up500_trans_cisbp.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/honeybee_cisbp_motif2TF.tbl")
MM_TFS_FNAME = os.path.join("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/TF_name.list")
SC_EXP_FNAME = os.path.join("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/ORN_exprMat.csv")

MODULES_FNAME = os.path.join(DATA_FOLDER, "modules.p")
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
ADJACENCIES_FNAME = os.path.join(DATA_FOLDER, "adjacencies.csv")
ex_matrix = pd.read_csv(SC_EXP_FNAME, sep=',', header=0, index_col=0).T
ex_matrix.head()

tf_names = load_tf_names(MM_TFS_FNAME)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

# Phase I: Inference of co-expression modules
adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=True)
adjancencies.head()
adjancencies.to_csv(ADJACENCIES_FNAME, index=False, sep='\t')


modules = list(modules_from_adjacencies(adjancencies, ex_matrix))
with open(MODULES_FNAME, 'wb') as f:
    pickle.dump(modules, f)

df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

df.head()

df.to_csv(MOTIFS_FNAME)

regulons = df2regulons(df)

with open(REGULONS_FNAME, 'wb') as f:
    pickle.dump(regulons, f)

#with open(REGULONS_FNAME, 'rb') as f:
#    regulons = pickle.load(f)

auc_mtx = aucell(ex_matrix, regulons, num_workers=1)

sns.clustermap(auc_mtx, figsize=(12,12))











