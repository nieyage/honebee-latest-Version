# Step1: Seurat 2 loom and scanpy (scRNA)
library(loomR)
# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
  library(SCopeLoomR)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)
  library(patchwork)
  set.seed(1234)
})
library(SeuratDisk)
ORN<- readRDS("./05_ORN_cluster2/02_second_cluster/06_rm_without_power/Unsupervised_ORN_remove_nopower_modify_the_tsne_recall_peak.rds")
DefaultAssay(ORN)<-"raw_RNA"
# seurat对象转换为loom文件
sdata.loom <- as.loom(x = ORN, filename = "/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus/01_scRNA/ORN_seurat2loom.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()

# sdata 2 scanpy
# in python
import scanpy as sc
adata = sc.read_loom("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus/01_scRNA/ORN_seurat2loom.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
adata.write('/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus/01_scRNA/adata.h5ad', compression='gzip')

# Step2: Signac 2 pycisTopic 
conda activate scenicplus
#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = '/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus/'
import pycisTopic
#make a directory for to store the processed scRNA-seq data.
#if not os.path.exists(os.path.join(work_dir, '02_scATAC')):
#    os.makedirs(os.path.join(work_dir, '02_scATAC'))
tmp_dir = '/md01/nieyg/tmp'

# in shell
# merge the fragment file 
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/fragments.tsv .
mv fragments.tsv Forager_fragments.tsv
cp /md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/fragments.tsv .
mv fragments.tsv NE_fragments.tsv
cp /md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/fragments.tsv .
mv fragments.tsv Nurse_fragments.tsv
# add sample info in barcode 
awk '$4="NE_"$4' NE_fragments.tsv > NE_fragments_2.tsv
awk '$4="Nurse_"$4' Nurse_fragments.tsv > Nurse_fragments_2.tsv
awk '$4="Forager_"$4' Forager_fragments.tsv > Forager_fragments_2.tsv
cat NE_fragments_2.tsv Nurse_fragments_2.tsv Forager_fragments_2.tsv > last_fragments.tsv
awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5}' last_fragments.tsv > last_fragments_2.tsv
cat last_fragments_2.tsv |sort -k1,1 -k2,2n | bgzip > last_fragments.tsv.gz
tabix -b 2 -e 3 -p bed last_fragments.tsv.gz

# in python 
fragments_dict = {'ORN': os.path.join(work_dir, '02_scATAC/fragments_file/last_fragments.tsv.gz')}

# 1.generate pseudobulk ATAC-seq profiles per cell type and call peaks

# load the cell type annotation we generated in the scRNA-seq analysis above.
import scanpy as sc
adata = sc.read_h5ad(os.path.join(work_dir, '01_scRNA/adata.h5ad'))
cell_data = adata.obs
cell_data['sample_id'] = 'ORN'
cell_data['subcluster'] = cell_data['subcluster'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.

# Get chromosome sizes (for hg38 here)
import pyranges as pr
import requests
import pandas as pd
chromsizes=pd.read_csv("/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt", sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                 variable = 'subcluster',                                                                     # variable by which to generate pseubulk profiles, in this case we want pseudobulks per celltype
                 sample_id_col = 'sample_id',
                 chromsizes = chromsizes,
                 bed_path = os.path.join(work_dir, '02_scATAC/consensus_peak_calling/pseudobulk_bed_files/'),  # specify where pseudobulk_bed_files should be stored
                 bigwig_path = os.path.join(work_dir, '02_scATAC/consensus_peak_calling/pseudobulk_bw_files/'),# specify where pseudobulk_bw_files should be stored
                 path_to_fragments = fragments_dict,                                                        # location of fragment fiels
                 n_cpu = 8,                                                                                 # specify the number of cores to use, we use ray for multi processing
                 normalize_bigwig = True,
                 remove_duplicates = True,
                 _temp_dir = os.path.join("/md01/nieyg/tmp", 'ray_spill'),
                 split_pattern = '-')

# save the bw and bed files 
import pickle
pickle.dump(bed_paths,
            open(os.path.join(work_dir, '02_scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'wb'))
pickle.dump(bw_paths,
           open(os.path.join(work_dir, '02_scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'wb'))

# Call peaks per pseudobulk profile
import pickle
bed_paths = pickle.load(open(os.path.join(work_dir, '02_scATAC/consensus_peak_calling/pseudobulk_bed_files/bed_paths.pkl'), 'rb'))
bw_paths =  pickle.load(open(os.path.join(work_dir, '02_scATAC/consensus_peak_calling/pseudobulk_bed_files/bw_paths.pkl'), 'rb'))
from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path='/data/R02/nieyg/ori/biosoft/conda/bin/macs2'
# Run peak calling
narrow_peaks_dict = peak_calling(macs_path,
                                 bed_paths,
                                 os.path.join(work_dir, '02_scATAC/consensus_peak_calling/MACS/'),
                                 genome_size='2.7e+09',
                                 n_cpu=8,
                                 input_format='BEDPE',
                                 shift=73,
                                 ext_size=146,
                                 keep_dup = 'all',
                                 q_value = 0.05,
                                 _temp_dir = os.path.join("/md01/nieyg/tmp", 'ray_spill'))

pickle.dump(narrow_peaks_dict,
            open(os.path.join(work_dir, '02_scATAC/consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'wb'))


# 2.merge these peaks into a consensus peak-set

from pycisTopic.iterative_peak_calling import *
# Other param
peak_half_width = 250
#path_to_blacklist= '../../pycisTopic/blacklist/hg38-blacklist.v2.bed'
# Get consensus peaks
consensus_peaks=get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes)
consensus_peaks.to_bed(
    path = os.path.join(work_dir, '02_scATAC/consensus_peak_calling/consensus_regions.bed'),
    keep=True,
    compression='infer',
    chain=False)


# 3.run topic modeling to find sets of co-accessible regions and to impute chromatin accessibility resolving the issue of drop outs

import scanpy as sc
adata = sc.read_h5ad(os.path.join(work_dir, '01_scRNA/adata.h5ad'))
scRNA_bc = adata.obs_names
cell_data = adata.obs
cell_data['sample_id'] = 'ORN'
cell_data['subcluster'] = cell_data['subcluster'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
del(adata)

import pickle

path_to_regions = {'ORN':os.path.join(work_dir, '02_scATAC/consensus_peak_calling/consensus_regions.bed')}
#path_to_blacklist= '../../pycisTopic/blacklist/hg38-blacklist.v2.bed'
metadata_bc = pickle.load(open(os.path.join(work_dir, '02_scATAC/quality_control/metadata_bc.pkl'), 'rb'))

from pycisTopic.cistopic_class import *
key = 'ORN'
cistopic_obj = create_cistopic_object_from_fragments(
                            path_to_fragments=fragments_dict[key],
                            path_to_regions=path_to_regions[key],
                            #path_to_blacklist=path_to_blacklist,
                            #metrics=metadata_bc[key],
                            valid_bc=list(set(scRNA_bc)),
                            n_cpu=1,
                            project=key,
                            split_pattern='-')
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
print(cistopic_obj)

pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, '02_scATAC/cistopic_obj.pkl'), 'wb'))

# Run topic modeling. The purpose of this is twofold:
#To find sets of co-accessible regions (topics), this will be used downstream as candidate enhancers (together with Differentially Accessible Regions (DARs)).
#To impute dropouts.

# Topic modeling can be computationaly intense!
import pickle
cistopic_obj = pickle.load(open(os.path.join(work_dir, '02_scATAC/cistopic_obj.pkl'), 'rb'))
from pycisTopic.cistopic_class import *
models=run_cgs_models(cistopic_obj,
                    n_topics=[2,4,10,16,32,48],
                    n_cpu=5,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = os.path.join(tmp_dir + 'ray_spill'))


if not os.path.exists(os.path.join(work_dir, '02_scATAC/models')):
    os.makedirs(os.path.join(work_dir, '02_scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, '02_scATAC/models/10x_ORN_models_500_iter_LDA.pkl'), 'wb'))

models = pickle.load(open(os.path.join(work_dir, '02_scATAC/models/10x_ORN_models_500_iter_LDA.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, '02_scATAC/cistopic_obj.pkl'), 'rb'))
from pycisTopic.lda_models import *
fig, ax = plt.subplots(figsize = (12,12))
import numpy as np
import matplotlib.pyplot as plt

model = evaluate_models(models,
                       select_model=16,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)
plt.savefig('evaluate_models.png')
plt.show()

#The metrics seem to stabelise with a model using 16 topics, so let’s choose that model.

cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, '02_scATAC/cistopic_obj.pkl'), 'wb'))


# Visualization
# umap:cell-topic probabilities to generate dimensionality reductions
from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)
plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['subcluster'])
plt.savefig('02_scATAC/metadata_umap.pdf')
plt.show()
plot_topic(cistopic_obj, reduction_name = 'UMAP', num_columns = 4)
plt.savefig('02_scATAC/topics_top_16_umap.pdf')
plt.show()

# Inferring candidate enhancer regions
# 1.First we will binarize the topics using the otsu method and by taking the top 3k regions per topic.

from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

# 2. Next we will calculate DARs per cell type
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='subcluster', var_features=variable_regions, split_pattern = '_',adjpval_thr=0.5,log2fc_thr=0.1,n_cpu=10)

if not os.path.exists(os.path.join(work_dir, '02_scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, '02_scATAC/candidate_enhancers'))
import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, '02_scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, '02_scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, '02_scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))

# Step3: Motif enrichment analysis using pycistarget

import pickle
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, '02_scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, '02_scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, '02_scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}

for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('Group')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('Group')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))

for DAR in ['1', '10', '11',  '13',  '16', '17',  '18',  '2', '23', '25', '27', '28', '3',
  '34', '35', '36', '37', '38', '39', '4', '40',  '42',  '43', '44',
  '46', '47', '49', '5', '51', '54', '55', '56',  '58',  '6', '60', '8', '9']:
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('Group')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    if len(regions) > 0:
        region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))


for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')


# in R and shell 
# create the custom annotation file 
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf')
gtf<- gtf[gtf$type=="transcript",]
gtf<- as.data.frame(gtf)
data<- data.frame(Chromosome=gtf$seqnames,Start=gtf$start,Strand=gtf$strand,Gene=gtf$gene_name,Transcript_type="protein_coding")
write.table(data,"/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_annot.tbl",row.names=F,sep="\t")
sed -i 's/"//g' custom_annot.tbl
sed -i '/MT/d;' custom_annot.tbl
sed -i 's/+/1/g' custom_annot.tbl
sed -i 's/-/-1/g' custom_annot.tbl


if not os.path.exists(os.path.join(work_dir, '03_motifs')):
    os.makedirs(os.path.join(work_dir, '03_motifs'))

from scenicplus.wrappers.run_pycistarget import run_pycistarget

from ctxcore.rnkdb import FeatherRankingDatabase
from pycistarget.utils import target_to_query

custom_annot =pd.read_csv("/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_annot.tbl", sep='\t')
our_rankings_db = os.path.join('/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/03_region/01_honeybee_cisBp.regions_vs_motifs.rankings.feather')
our_scores_db = os.path.join('/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/03_region/01_honeybee_cisBp.regions_vs_motifs.scores.feather')
our_motif_annotation = os.path.join('/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/honeybee_cisbp_motif2TF.tbl')
run_pycistarget(
    region_sets = region_sets,
    species = 'custom',
    custom_annot=custom_annot,
    save_path = os.path.join(work_dir, '03_motifs'),
    ctx_db_path = our_rankings_db,
    dem_db_path = our_scores_db,
    promoter_space=20,
    path_to_motif_annotations = our_motif_annotation,
    run_without_promoters = True,
    n_cpu = 1,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10',
    )

import dill
menr = dill.load(open(os.path.join(work_dir, '03_motifs/menr_honeybee_cisBp_motif.pkl'), 'rb'))

pb_fly_ranking_db =  os.path.join('/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/03_region/01_honeybee_cisBp.regions_vs_motifs.rankings.feather')
pb_fly_scores_db =  os.path.join('/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/03_region/01_honeybee_cisBp.regions_vs_motifs.rankings.feather')

pb_fly_motif_annotation = os.path.join('/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/motifs-v10nr_clust-nr.flybase-m0.001-o0.0-trans2honeybee.tbl')

run_pycistarget(
    region_sets = region_sets,
    species = 'custom',
    custom_annot=custom_annot,
    save_path = os.path.join(work_dir, '03_motifs'),
    ctx_db_path = pb_fly_ranking_db,
    dem_db_path = pb_fly_scores_db,
    path_to_motif_annotations = pb_fly_motif_annotation,
    run_without_promoters = True,
    n_cpu = 1,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust',
    )

import dill
menr = dill.load(open(os.path.join(work_dir, '03_motifs/menr_database_fly.pkl'), 'rb'))


# We now have completed all the steps necessary for starting the SCENIC+ analysis
# Step4: inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = '/data/R02/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus'

adata = sc.read_h5ad(os.path.join(work_dir, '01_scRNA/adata.h5ad'))
cistopic_obj = dill.load(open(os.path.join(work_dir, '02_scATAC/cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, '03_motifs/menr_honeybee_cisBp_motif.pkl'), 'rb'))

from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr,
    bc_transform_func = lambda x: f'{x}-ORN' #function to convert scATAC-seq barcodes to scRNA-seq ones
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
scplus_obj


#only keep the first two columns of the PCA embedding in order to be able to visualize this in SCope
#scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]
#scplus_obj.dr_cell['GEX_rep'] = scplus_obj.dr_cell['GEX_rep'].iloc[:, 0:2]

motifs_list_filename="/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/01_creat_cisTarget_databases/01_we_collect_motifs.lst"

# skip biomart_host 
# https://github.com/aertslab/scenicplus/issues/48

# Get chromosome sizes (for hg38 here)
import pyranges as pr
import requests
import pandas as pd
chromsizes=pd.read_csv("/md01/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/star/chrNameLength.txt", sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)

# create the pr_annot 
library(Biostrings)
gtf <- rtracklayer::import('/md01/nieyg/ref/10X/Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic_UPDATED_chemoreceptors_CORRECTED_geneid.sorted.gtf')
gtf<- gtf[gtf$type=="transcript",]
gtf<- as.data.frame(gtf)

data<- data.frame(Chromosome=gtf$seqnames,Start=gtf$start,End=gtf$end,
	Strand=gtf$strand,Gene=gtf$gene_name,
	Transcription_Start_Site=gtf$start,Transcript_type="protein_coding")
data$Strand <- as.character(data$Strand)
data$Chromosome <- as.character(data$Chromosome)

write.table(data,"/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_pr_annot.tbl",row.names=F,sep="\t")
sed -i 's/"//g' /md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_pr_annot.tbl
sed -i '/MT/d;' /md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_pr_annot.tbl
#sed -i 's/+/1/g' /md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_pr_annot.tbl
#sed -i 's/-/-1/g' /md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_pr_annot.tbl

pr_annot=pd.read_csv("/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/custom_pr_annot.tbl", sep='\t')
pr_annot=pr.PyRanges(pr_annot)

from scenicplus.enhancer_to_gene import get_search_space
get_search_space(
        scplus_obj,
        species = None,
        assembly = None,
        pr_annot = pr_annot,
        pr_chromsizes = chromsizes,
        upstream = [500, 150000],
        downstream = [500, 150000])

from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_subcluster'],
        species = None,
        assembly = None,
        tf_file = motifs_list_filename,
        save_path = os.path.join(work_dir, '04_scenicplus'),
        biomart_host = 'http://www.ensembl.org',
        upstream = [500, 10000],
        downstream = [500, 10000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = '/md01/nieyg/software/bedToBigBed',
        n_cpu = 12,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, '04_scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)
import pickle
pickle.dump(scplus_obj, open(os.path.join(work_dir, '04_scenicplus/scplus_obj_default_honeybee_cisBPmotif.pkl'), 'wb'))

# save html file 

with open('output.html', 'w') as f:
   f.write(menr['DEM_topics_otsu_All']).DEM_results('Topic8').data)


# Step5: Note on the output of SCENIC+

scplus_obj

scplus_obj.to_df('EXP').head()

scplus_obj.to_df('ACC').head()

scplus_obj.metadata_cell.head()
scplus_obj.metadata_regions.head()
scplus_obj.metadata_genes.head()
scplus_obj.menr.keys()
scplus_obj.dr_cell.keys()
scplus_obj.uns.keys()

scplus_obj.uns['eRegulons'][0:5]
for attr in dir(scplus_obj.uns['eRegulons'][0]):
    if not attr.startswith('_'):
        print(f"{attr}: {getattr(scplus_obj.uns['eRegulons'][0], attr) if not type(getattr(scplus_obj.uns['eRegulons'][0], attr)) == list else getattr(scplus_obj.uns['eRegulons'][0], attr)[0:5]}")
scplus_obj.uns['eRegulon_metadata'].head()

# Step6: Downstream analysis

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
import dill
work_dir = 'pbmc_tutorial'
scplus_obj = dill.load(open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'rb'))

from scenicplus.preprocessing.filtering import apply_std_filtering_to_eRegulons
apply_std_filtering_to_eRegulons(scplus_obj)

scplus_obj.uns['eRegulon_metadata_filtered'].head()

from scenicplus.eregulon_enrichment import score_eRegulons
region_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/region_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
gene_ranking = dill.load(open(os.path.join(work_dir, 'scenicplus/gene_ranking.pkl'), 'rb')) #load ranking calculated using the wrapper function
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 5)
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures_filtered',
                key_added = 'eRegulon_AUC_filtered',
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 5)
from scenicplus.dimensionality_reduction import run_eRegulons_tsne, run_eRegulons_umap
run_eRegulons_umap(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_UMAP', #overwrite previously calculated UMAP
)
run_eRegulons_tsne(
    scplus_obj = scplus_obj,
    auc_key = 'eRegulon_AUC_filtered',
    reduction_name = 'eRegulons_tSNE', #overwrite previously calculated tSNE
)

from scenicplus.dimensionality_reduction import plot_metadata_given_ax
import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline

#specify color_dictionary

color_dict = {
    'B_cells_1': "#065143",
    'B_cells_2': "#70B77E",
    'CD4_T_cells': "#E0A890",
    'CD8_T_cells': "#F56476",
    'NK_cells': "#CE1483",
    'Dendritic_cells': "#053C5E" ,
    'FCGR3A+_Monocytes': "#38A3A5",
    'CD14+_Monocytes': "#80ED99"
}

fig, axs = plt.subplots(ncols=2, figsize = (16, 8))
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[0],
    reduction_name = 'eRegulons_UMAP',
    variable = 'GEX_celltype', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    color_dictionary={'GEX_celltype': color_dict}
)
plot_metadata_given_ax(
    scplus_obj=scplus_obj,
    ax = axs[1],
    reduction_name = 'eRegulons_tSNE',
    variable = 'GEX_celltype', #note the GEX_ prefix, this metadata originated from the gene expression metadata (on which we did the cell type annotation before)
    color_dictionary={'GEX_celltype': color_dict}
)
fig.tight_layout()
sns.despine(ax = axs[0]) #remove top and right edge of axis border
sns.despine(ax = axs[1]) #remove top and right edge of axis border
plt.show()

from scenicplus.dimensionality_reduction import plot_eRegulon
plot_eRegulon(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    selected_regulons = ['EOMES_+', 'GATA3_+', 'TCF7_+', 'CEBPA_+', 'PAX5_+'],
    scale = True,
    auc_key = 'eRegulon_AUC_filtered')

from scenicplus.dimensionality_reduction import plot_AUC_given_ax

fig, ax = plt.subplots(figsize = (8,8))
plot_AUC_given_ax(
    scplus_obj = scplus_obj,
    reduction_name = 'eRegulons_tSNE',
    feature = 'PAX5_+_(119g)',
    ax = ax,
    auc_key = 'eRegulon_AUC_filtered',
    signature_key = 'Gene_based')
sns.despine(ax = ax)
plt.show()

from scenicplus.cistromes import TF_cistrome_correlation, generate_pseudobulks

generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Gene_based')
generate_pseudobulks(
        scplus_obj = scplus_obj,
        variable = 'GEX_celltype',
        auc_key = 'eRegulon_AUC_filtered',
        signature_key = 'Region_based')

TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Gene_based',
            out_key = 'filtered_gene_based')
TF_cistrome_correlation(
            scplus_obj,
            use_pseudobulk = True,
            variable = 'GEX_celltype',
            auc_key = 'eRegulon_AUC_filtered',
            signature_key = 'Region_based',
            out_key = 'filtered_region_based')

scplus_obj.uns['TF_cistrome_correlation']['filtered_region_based'].head()








































