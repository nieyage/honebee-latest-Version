load("/data/R03/chenxy957/project/Miss.nie/bee/eD_table.RData")
load("/data/R03/chenxy957/project/Miss.nie/bee/cR_table.RData")
eD_table
cR_table
library(Signac)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(sctransform)
library(scater)
library(celda)
library(DropletUtils)
library(Matrix)
library(DoubletFinder)
set.seed(1234)
eD_barcode<- rownames(eD_table)
cR_barcode<- rownames(cR_table)
eD_barcode_unique<- setdiff(eD_barcode,cR_barcode)
cR_barcode_unique<- setdiff(cR_barcode,eD_barcode)
overlap_barcode<- intersect(eD_barcode,cR_barcode)


gtf <- rtracklayer::import('/data/R02/nieyg/ref/10X/Amel_HAv3.1/Amel_HAv3_1/genes/genes.gtf.gz')
gtf$gene_biotype<-"protein_coding"
gene.coords <- gtf[gtf$gene_biotype == 'protein_coding']

# file path 
Forager_counts <- Read10X_h5("/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/raw_feature_bc_matrix.h5")
Forager_fragpath <- "/md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/fragments.tsv.gz"
F_Peaks<-Forager_counts$Peaks;


chrom_assay <- CreateChromatinAssay(
  counts = F_Peaks,
  sep = c(":", "-"),
  fragments = Forager_fragpath,
  annotation = gene.coords,
  min.cells = 12,
  min.features = 200
)

Forager <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

all_barcode<- unique(c(eD_barcode,cR_barcode))


Forager$group<- NA;
data<- data.frame(barcode=c(eD_barcode_unique,cR_barcode_unique,overlap_barcode),
	type=c(rep("eD_unique",length(eD_barcode_unique)),rep("cR_unique",length(cR_barcode_unique)),
		rep("overlap",length(overlap_barcode))))

Forager$group<-data[match(rownames(Forager@meta.data),data$barcode),2]

   Forager <- NucleosomeSignal(Forager)
   Forager <- TSSEnrichment(Forager,fast=FALSE)
Idents(Forager)<- Forager$group
#Forager<- subset(Forager,idents=c("eD_unique","overlap"))
  pdf("./00_Figure/MBE/TSS_distribution_Forager.pdf")
  TSS<-TSSPlot(Forager, group.by = 'group') + NoLegend()+ labs(title = "Forager")
  Frag<-FragmentHistogram(object = Forager, group.by = 'group')+ labs(title = "Forager")
  print(TSS);
  print(Frag);
  dev.off();


TSSPlot <- function(
  object,
  assay = "peaks",
  group.by = NULL,
  idents = NULL
)

assay <- DefaultAssay(object = Forager)
  
  # get the normalized TSS enrichment matrix
  positionEnrichment <- GetAssayData(
    object = Forager,
    assay = assay,
    layer = "positionEnrichment"
  )
  enrichment.matrix <- positionEnrichment[["TSS"]]
object<- Forager
  # average the signal per group per base
  obj.groups <- GetGroups(
    object = object,
    group.by = "group",
    idents = "group"
  )
  groupmeans <- ApplyMatrixByGroup(
    mat = enrichment.matrix,
    groups = obj.groups,
    fun = colMeans,
    normalize = FALSE
  )

  p <- ggplot(
    data = groupmeans,
    mapping = aes(x = position, y = norm.value, color = group)
  ) +
    geom_line(stat = "identity", size = 0.2) +
    facet_wrap(facets = ~group) +
    xlab("Distance from TSS (bp)") +
    ylab(label = "Mean TSS enrichment score") +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank()
    ) +
    ggtitle("TSS enrichment")
  return(p)
}

FragmentHistogram <- function(
  object,
  assay = NULL,
  region = "chr1-1-2000000",
  group.by = NULL,
  cells = NULL,
  log.scale = FALSE,
  ...
) {
  cells <- SetIfNull(x = cells, y = colnames(x = object))
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))

  reads <- MultiGetReadsInRegion(
    object = object,
    assay = assay,
    region = region,
    cells = cells,
    verbose = FALSE,
    ...
  )
  # add group information
  if (is.null(x = group.by)) {
    groups <- Idents(object = object)
  } else {
    md <- object[[]]
    groups <- object[[group.by]]
    groups <- groups[, 1]
    names(x = groups) <- rownames(x = md)
  }
  reads$group <- groups[reads$cell]
  if (length(x = unique(x = reads$group)) == 1) {
    p <- ggplot(data = reads, aes(length)) +
      geom_histogram(bins = 200)
  } else {
    p <- ggplot(data = reads, mapping = aes(x = length, fill = group)) +
      geom_histogram(bins = 200) +
      facet_wrap(~group, scales = "free_y")
  }
  p <- p + xlim(c(0, 800)) +
    theme_classic() +
    theme(
      legend.position = "none",
      strip.background = element_blank()
    ) +
    xlab("Fragment length (bp)") +
    ylab("Count")
  if (log.scale) {
    p <- p + scale_y_log10()
  }
  return(p)
}


