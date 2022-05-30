###########################
#Prepare enviroment
set.seed(1234567)
options(stringsAsFactors = FALSE)

#DLLs controls and changed for >100 DLLs :)
Sys.getenv()
length(getLoadedDLLs()) 

#Set working dir
setwd("/R_functions_path")

###########################
#LOAD -> Libraries
###########################
source("seurat_melanoma_timeseries_libraries.R")
###########################
#LOAD -> FUNCTION
###########################
source("seurat_melanoma_timeseries_functions.R")
###########################
#LOAD -> CHARTS
###########################
source("seurat_melanoma_timeseries_charts.R")


###########################
#Working directory
setwd("/your_working_directory")

load("filtered_rna_atac_data.RData")

#######################################################
#PROCESSING -> COMBINE DATASETS SEURAT CCA RNAseq + ATACseq
#######################################################
#As symplify version of Seurat propouse
#Generate the Gene Activity Matrix from ATACseq data / seq.levels tune for the number of chromosomes of the specie (mouse=19, human=21)
activity.matrix = CreateGeneActivityMatrix(peak.matrix=count_atac,
                                            annotation.file="/genomes/genes.gtf", 
                                            seq.levels=c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE
                                            )


#######################################################
#Creating Seurat objects
#######################################################

#######################################################
#ATACseq
pbmc.atac = CreateSeuratObject(
  counts = count_atac, assay = "ATAC", project = "10x_ATAC"
)

#ADD the GeneActivityMatrix and process it
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
pbmc.atac$tech <- "atac"

rm(list=setdiff(ls(), c("pbmc.atac","count_rna","seurat_data_plot","seurat_data_plot_ATAC","seurat_data_plot_RNA","blankPlot")))

#######################################################
#Adding % mit and % rib information for each cell
pbmc.atac[["percent.mt"]] <- PercentageFeatureSet(pbmc.atac, pattern = "^MT-")
pbmc.atac[["percent.rb"]] <- PercentageFeatureSet(pbmc.atac, pattern = c("^RPS", "^RPL"))

#######################################################
#Visualize data for filtering propouse
seurat_data_plot_ATAC(paste0("10x_ATAC_Seurat_data_filtered.pdf"))

#REMOVE CELLS HAVING MORE THAN 5% MIT GENES
pbmc.atac = subset(x = pbmc.atac, percent.mt < 5)

#######################################################
#Visualize data for filtering propouse
seurat_data_plot_ATAC(paste0("10x_ATAC_Seurat_data_filtered_2.pdf"))

DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac)
pbmc.atac <- FindVariableFeatures(pbmc.atac)

#Process the ATAC data
DefaultAssay(pbmc.atac) <- "ATAC"
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q0') 
pbmc.atac <- RunSVD(
  object = pbmc.atac,
  assay = 'ATAC',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

pbmc.atac = RunUMAP(pbmc.atac, reduction = "lsi", dims = 1:30)


#######################################################
#RNAseq
pbmc.rna = CreateSeuratObject(
  counts = count_rna, min.cells = 3, min.features = 200
)

pbmc.rna$tech = "rna"
rm(list=setdiff(ls(), c("pbmc.atac","pbmc.rna","seurat_data_plot","seurat_data_plot_ATAC","seurat_data_plot_RNA","blankPlot")))

#######################################################
#Adding % mit and % rib information for each cell
pbmc.rna[["percent.mt"]] <- PercentageFeatureSet(pbmc.rna, pattern = "^MT-")
pbmc.rna[["percent.rb"]] <- PercentageFeatureSet(pbmc.rna, pattern = c("^rpRPSs", "^RPL"))

#######################################################
#Visualize data for filtering propouse
seurat_data_plot_RNA(paste0("10x_RNA_Seurat_data_filtered.pdf"))

#REMOVE CELLS HAVING MORE THAN 5% MIT GENES
pbmc.rna = subset(x = pbmc.rna, percent.mt < 5)

#######################################################
#Visualize data for filtering propouse
seurat_data_plot_RNA(paste0("10x_RNA_Seurat_data_filtered_2.pdf"))

pbmc.rna = NormalizeData(pbmc.rna)
pbmc.rna = ScaleData(pbmc.rna)
pbmc.rna = FindVariableFeatures(pbmc.rna)
pbmc.rna = RunPCA(object = pbmc.rna, npcs = 30,verbose = FALSE)
pbmc.rna = RunUMAP(object = pbmc.rna, reduction = "pca",
                          dims = 1:15)

#######################################################
#Plot resulting analysis
pdf(file="RNA_and_ATAC_post_analysis.pdf", width=12, height=6)
p1 = DimPlot(object = pbmc.atac, reduction = "umap", label = TRUE) + ggtitle("scATAC-seq")
p2 = DimPlot(object = pbmc.rna, reduction = "umap", label = TRUE) + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))
dev.off()

#######################################################
#Transfer and integrate RNAseq and ATACseq using RNAseq as reference
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
    
genes.use = VariableFeatures(pbmc.rna)

refdata = GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]

imputation = TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]])

pbmc.atac[["RNA"]] = imputation
coembed = merge(x = pbmc.rna, y = pbmc.atac)

coembed = ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed = RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed = RunUMAP(coembed, dims = 1:30)

coembed = FindNeighbors(object = coembed, dims = 1:15)
coembed = FindClusters(object = coembed, reduction.type = "pca",dims = 1:15, resolution = 0.6)

pdf(file="RNA_and_ATAC_UMAP_2.pdf", width=12, height=6)
p1 = DimPlot(coembed, reduction = "umap", group.by = "tech")
p2 = DimPlot(coembed, reduction = "umap", label = TRUE)
CombinePlots(list(p1, p2))
dev.off()


#########################################################################################
#########################################################################################
#FINDMARKERS all clusters agains rest of cells
DefaultAssay(coembed) <- "RNA"

list_RNAseq_data_raw.markers = FindAllMarkers(object = coembed, only.pos = TRUE, min.pct = 0.25)
print(x = head(x = list_RNAseq_data_raw.markers, n = 5))
list_RNAseq_data_raw.markers = (subset(list_RNAseq_data_raw.markers,(list_RNAseq_data_raw.markers['p_val_adj'] <= 0.01 )))

write.csv(list_RNAseq_data_raw.markers, file = "Markers.csv")

#Plot markers for each cluster
library(dplyr)
list_RNAseq_data_raw.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10

pdf(file="markers.pdf", width=16, height=20)
DoHeatmap(object = coembed, features = top10$gene) + NoLegend()
dev.off()


save.image("your_standard_analysis_Seurat_output.RData")

