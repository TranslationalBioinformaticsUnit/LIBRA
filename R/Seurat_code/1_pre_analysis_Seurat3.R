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
source("libraries.R")
###########################
#LOAD -> FUNCTION
###########################
source("functions.R")
###########################
#LOAD -> CHARTS
###########################
source("charts.R")


###########################
#LOAD DATA__ 2 (Dataset_1)
###########################
#RNA
indata <- Matrix::readMM("GSE126074_P0_BrainCortex_SNAREseq_cDNA.counts.mtx") 
# format cell info
cellinfo <- read.table("GSE126074_P0_BrainCortex_SNAREseq_cDNA.barcodes.tsv")
# format gen info
geninfo <- read.table("GSE126074_P0_BrainCortex_SNAREseq_cDNA.genes.tsv")

indata = as.matrix(indata)
colnames(indata) = as.character(unlist(cellinfo))
rownames(indata) = as.character(unlist(geninfo))
 
data_rna = indata

#RNA_old
indata <- Matrix::readMM("GSE126074_AdBrainCortex_SNAREseq_cDNA.counts.mtx") 
# format cell info
cellinfo <- read.table("GSE126074_AdBrainCortex_SNAREseq_cDNA.barcodes.tsv")
# format gen info
geninfo <- read.table("GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv")

indata = as.matrix(indata)
colnames(indata) = as.character(unlist(cellinfo))
rownames(indata) = as.character(unlist(geninfo))
 
data_rna_old = indata


#ATAC
indata <- Matrix::readMM("GSE126074_P0_BrainCortex_SNAREseq_chromatin.counts.mtx") 
# format cell info
cellinfo <- read.table("GSE126074_P0_BrainCortex_SNAREseq_chromatin.barcodes.tsv")
# format gen info
geninfo <- read.table("GSE126074_P0_BrainCortex_SNAREseq_chromatin.peaks.tsv")

indata = as.matrix(indata)
colnames(indata) = as.character(unlist(cellinfo))
rownames(indata) = as.character(unlist(geninfo))
 
data_atac = indata

#ATAC_old
indata <- Matrix::readMM("GSE126074_AdBrainCortex_SNAREseq_chromatin.counts.mtx") 
# format cell info
cellinfo <- read.table("GSE126074_AdBrainCortex_SNAREseq_chromatin.barcodes.tsv")
# format gen info
geninfo <- read.table("GSE126074_AdBrainCortex_SNAREseq_chromatin.peaks.tsv")

indata = as.matrix(indata)
colnames(indata) = as.character(unlist(cellinfo))
rownames(indata) = as.character(unlist(geninfo))
 
data_atac_old = indata

######################################################
#Start downstream analysis
######################################################
#Working directory
setwd("/your_working_directory")
###########################################################
#ALL together
filename="Niche_Nature_Mm"
file_list = c("data_atac","data_atac_old","data_rna","data_rna_old")
file_list_mat = c("data_atac_mat","data_atac_mat_old", "data_rna_mat","data_rna_mat_old")
list_all = c(data_atac=SingleCellExperiment(assays = list(counts = data_atac)),
             data_atac_old=SingleCellExperiment(assays = list(counts = data_atac_old)),
             data_rna=SingleCellExperiment(assays = list(counts = data_rna)),
             data_rna_old=SingleCellExperiment(assays = list(counts = data_rna_old))
)

#ADD rowData the rownames for genes
list_all_Seurat = c()
rm(list = file_list)
rm(list = file_list_mat)
rm(file_list_mat)

###########################################################
#PRE-filtering
###########################################################
#Replace possible NAs by 0 in count matrix
for (n in 1:length(file_list)){
  counts(list_all[[n]])[is.na(counts(list_all[[n]]))] <- 0
}

#Remove Feature(genes,rna,etc) not expressed in any cell
for (n in 1:length(file_list)){
  keep_feature = rowSums(as.matrix(counts(list_all[[n]])) > 0) > 0
  list_all[[n]] = list_all[[n]][keep_feature, ]
}
rm(keep_feature)
#######################################################
#Define with features(genes,rna,etc) are the ERCC spike-ins and mitochondrial genes
for (n in 1:length(file_list)){
  isSpike(list_all[[n]], "ERCC") = grepl("ERCC-", rownames(list_all[[n]]))
  isSpike(list_all[[n]], "MT") = grepl("MT-", rownames(list_all[[n]]))
}
gc()

#######################################################
#Calculate de quality metrics
for (n in 1:length(file_list)){
  list_all[[n]]  = calculateQCMetrics(list_all[[n]],
                                      feature_controls = list(
                                        ERCC = isSpike(list_all[[n]], "ERCC"), 
                                        MT = isSpike(list_all[[n]], "MT")
                                      )
  )
  print (sprintf("QC added to sample -> %s", file_list[n]))
}
###########################################################
#PRE-PROCESSING -> FILTERING BAD QUALITY INSTANCES AND VARIABLES
###########################################################
### total_counts / total_features
lista_dataframes_counts = filter_instances_variables("total_counts", "Filtering_total_counts_plots_Histograms.pdf", "Filtering_total_counts_plots_Tables.pdf")
lista_dataframes_features = filter_instances_variables("total_features_by_counts", "Filtering_total_features_plots_Histograms.pdf", "Filtering_total_features_plots_Tables.pdf")

#######################################################
###ERCC/ MT 
lista_dataframes_ERCC = filter_ERCC_MT("is_spike_ERCC")
#lista_dataframes_MT = filter_ERCC_MT("is_spike_MT")

#######################################################
###Apply all previous created filters and filter umi based on those filters
for (n in 1:length(file_list)){
  #unlist(lista_dataframes_features[n]) &
  list_all[[n]]$use = (unlist(lista_dataframes_features[n]) & unlist(lista_dataframes_counts[n]))
  print (table(list_all[[n]]$use))
}
gc()

#######################################################
#Gene filtering after cell filtering !!!
#Low quality genes + MT + ERCC
#[x > 0] 1 before
lista_filter_genes = list()
for (n in 1:length(file_list)){
  #
  filter_genes = apply(
    counts(list_all[[n]][ , colData(list_all[[n]])$use]), 
    1, 
    function(x) length(x[x > 0]) >= 2
  )
  lista_filter_genes = append(lista_filter_genes, list(filter_genes))
  #Filter rows (genes)
  #& !unlist(lista_dataframes_MT[n]) 
  rowData(list_all[[n]])$use = (unlist(lista_filter_genes[n]) & !unlist(lista_dataframes_ERCC[n]))
  print(table(rowData(list_all[[n]])$use))
}

#lista_dataframes_MT
rm(list=setdiff(ls(), "list_all"))
gc()


#Second dataset
count_atac=counts(list_all[[1]][rowData(list_all[[1]])$use , colData(list_all[[1]])$use])
head(count_atac[1:10,1:10])
count_atac_old=counts(list_all[[2]][rowData(list_all[[2]])$use , colData(list_all[[2]])$use])
head(count_atac_old[1:10,1:10])

count_rna=counts(list_all[[3]][rowData(list_all[[3]])$use , colData(list_all[[3]])$use])
head(count_rna[1:10,1:10])
count_rna_old=counts(list_all[[4]][rowData(list_all[[4]])$use , colData(list_all[[4]])$use])
head(count_rna_old[1:10,1:10])

rm(list=setdiff(ls(), c("count_atac","count_atac_old","count_rna","count_rna_old","seurat_data_plot","seurat_data_plot_ATAC","seurat_data_plot_RNA","blankPlot")))
gc()

save.image("filtered_rna_atac_data.RData")
