###########################
#Prepare enviroment
set.seed(1234567)
options(stringsAsFactors = FALSE)

library("stringr")
library("scclusteval")
library("Seurat")

#SET WORKING DIRECTORY WHERE ALL LIBRA fine-tune ARE STORED
setwd("set_your_path_here")

#IMPORT ALL OUTCOMES GENERATED
temp = list.files(pattern=paste0("test_n_"))
myfiles = lapply(temp, FUN = read.csv, row.names = 1)

names = c()
scores_dt1 = c() 
scores_dt2 = c() 
n_clusters = c()                           
for (i in (1:length(myfiles))){
  names = append(names, temp[i]) 
  xx= str_sub(temp[i], -6, -5)
  
  LIBRA =  myfiles[[i]]
  integrated_idents@reductions$libra = integrated_idents@reductions$pca
  rownames(LIBRA) = rownames(integrated_idents@reductions$libra@cell.embeddings)
  integrated_idents@reductions$libra@cell.embeddings = as.matrix(LIBRA)
  integrated_idents@reductions$pca@cell.embeddings = integrated_idents@reductions$libra@cell.embeddings
  
  DefaultAssay(integrated_idents) <- "RNA"
  
  integrated_idents = FindNeighbors(object = integrated_idents, reduction = "pca", dims = 1:xx)
  integrated_idents = FindClusters(object = integrated_idents, resolution = 0.8)
  
  jacard_vector_rna_score = c()
  jacard=PairWiseJaccardSetsHeatmap(RNA_idents@active.ident, integrated_idents@active.ident,
                         show_row_dend = F, show_column_dend = F,
                         cluster_row = F, cluster_column =F)
  jacard_distance = jacard@matrix
  is.na(jacard_distance) <- jacard_distance==0
  jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]
  jacard_vector_rna_score = append(jacard_vector_rna_score, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))
  print(jacard_vector_rna_score)

  scores_dt1 = append(scores_dt1, jacard_vector_rna_score) 

  jacard_vector_rna_score = c()
  jacard=PairWiseJaccardSetsHeatmap(ADT_idents@active.ident, integrated_idents@active.ident,
                         show_row_dend = F, show_column_dend = F,
                         cluster_row = F, cluster_column =F)
  jacard_distance = jacard@matrix
  is.na(jacard_distance) <- jacard_distance==0
  jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]
  jacard_vector_atac_score = append(jacard_vector_rna_score, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))
  print(jacard_vector_atac_score)
  
  scores_dt2 = append(scores_dt2, jacard_vector_atac_score) 
  
  clusters = (max(as.numeric(levels(integrated_idents@meta.data$seurat_clusters)))+1)
  n_clusters = append(n_clusters, clusters)
} 
names
scores_dt1
scores_dt2
n_clusters