###########################
#Prepare enviroment
set.seed(1234567)
options(stringsAsFactors = FALSE)

#Uncomment for test // Comment for deploy
#tensorflow::tf$random$set_seed(1234567)
#use_session_with_seed(1234567, disable_gpu = TRUE, disable_parallel_cpu = TRUE, quiet = FALSE)

#Required libraries
library("keras")
library("stringr")
library("scclusteval")
library("keras")
library("Seurat")
library("ggplot2")
library("Signac")

K <- keras::backend()
tensorflow::tf$compat$v1$disable_eager_execution()


#--------------------------------------------------------------
#--------------------------------------------------------------
#INPUT DATA
#Seurat slots required may change base on omics used, in this example we make use of scRNA and scATAC for "pbmc" dataset
setwd("/your_path")
load("your_standard_analysis_Seurat_3_output.RData")

DefaultAssay(pbmc) <- 'RNA'
rna = t(data.frame(pbmc@assays$RNA@data))
atac = t(data.frame(pbmc@assays$ATAC@data))

#Set of overlapping cells for evade possible unpaired cells
overlap_cells = rownames(rna)[rownames(rna) %in% rownames(atac)]

#Subset data base on paired cells only
rna = rna[overlap_cells,]
atac = atac[overlap_cells,]

#Train data x_train(input of LIBRA) x_train2(output of LIBRA)
x_train=(atac)
x_train2=(rna)


#--------------------------------------------------------------
#--------------------------------------------------------------
#Store all results from 10 iterations for desired middle layer sizes
latent_dim_val = c(5,10,20,40,60)

#--------------------------------------------------------------
#--------------------------------------------------------------
#Store all results from 10 iterations 
setwd("/your_output_path_from_LIBRA")

for (m in latent_dim_val){
  all_mean_square_errors_rna = c()
  loss_mean_square_errors_rna = c()
  all_mean_square_errors_atac = c()
  loss_mean_square_errors_atac = c()
  
  all_mean_square_errors_all_rna = c()
  loss_mean_square_errors_all_rna = c()
  all_mean_square_errors_all_atac = c()
  loss_mean_square_errors_all_atac = c()
  
  for (n in 1:10){ #Selected 10 as the number of times to run LIBRA for each middle layer size and configuration
    # Parameters --------------------------------------------------------------
    batch_size <- 7000L #Depends on dataset number of cells
    original_dim <- ncol(x_train) 
    latent_dim <- m 
    intermediate_dim_1 <- 512L 
    intermediate_dim_2 <- 256L 
    epochs <- 1500L 
    dropout <- 0.2
    
    # Model definition --------------------------------------------------------
    model <- keras_model_sequential() 
    model %>% 
      layer_dense(units = intermediate_dim_1, input_shape = ncol(x_train)) %>% 
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout) %>%  
      layer_dense(units = intermediate_dim_2, input_shape = ncol(x_train)) %>% 
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout) %>%  
      layer_dense(units = m, name = "bottleneck") %>%
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout) %>% 
      layer_dense(units = intermediate_dim_2) %>%
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout) %>% 
      layer_dense(units = intermediate_dim_1) %>%
      layer_activation_leaky_relu() %>%
      layer_dense(units = ncol(x_train2), name = "output")
    
    
    model %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(),
      metrics = list("mean_squared_error")
    )
    
    #Store the data for each interations in the training process
    checkpoint <- callback_model_checkpoint(
      monitor = "loss",
      filepath = paste0(m,"_nodes",n,"_interation","_ae_1_model.hdf5"), 
      save_best_only = TRUE, 
      period = 1,
      verbose = 1
    )
    
    early_stopping <- callback_early_stopping(monitor = "loss", patience = 20)
    
    history_atac <- model %>% fit(
      x = x_train, 
      y = x_train2, 
      epochs = epochs, 
      batch_size = batch_size,
      callbacks = list(checkpoint, early_stopping),
      validation_split = 0 #Change to 0.2 for optimize (0.2 used for first try quality measurement)
    )
    
    
    #####################################################################
    # Outcomes ----------------------------------------------------------
    #Output layer output
    layer_name = 'output'
    output_layer_model <- keras_model(inputs = model$input,
                                      outputs = get_layer(model, layer_name)$output)
    output_output <- predict(output_layer_model, x_train)
    
    #Intermediate layer output
    layer_name = 'bottleneck'
    intermediate_layer_model <- keras_model(inputs = model$input,
                                            outputs = get_layer(model, layer_name)$output)
    intermediate_output <- predict(intermediate_layer_model, x_train)
    
    #Set correct names if not done before
    rownames(intermediate_output) = rownames(x_train)
    rownames(intermediate_output) = gsub(".", '-', rownames(intermediate_output), fixed = T)
    
    #FOR ATAC (input_omic)      
    pbmc.atac = pbmc
    
    #Add computed shared components on new object "_libra" without overwriting original one
    pbmc.atac_libra = pbmc.atac
    pbmc.atac_libra@reductions$lsi@cell.embeddings = intermediate_output[match(rownames(pbmc.atac_libra@reductions$lsi@cell.embeddings),rownames(intermediate_output)),]
    
    pbmc.atac_libra@reductions$lsi@assay.used = "ATAC"
    pbmc.atac_libra = FindNeighbors(object = pbmc.atac_libra, reduction = "lsi", dims = 1:m)
    pbmc.atac_libra= FindClusters(object = pbmc.atac_libra, resolution = 0.8, n.start = 1000)
    
    #Create plot base on original components
    pdf(file=paste0("ATAC_on_ATACproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.atac_libra, reduction = "umap.atac", group.by = "celltype")
    p2 = DimPlot(pbmc.atac_libra, reduction = "umap.atac", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    pdf(file=paste0("ATAC_on_RNAproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.atac_libra, reduction = "umap.rna", group.by = "celltype")
    p2 = DimPlot(pbmc.atac_libra, reduction = "umap.rna", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    pdf(file=paste0("ATAC_on_INTEGRATEDproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.atac_libra, reduction = "wnn.umap", group.by = "celltype")
    p2 = DimPlot(pbmc.atac_libra, reduction = "wnn.umap", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    #Obtain UMAP base on new shared components
    pbmc.atac_libra = RunUMAP(pbmc.atac_libra, reduction = "lsi", dims = 1:m)
    
    #Create plot base on shared components
    pdf(file=paste0("ATAC_on_new_INTEGRATEDproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.atac_libra, reduction = "umap", group.by = "celltype")
    p2 = DimPlot(pbmc.atac_libra, reduction = "umap", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off() 
    
    #Compare to original clustering using JACARD
    jacard=PairWiseJaccardSetsHeatmap(pbmc.atac@active.ident, pbmc.atac_libra@active.ident,
                                      show_row_dend = F, show_column_dend = F,
                                      cluster_row = F, cluster_column =F)
    pdf(file=paste0("ATAC_clusters_JACARD_",m,"_nodes",n,"_interation_ALL_DATA",".pdf"), width=12, height=6)
    print(jacard)
    dev.off()
    
    
    #####################################################################
    #####################################################################
    
    #Map the RNA to shared space generated by LIBRA
    ae <- keras_model_sequential() 
    ae %>% 
      layer_dense(units = intermediate_dim_1, input_shape = ncol(x_train2)) %>% 
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout) %>%  
      layer_dense(units = intermediate_dim_2, input_shape = ncol(x_train2)) %>% 
      layer_activation_leaky_relu() %>%
      layer_dropout(rate = dropout) %>% 
      layer_dense(units = m, name = "bottleneck") %>%
      layer_activation_leaky_relu() %>%
      layer_dense(units = ncol(intermediate_output), name = "output")
    
    #Compile the model with appropriate loss function, optimizer, and metrics
    ae %>% compile(
      loss = "mean_squared_error",
      optimizer = optimizer_adam(),
      metrics = list("mean_squared_error")
    )
    
    #Store the data for each interations in the training process
    checkpoint <- callback_model_checkpoint(
      monitor = "loss",
      filepath = paste0(m,"_nodes",n,"_interation","_ae_2_model.hdf5"), 
      save_best_only = TRUE, 
      period = 1,
      verbose = 1
    )
    
    early_stopping <- callback_early_stopping(monitor = "loss", patience = 20)
    
    history_rna <- ae %>% fit(
      x = x_train2, 
      y = intermediate_output, 
      epochs = 1500, 
      batch_size = 7000, #Depends on dataset number of cells
      callbacks = list(checkpoint, early_stopping),
      validation_split = 0 #Change to 0.2 for optimize (0.2 used for first try quality measurement)
    )
    
    
    #####################################################################
    # Outcomes ----------------------------------------------------------
    #Output layer output
    layer_name = 'output'
    output_layer_model <- keras_model(inputs = ae$input,
                                      outputs = get_layer(ae, layer_name)$output)
    output_output <- predict(output_layer_model, x_train2)
    
    #Set correct names
    rownames(output_output) = rownames(x_train) 
    
    #FOR RNA (output_omic)
    pbmc.rna = pbmc
    
    #Add computed shared components
    pbmc.rna_libra = pbmc.rna
    pbmc.rna_libra@reductions$pca@cell.embeddings = intermediate_output[match(rownames(pbmc.rna_libra@reductions$pca@cell.embeddings),rownames(intermediate_output)),]
    
    pbmc.rna_libra@reductions$pca@assay.used = "RNA"
    pbmc.rna_libra = FindNeighbors(object = pbmc.rna_libra, reduction = "pca", dims = 1:m)
    pbmc.rna_libra = FindClusters(object = pbmc.rna_libra, resolution = 0.8, n.start = 1000)
    
    #Create plot base on original components
    pdf(file=paste0("RNA_on_ATACproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.rna_libra, reduction = "umap.atac", group.by = "celltype")
    p2 = DimPlot(pbmc.rna_libra, reduction = "umap.atac", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    pdf(file=paste0("RNA_on_RNAproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.rna_libra, reduction = "umap.rna", group.by = "celltype")
    p2 = DimPlot(pbmc.rna_libra, reduction = "umap.rna", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    pdf(file=paste0("RNA_on_INTEGRATEDproj_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=18, height=6)
    p1 = DimPlot(pbmc.rna_libra, reduction = "wnn.umap", group.by = "celltype")
    p2 = DimPlot(pbmc.rna_libra, reduction = "wnn.umap", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    #Same clusters as first one
    pbmc.rna_libra@graphs$RNA_snn = pbmc.atac_libra@graphs$RNA_snn
    #pbmc.rna_libra@graphs$RNA_snn@assay.used = "RNA"
    pbmc.rna_libra = FindClusters(object = pbmc.rna_libra, resolution = 2, n.start = 1000)
    pdf(file=paste0("ssssRNA_clusters_original_UMAP_",m,"_nodes",n,"_interation",".pdf"), width=12, height=6)
    p1 = DimPlot(pbmc.rna_libra, reduction = "wnn.umap", group.by = "celltype")
    p2 = DimPlot(pbmc.rna_libra, reduction = "wnn.umap", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    #Obtain UMAP base on new shared components
    pbmc.rna_libra = RunUMAP(pbmc.rna_libra, reduction = "pca", dims = 1:m)
    
    ##Create plot base on shared components
    pdf(file=paste0("RNA_clusters_shared_UMAP_2network",m,"_nodes",n,"_interation",".pdf"), width=12, height=6)
    p1 = DimPlot(pbmc.rna_libra, reduction = "umap", group.by = "celltype")
    p2 = DimPlot(pbmc.rna_libra, reduction = "umap", label = TRUE)
    print(CombinePlots(list(p1, p2)))
    dev.off()
    
    #Compare to original clustering using JACARD
    jacard=PairWiseJaccardSetsHeatmap(pbmc.rna@active.ident, pbmc.rna_libra@active.ident,
                                      show_row_dend = F, show_column_dend = F,
                                      cluster_row = F, cluster_column =F)
    pdf(file=paste0("RNA_clusters_JACARD_2network",m,"_nodes",n,"_interation_ALL_DATA",".pdf"), width=12, height=6)
    print(jacard)
    dev.off()
    
    
    #####################################################################
    # Statistic for models ----------------------------------------------
    #Store all results from 10 iterations 
    
    #FOR ATAC (Input omic)
    all_mean_square_errors_atac = append(all_mean_square_errors_atac, min(history_atac$metrics$mean_squared_error))
    
    loss_atac = evaluate(model, x = x_train, y = x_train2)
    all_mean_square_errors_atac = append(all_mean_square_errors_atac, loss_atac$mean_squared_error)
    
    #FOR RNA (Output omic)
    all_mean_square_errors_rna = append(all_mean_square_errors_rna, min(history_rna$metrics$mean_squared_error))
    
    loss_rna <- evaluate(ae, x = x_train2, y = intermediate_output)
    loss_mean_square_errors_rna = append(loss_mean_square_errors_rna, loss_rna$mean_squared_error)
    
    #####################################################################
    # Distance between cells RNA and ATAC for each iteration 
    #Euclidean distance measure distribution for cell distances
    library(rdist)
    euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
    distance_vector = c()
    for (o in 1:nrow(rna)){
      distance_vector = append(distance_vector, euc.dist(intermediate_output[o,], output_output[o,]))
      print(n)
    }
    
    pdf(file=paste0("Distances_eu_RNA_vs_ATAC_",m,"_nodes",n,"_interation.pdf"), width=12, height=6)
    myplot = ggplot(data.frame(distance_vector), aes(x=distance_vector)) + geom_density() +
      geom_vline(aes(xintercept=median(distance_vector), color="median"), linetype="dashed", size=1) +
      geom_vline(aes(xintercept=mean(distance_vector), color="mean"), linetype="dashed", size=1) +
      scale_color_manual(name = "statistics", values = c(median = "blue", mean = "red"))
    print(myplot)
    dev.off()
    
    #Export Rdata
    save.image(paste0("RNA_and_ATAC_",m,"_nodes_",n,"interaction.RData"))
  }
  #FOR RNA
  all_mean_square_errors_all_rna = append(all_mean_square_errors_all_rna, list(all_mean_square_errors_rna))
  assign(paste0("RNA_nodes_", m,"_training_all_mean_square_errors"), all_mean_square_errors_all_rna)
  write.table(get(paste0("RNA_nodes_", m,"_training_all_mean_square_errors")), file=paste0("RNA_nodes_", m,"_training_all_mean_square_errors.txt"), sep="\t")
  
  loss_mean_square_errors_all_rna = append(loss_mean_square_errors_all_rna, list(loss_mean_square_errors_rna))
  assign(paste0("RNA_nodes_", m,"_test_all_mean_square_errors"), loss_mean_square_errors_all_rna)
  write.table(get(paste0("RNA_nodes_", m,"_test_all_mean_square_errors")), file=paste0("RNA_nodes_", m,"_test_all_mean_square_errors.txt"), sep="\t")
  
  #FOR ATAC
  all_mean_square_errors_all_atac = append(all_mean_square_errors_all_atac, list(all_mean_square_errors_atac))
  assign(paste0("ATAC_nodes_", m,"_training_all_mean_square_errors"), all_mean_square_errors_all_atac)
  write.table(get(paste0("ATAC_nodes_", m,"_training_all_mean_square_errors")), file=paste0("ATAC_nodes_", m,"_training_all_mean_square_errors.txt"), sep="\t")
  
  loss_mean_square_errors_all_atac = append(loss_mean_square_errors_all_atac, list(loss_mean_square_errors_atac))
  assign(paste0("ATAC_nodes_", m,"_test_all_mean_square_errors"), loss_mean_square_errors_all_atac)
  write.table(get(paste0("ATAC_nodes_", m,"_test_all_mean_square_errors")), file=paste0("ATAC_nodes_", m,"_test_all_mean_square_errors.txt"), sep="\t")
}
