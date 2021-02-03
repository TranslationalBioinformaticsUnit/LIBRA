#Summarizations of metrics computed on LIBRA 

################################
setwd("/your_output_path_from_LIBRA")
#Summarize JACARD distances into 1 value (PPJI)
latent_dim_val = c(5,10,20,40,60)
#latent_dim_val = c(10) for explicit case
for (m in latent_dim_val){
  jacard_vector_rna = c()
  jacard_vector_rna_test = c()
  jacard_vector_atac = c()
  for (n in 1:10){
    load(paste0("RNA_and_ATAC_",m,"_nodes_",n,"interaction.RData")) 
    
    #RNA
    jacard=PairWiseJaccardSetsHeatmap(pbmc_rna@active.ident, pbmc.atac_libra@active.ident, #pbmc.atac_libra and pbmc.rna_libra are the same
                                      show_row_dend = F, show_column_dend = F,
                                      cluster_row = F, cluster_column =F)
    
    #PPJI summarization
    jacard_distance = jacard@matrix
    is.na(jacard_distance) <- jacard_distance==0
    jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]
    jacard_vector_rna = append(jacard_vector_rna, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))
    print(jacard_vector_rna)
    
    #Validation split data
    validation_split_cells = tail(names(pbmc.rna@active.ident),floor((length(names(pbmc.rna@active.ident))*20)/100)) #If validation_split = 0.2
    
    #Compare to original clustering using JACARD ONLY ON TEST DATA
    A = pbmc.rna@active.ident[names(pbmc.rna@active.ident) %in% validation_split_cells]
    B = pbmc.rna_libra@active.ident[names(pbmc.rna_libra@active.ident) %in% validation_split_cells]
    
    jacard=PairWiseJaccardSetsHeatmap(A, B,
                                      show_row_dend = F, show_column_dend = F,
                                      cluster_row = F, cluster_column =F)
    #Mean for non 0 vales
    jacard_distance = jacard@matrix
    is.na(jacard_distance) <- jacard_distance==0
    jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]
    jacard_vector_rna_test = append(jacard_vector_rna_test, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))
    
    #ATAC
    jacard=PairWiseJaccardSetsHeatmap(pbmc_atac@active.ident, pbmc.atac_libra@active.ident, #pbmc.atac_libra and pbmc.rna_libra are the same
                                      show_row_dend = F, show_column_dend = F,
                                      cluster_row = F, cluster_column =F)
    
    #Meand for non 0 vales
    jacard_distance = jacard@matrix
    is.na(jacard_distance) <- jacard_distance==0
    jacard_distance = jacard_distance[,colSums(is.na(jacard_distance))<nrow(jacard_distance)]
    jacard_vector_atac = append(jacard_vector_atac, mean(apply(jacard_distance,1,sum, na.rm=TRUE)))    
    
  }
  write.table(jacard_vector_rna, file=paste0("RNA_nodes_", m,"Jacard_RNA_summary.txt"), sep="\t")
  write.table(jacard_vector_rna_test, file=paste0("RNA_test_nodes_", m,"Jacard_RNA_summary.txt"), sep="\t")
  write.table(jacard_vector_atac, file=paste0("ATAC_nodes_", m,"Jacard_ATAC_summary.txt"), sep="\t")
}  


#Highlight bimodal cells where are there? (test cells)
for (m in latent_dim_val){
  percentage_bimodal_test_cells = c()
  for (n in 1:10){
    load(paste0("RNA_and_ATAC_",m,"_nodes_",n,"interaction.RData")) 
    
    #Get distances RNA to ATAC in shared space / OBTAIN THE LOCAL MINIMUN
    minimun_local = distance_vector[which.min(abs(diff(distance_vector)))]
    
    #Get the local minimun to the min distances
    distance_vector[distance_vector < minimun_local] <- 0
    
    #Get the local minimun to the max distances
    distance_vector[distance_vector >= minimun_local] <- 1 
    
    pbmc.atac@meta.data$distance_vector = distance_vector
    bimodal_state = rownames(pbmc.atac@meta.data)[pbmc.atac@meta.data$distance_vector == 1]
    
    percentage_bimodal_test_cells = append(percentage_bimodal_test_cells,(sum(bimodal_state %in% validation_split_cells) * 100) / length(validation_split_cells))
  }
  percentage_bimodal_test_cells_mean = mean(percentage_bimodal_test_cells)
  percentage_bimodal_test_cells_median = median(percentage_bimodal_test_cells)
  percentage_bimodal_test_cells_sd = sd(percentage_bimodal_test_cells)
  
  percentage_bimodal_test_cells = append(percentage_bimodal_test_cells, percentage_bimodal_test_cells_mean)
  percentage_bimodal_test_cells = append(percentage_bimodal_test_cells, percentage_bimodal_test_cells_median)
  percentage_bimodal_test_cells = append(percentage_bimodal_test_cells, percentage_bimodal_test_cells_sd)
  
  write.table(percentage_bimodal_test_cells, file=paste0("BIMODAL_CELLS", m,"_explanation_test_cells.txt"), sep="\t")
  percentage_bimodal_test_cells_mean
}  

#Summarize distances of bimodal_1 and bimodal_2
latent_dim_val = c(5,10,20,40,60)
#latent_dim_val = c(10) for explicit case
for (m in latent_dim_val){
  bimodal_1_all = c()
  bimodal_2_all = c()
  bimodal_1_mean_all = c()
  bimodal_1_mean_all = c()
  bimodal_1_mean_all = c()
  bimodal_2_mean_all = c()
  bimodal_2_mean_all = c()
  bimodal_2_mean_all = c()
  bimodal_1_mean = c()
  bimodal_1_median = c()
  bimodal_1_sd = c()
  bimodal_2_mean = c()
  bimodal_2_median = c()
  bimodal_2_sd = c()
  for (n in 1:10){
    load(paste0("RNA_and_ATAC_",m,"_nodes_",n,"interaction.RData")) 
    
    #Get distances RNA to ATAC in shared space / OBTAIN THE LOCAL MINIMUN
    d <- density(distance_vector)
    minimun_local = optimize(approxfun(d$x,d$y),interval=c(((min(distance_vector)*125)/100),((max(distance_vector)*60)/100)))$minimum
    
    #Min distance
    min_distance = min(distance_vector)
    modal_1 = distance_vector[distance_vector < minimun_local & distance_vector > min_distance]
    
    #Max distance
    max_distance = max(distance_vector)
    modal_2 = distance_vector[distance_vector > minimun_local & distance_vector < max_distance]
    
    mean_range = (abs(min_distance - minimun_local) / 4)
    
    #modal_1
    h = hist(modal_1, breaks=100000)
    i = which.max(h$counts)
    max_peak_modal_1 = h$mids[i]
    
    peak_upper = max_peak_modal_1 + mean_range
    peak_botton = max_peak_modal_1 - mean_range
    
    if (peak_upper < max(modal_1)) {
      peak_upper = peak_upper
    } else {
      peak_upper = max(modal_1)
    }
    if (peak_botton > min(modal_1)) {
      peak_botton = peak_botton
    } else {
      peak_botton = min(modal_1)
    }   
    
    modal_1_peak_zone = modal_1[(modal_1 > peak_botton) & (modal_1 < peak_upper)]
    
    bimodal_1_mean = append(bimodal_1_mean, mean(modal_1_peak_zone))
    bimodal_1_median = append(bimodal_1_median, median(modal_1_peak_zone))
    bimodal_1_sd = append(bimodal_1_sd, sd(modal_1_peak_zone))
    
    #modal_2
    h = hist(modal_2, breaks=100000)
    i = which.max(h$counts)
    max_peak_modal_2 = h$mids[i]
    
    peak_upper = max_peak_modal_2 + mean_range
    peak_botton = max_peak_modal_2 - mean_range
    
    if (peak_upper < max(modal_2)) {
      peak_upper = peak_upper
    } else {
      peak_upper = max(modal_2)
    }
    if (peak_botton > min(modal_2)) {
      peak_botton = peak_botton
    } else {
      peak_botton = min(modal_2)
    }   
    
    modal_2_peak_zone = modal_2[(modal_2 > peak_botton) & (modal_2 < peak_upper)]
    
    bimodal_2_mean = append(bimodal_2_mean, mean(modal_2_peak_zone))
    bimodal_2_median = append(bimodal_2_median, median(modal_2_peak_zone))
    bimodal_2_sd = append(bimodal_2_sd, sd(modal_2_peak_zone)) 
  }
  bimodal_1_mean_all = mean(bimodal_1_mean)
  bimodal_1_median_all = mean(bimodal_1_median)
  bimodal_1_sd_all = mean(bimodal_1_sd)
  
  bimodal_2_mean_all = mean(bimodal_2_mean)
  bimodal_2_median_all = mean(bimodal_2_median)
  bimodal_2_sd_all = mean(bimodal_2_sd)
  
  bimodal_1_all = append(bimodal_1_all, bimodal_1_mean_all)
  bimodal_1_all = append(bimodal_1_all, bimodal_1_median_all)
  bimodal_1_all = append(bimodal_1_all, bimodal_1_sd_all)
  print(bimodal_1_all)
  
  bimodal_2_all = append(bimodal_2_all, bimodal_2_mean_all)
  bimodal_2_all = append(bimodal_2_all, bimodal_2_median_all)
  bimodal_2_all = append(bimodal_2_all, bimodal_2_sd_all)
  print(bimodal_2_all)
  
  write.table(bimodal_1_all, file=paste0("BIMODAL_CELLS", m,"modal_1_peak_zone.txt"), sep="\t")
  write.table(bimodal_2_all, file=paste0("BIMODAL_CELLS", m,"modal_2_peak_zone.txt"), sep="\t")
}  
