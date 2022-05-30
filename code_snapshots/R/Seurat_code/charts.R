##########################################################################
#Blank chart
blankPlot <- ggplot()+geom_blank(aes(1,1)) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

##########################################################################
#Visualize data for filtering propouse
seurat_data_plot<- function(output_pdf_name) {
  pdf(file=output_pdf_name, width=20, height=6)
  #
  for (n in 1:length(file_list)){
    vln=VlnPlot(
      object = list_all_Seurat[[n]], 
      #, "percent.mt", "percent.rb"
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
    )
    
    scatterPlot = ggplot(list_all_Seurat[[n]]@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point() + geom_rug()
    
    xdensity = ggplot(list_all_Seurat[[n]]@meta.data, aes(nCount_RNA)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none") +
      geom_vline(xintercept = mean(list_all_Seurat[[n]]@meta.data$nCount_RNA), color='red') +
      geom_vline(xintercept = median(list_all_Seurat[[n]]@meta.data$nCount_RNA), color='blue')
    
    ydensity = ggplot(list_all_Seurat[[n]]@meta.data, aes(nFeature_RNA)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none")+
      geom_vline(xintercept = mean(list_all_Seurat[[n]]@meta.data$nFeature_RNA), color='red') +
      geom_vline(xintercept = median(list_all_Seurat[[n]]@meta.data$nFeature_RNA), color='blue')
    
    scater_gen_umi = gridExtra::arrangeGrob(xdensity, blankPlot, scatterPlot, ydensity, 
                                  ncol=2, nrow=2, widths=c(4, 2), heights=c(1.4, 4))
    
    plot(ggarrange(vln, scater_gen_umi,
              labels = c("A", "B"),
              ncol = 2, nrow = 1,
              widths = c(1, 1.5))
         )
  }
  dev.off()
}

##########################################################################
#Visualize data for filtering propouse FOR ATAC 1 SAMPLE
seurat_data_plot_ATAC<- function(output_pdf_name) {
  pdf(file=output_pdf_name, width=20, height=6)
  #
    vln=VlnPlot(
      object = pbmc.atac, 
      #, "percent.mt", "percent.rb"
      features = c("nFeature_ATAC", "nCount_ATAC", "percent.mt", "percent.rb")
    )
    
    scatterPlot = ggplot(pbmc.atac@meta.data, aes(x=nCount_ATAC, y=nFeature_ATAC)) +
      geom_point() + geom_rug()
    
    xdensity = ggplot(pbmc.atac@meta.data, aes(nCount_ATAC)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none") +
      geom_vline(xintercept = mean(pbmc.atac@meta.data$nCount_ATAC), color='red') +
      geom_vline(xintercept = median(pbmc.atac@meta.data$nCount_ATAC), color='blue')
    
    ydensity = ggplot(pbmc.atac@meta.data, aes(nFeature_ATAC)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none")+
      geom_vline(xintercept = mean(pbmc.atac@meta.data$nFeature_ATAC), color='red') +
      geom_vline(xintercept = median(pbmc.atac@meta.data$nFeature_ATAC), color='blue')
    
    scater_gen_umi = gridExtra::arrangeGrob(xdensity, blankPlot, scatterPlot, ydensity, 
                                  ncol=2, nrow=2, widths=c(4, 2), heights=c(1.4, 4))
    
    plot(ggarrange(vln, scater_gen_umi,
              labels = c("A", "B"),
              ncol = 2, nrow = 1,
              widths = c(1, 1.5))
         )
  dev.off()
}
##########################################################################
#Visualize data for filtering propouse FOR RNA 1 SAMPLE
seurat_data_plot_RNA<- function(output_pdf_name) {
  pdf(file=output_pdf_name, width=20, height=6)
  #
    vln=VlnPlot(
      object = pbmc.rna, 
      #, "percent.mt", "percent.rb"
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
    )
    
    scatterPlot = ggplot(pbmc.rna@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) +
      geom_point() + geom_rug()
    
    xdensity = ggplot(pbmc.rna@meta.data, aes(nCount_RNA)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none") +
      geom_vline(xintercept = mean(pbmc.rna@meta.data$nCount_RNA), color='red') +
      geom_vline(xintercept = median(pbmc.rna@meta.data$nCount_RNA), color='blue')
    
    ydensity = ggplot(pbmc.rna@meta.data, aes(nFeature_RNA)) + 
      geom_density(alpha=.5) + 
      theme(legend.position = "none")+
      geom_vline(xintercept = mean(pbmc.rna@meta.data$nFeature_RNA), color='red') +
      geom_vline(xintercept = median(pbmc.rna@meta.data$nFeature_RNA), color='blue')
    
    scater_gen_umi = gridExtra::arrangeGrob(xdensity, blankPlot, scatterPlot, ydensity, 
                                  ncol=2, nrow=2, widths=c(4, 2), heights=c(1.4, 4))
    
    plot(ggarrange(vln, scater_gen_umi,
              labels = c("A", "B"),
              ncol = 2, nrow = 1,
              widths = c(1, 1.5))
         )
  dev.off()
}
















