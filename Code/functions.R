##########################################################################
#1. Library size in terms of umis
#2. by counting unique genes detected in each sample
filter_instances_variables<- function(total_values,output_pdf_name,output_pdf_name_2) {
  lista_dataframes = list()
  #
  pdf(file=output_pdf_name)
  for (n in 1:length(file_list)){
    min<-0
    max<-max(list_all[[n]][[total_values]])
    gup<-quantile(list_all[[n]][[total_values]],probs=(.9))
    gdown<-quantile(list_all[[n]][[total_values]],probs=(.1))
    hist(
      main=file_list[n],
      list_all[[n]][[total_values]],
      breaks = 100,
      xlim=c(min,max),
      xlab="Total counts"
    )
    abline(v = c(gdown,gup), col = "red")
    #Filter library size
    filter_by <- (list_all[[n]][[total_values]] > gdown & list_all[[n]][[total_values]] < gup)
    assign(paste(file_list[n],"filter_by",sep = "_"), filter_by)
    
    #Intialize
    lista_dataframes = append(lista_dataframes, list(get(paste(file_list[n],"filter_by",sep = "_"))))
    
    print (sprintf("Filtering %s plots_Histograms added to pdf -> %s", file_list[n],output_pdf_name))
  }
  dev.off()
  
  #Export tables with results
  pdf(file=output_pdf_name_2)
  for (n in 1:length(file_list)){
    #Plot table with amount of elements filtered
    df=knitr::kable(
      as.data.frame(table(get(paste(file_list[n],"filter_by",sep = "_")))),
      booktabs = TRUE,
      row.names = FALSE,
      caption = 'Filtered out (FALSE)',
      col.names = c(paste(file_list[n],"filter_by",sep = "_"), "Freq")
    )
    #Export tables
    p<-tableGrob(df)
    grid.arrange(p)
    
    print (sprintf("Saving %s table of filtered out elements to pdf -> %s", file_list[n],output_pdf_name))
  }
  dev.off()
  
  return(lista_dataframes)
}

##########################################################################
#Fuction for filter ERCC and MT
filter_ERCC_MT<- function(is_spike) {
  lista_filtrada = list()
  #
  for (n in 1:length(file_list)){
    #ERCC
    filtrada_by <- list_all[[n]]@int_elementMetadata@listData[[is_spike]]
    assign(paste(file_list[n],"filtrada_by",sep = "_"), filtrada_by)
    lista_filtrada = append(lista_filtrada, list(get(paste(file_list[n],"filtrada_by",sep = "_"))))
  }
  
  return(lista_filtrada)
}

