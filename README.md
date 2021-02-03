# Machine Translation between paired Single Cell Multi Omics Data
This repository contains online data and LIBRA code to analyze and visualize paired Single-cell multiomics integration analysis and it's downstream analysis outputs. Metrics are also available for quantify outputs quality. 

For details, please visit https://www.biorxiv.org/content/10.1101/2021.01.27.428400v1.


## Material of interest


### Datasets:
For developing and testing the performance and quality of LIBRA, three sets of paired multimodality datasets were used. The original datasets are stored under "Data" folder.

- DataSet1, SNARE-seq1: GSE126074. SNARE-seq is a droplet-based single nucleus protocol which allows the simultaneous profiling of mRNA expression and chromatin accessibility. DataSet1 contains wild type Mouse brain cortex cells from neonatal (5.081 nuclei) and adult (10.309 nuclei) mouse brain for both paired scRNA and scATAC profiles.
- DataSet2, PBMC2: 10X Genomics website repository. 10X Genomics "Multiome ATAC + Gene Expression" protocol. The data set consists of 10.412 Human PBMC cells from a healthy donor without granulocytes removed by cell sorting.
- DataSet3, CITE-seq2,3: GSE128639. This data has been generated using the CITE-seq protocol, which provides paired profiles of (i) scRNAseq and an (ii) ADT panel for 25 antibodies. This set contains a total of 33.454 Human Bone Marrow cells for both data-modalities.


### Required libraries:
  - library(devtools)
  - library("keras")
  - library("stringr")
  - library("scclusteval")
  - library("keras")
  - library("Seurat")
  - library("ggplot2")
  - library("Signac") 


### LIBRA Workflow
![Test Image 1](F1_2.pdf)
  

  https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/F1_2.pdf
  
