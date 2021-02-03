# Machine Translation between paired Single Cell Multi Omics Data
This repository contains online data and LIBRA code to analyze and visualize paired Single-cell multiomics integration analysis and it's downstream analysis outputs. Metrics are also available for quantify outputs quality. 

For details, please visit https://www.biorxiv.org/content/10.1101/2021.01.27.428400v1.


## Material of interest


### Datasets:
For developing and testing the performance and quality of LIBRA, three sets of paired multimodality datasets were used.

- DataSet1, SNARE-seq1: GSE126074. SNARE-seq is a droplet-based single nucleus protocol which allows the simultaneous profiling of mRNA expression and chromatin accessibility. DataSet1 contains wild type Mouse brain cortex cells from neonatal (5.081 nuclei) and adult (10.309 nuclei) mouse brain for both paired scRNA and scATAC profiles.
- DataSet2, PBMC2: 10X Genomics website repository. 10X Genomics "Multiome ATAC + Gene Expression" protocol. The data set consists of 10.412 Human PBMC cells from a healthy donor without granulocytes removed by cell sorting.
- DataSet3, CITE-seq2,3: GSE128639. This data has been generated using the CITE-seq protocol, which provides paired profiles of (i) scRNAseq and an (ii) ADT panel for 25 antibodies. This set contains a total of 33.454 Human Bone Marrow cells for both data-modalities.


### Required libraries:
Please install the following R libraries before using LIBRA: devtools, keras, stringr, scclusteval, Seurat, ggplot2, Signac, scater, gridExtra, ggpubr, biomaRt, scran, cowplot, Matrix, data.table, GenomeInfoDb, EnsDb.Hsapiens.v75, patchwork, rhdf5


### LIBRA Workflow:
![workflow.png](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/workflow.png)


### Usage:
- Preprocessing
  - Removing low quality features and cells
    - Use "1_pre_analysis_Seurat3.R" under Code folder (or your own QC pipeline)
    > **Preprocessing: Feel free to add additional filtering to proposed ones as doublets removal specific libraries pipelines** 
    
- Analysis
  - Normalize and visualize with Seurat for independet samples (or other analysis pipeline that generates at least normalized matrix as output). If more than one sample is present for a given omic first integrate them as usual pipelines requires.
    - Use "2_analysis_Seurat3.R" under Code folder (or your own Normalization pipeline)
    > **Normalization method: Different normalization methods can enhance different LIBRA applications shuch as integration or prediction.** 

- LIBRA 
  - Use analysis_output.RData as input for LIBRA neural network
    - Use "LIBRA.R" under Code folder for networks training
    - Use "Metrics_LIBRA.R" for additional quality metrix computation
    > **Outputs: Different outputs generated during the training will be stored in the working directory.**
