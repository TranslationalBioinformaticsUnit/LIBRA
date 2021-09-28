# LIBRA
*Machine Translation between paired Single-Cell Multi-Omics Data*

This repository contains:
- [Online data](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/data/) employed on LIBRA manuscript.
- [Seurat code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/) employed to analyze and visualize integration analysis performed by Seurat.
- [LIBRA code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/) to analyze and visualize paired Single-cell multi-omics integration and prediction analysis.
- [LIBRA metrics](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/) are also available for quantifying outputs quality base on PPJI preservation measurement.

For more information , check out the [online manuscript](https://www.biorxiv.org/content/10.1101/2021.01.27.428400v1) currently at biorxiv repository (will be updated asap).

## Material of interest

### Tools against which it was compared
For validating LIBRA performance we compared it agains other:

- **Integration performance compred to - non-published/available**: [Seurat4](https://github.com/satijalab/seurat).

- **Prediction performance compared to - published/available**: [Seurat3](https://satijalab.org/seurat/articles/integration_mapping.html), [MOFA+](https://biofam.github.io/MOFA2/index.html), [totalVI](https://github.com/YosefLab/scvi-tools), [BABEL](https://github.com/wukevin/babel).

### Required time to run the tool:
All test have been executed on a CPU based server (Intel Corporation Xeon E3-1200 v3 Processor) with reasonable times.

### Required libraries and environmental versioning:
-Please install the following R v3.5.2 libraries before using LIBRA:  
devtools_2.3.1, keras_2.2.5.0, stringr_1.4.0, scclusteval_0.0.0.9000, Seurat_3.2.0, ggplot2_3.3.2, Signac_0.2.1, scater_1.10.1, gridExtra_2.3, ggpubr_0.2.3, biomaRt_2.38.0, scran_1.10.2, cowplot_1.1.1, Matrix_1.2-18, data.table_1.12.4, GenomeInfoDb_1.18.2, EnsDb.Hsapiens.v75_2.99.0, patchwork_1.0.0 and rhdf5_2.26.2.  

-Following libraries requires R v4.0.2:  
Seurat_4.0.0, MOFA2_1.0.1.  

-Please install the following Python libraries:  
scanpy_1.5.0, scvi_0.8.1 and anndata_0.7.5.  


### LIBRA visual workflow:
![workflow.png](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/figures/workflow.png)

### LIBRA code workflow:
Find in "Jupyter_notebook folder" two Jupyter files containing a general example pipeline from LIBRA as well as most relevant QC for measuring performance over integration and prediction outcomes.

### Raw code usage:
- Preprocessing
  - Removing low-quality features and cells
    - Use "1_pre_analysis_Seurat3.R" under Code folder (or your own QC pipeline)
    
- Analysis
  - Normalize and visualize with Seurat for independet samples (or other analysis pipeline that generates at least normalized matrix as output). If more than one sample is present for a given omic first integrate them as usual pipelines requires.
    - Use "2_analysis_Seurat3.R" under Code folder (or your own Normalization pipeline)
    > **Normalization method: Different normalization methods can enhance different LIBRA applications such as - by now - integration or prediction.** As an example Seurat SCT method has shown similar integration performance but a significant decrease on prediction power.

- LIBRA 
  - Use analysis_output.RData as input for LIBRA neural network
    - Use "LIBRA.R" under Code folder for networks training
    - Use "Metrics_LIBRA.R" for additional quality metrix computation
    > **Outputs: Different outputs generated during the training will be stored in the working directory.**

LIBRA was develop in R but feel free to use "rpy2" Python library (https://rpy2.github.io/) for running LIBRA on R snippet through Python console otherwise if your preprocessing was performed using Python but you are interested in running libra in its R implementation you are able to move it to R by using "reticulate" package (https://rstudio.github.io/reticulate/).
