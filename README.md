# LIBRA
*Machine Translation between paired Single-Cell Multi-Omics Data*

This repository contains:
- [Online data](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/data/) employed on LIBRA manuscript.
- [Seurat code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/) employed to analyze and visualize integration analysis performed by Seurat.
- [LIBRA code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/) to analyze and visualize paired Single-cell multi-omics integration and prediction analysis.
- [LIBRA metrics](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/) are also available for quantifying outputs quality base on PPJI preservation measurement.

For more information , check out the [online manuscript](https://www.biorxiv.org/content/10.1101/2021.01.27.428400v1) currently at biorxiv repository (will be updated asap).

# Installation:

To run LIBRA pipeline or any other metric generated in the manuscript the following enviromental settings are required:

- Please run [R3_requirements.R](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/files/R3_requirements.R) under Rv3.5.2 or higher Rv3 for automatically install all dependencies required before using LIBRA.

- Please run [R4_requirements.R](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/files/R4_requirements.R) under Rv4.0.3 or higher Rv4 for automatically install all dependencies required before using LIBRA.

-Following libraries are **not supported under Rv3 enviroment: Seurat_4.0.0, MOFA2_1.0.1**.  

-Please install the following **Python libraries**: scanpy_1.5.0, **scvi_0.8.1 ([for totalVI](https://docs.scvi-tools.org/en/stable/installation.html))** and anndata_0.7.5. 

# Usage:

- **LIBRA pipeline is made easy** to be run **especially for any Seurat package user**. 
- The code is executed/stored in **Seurat R objects**, this allows the user to **benefit from the long ecosystem of functions and structures present in Seurat**, working under LIBRA modeling. 
- **Either Seurat 3 in Rv3 enviroment or Seurat 4 in Rv4** enviroment can by used **by hand of LIBRA**.

## Start running LIBRA

### Basic vignettes:
- Model training and integration/prediction [vignette](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/vignettes/Jupyter_notebook/LIBRA_main_pipeline_v1.0.1.ipynb) for a quick example. 
- PPJI preservation metric computation [vignette](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/vignettes/Jupyter_notebook/LIBRA_ppji_metric_v1.0.1.ipynb) for a quick example. 

### Vignettes repository:
- All vignettes can be found [here](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/vignettes/).

# Material of interest

### Tools against which it was compared
For validating LIBRA performance we compared it agains other:

- **Integration performance compred to - non-published/available**: [Seurat4](https://github.com/satijalab/seurat).

- **Prediction performance compared to - published/available**: [Seurat3](https://satijalab.org/seurat/articles/integration_mapping.html), [MOFA+](https://biofam.github.io/MOFA2/index.html), [totalVI](https://github.com/YosefLab/scvi-tools), [BABEL](https://github.com/wukevin/babel).

### LIBRA visual workflow:
![workflow.png](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/figures/workflow.png)

> **LIBRA was develop in R but feel free to use "rpy2" Python library (https://rpy2.github.io/) for running LIBRA on R snippet through Python console otherwise if your preprocessing was performed using Python but you are interested in running libra in its R implementation you are able to move it to R by using "reticulate" package (https://rstudio.github.io/reticulate/).**
