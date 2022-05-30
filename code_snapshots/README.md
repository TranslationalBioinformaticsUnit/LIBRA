# Prerequisites (for running manuscript codes)
 
To run LIBRA pipeline form R or Python source code or any other metric generated in the manuscript the following environmental settings are required:

- Run [R3_requirements.R](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/files/R3_requirements.R) under R 3.5.2 or higher R 3.X.X for automatically install all dependencies required before using LIBRA.

- Run [R4_requirements.R](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/files/R4_requirements.R) under R 4.0.3 or higher R 4.X.X for automatically install all dependencies required before using LIBRA.

-Following libraries are **not supported under R 3.X.X environment: Seurat_4.0.0, MOFA2_1.0.1**.  

To run LIBRA fine-tune pipeline generated in the manuscript the following environmental settings are required:

-Please install the following **Python libraries** under Python v3.7.1 or higher: scanpy_1.5.0, **scvi_0.8.1 ([for totalVI](https://docs.scvi-tools.org/en/stable/installation.html))**, anndata_0.7.5, pandas_1.3.4, numpy_1.18.1, scipy_1.7.1, keras_2.7.0 and multiprocessing_2_6_2.

### Versatility:
> **LIBRA manuscript code was develop in R but feel free to use "rpy2" Python library (https://rpy2.github.io/) for running LIBRA on R snippet through Python console otherwise if your preprocessing was performed using Python but you are interested in running libra in its R implementation you are able to move it to R by using "reticulate" package (https://rstudio.github.io/reticulate/).**

# Usage (manuscript codes)

- **LIBRA pipeline is made easy** to be run especially for any Seurat package user. 
- The code is executed/stored in **Seurat R objects**, this allows the user to **benefit from the long ecosystem of functions and structures present in Seurat**. 
- Either **Seurat3** in R 3.X.X environment or **Seurat4** in R 4.X.X environment can be used **by hand of LIBRA**.
- The **valid input for LIBRA** is any pair of omic matrices assigning the cell information in the rows and the feature information in the columns.
- Easiest way of running LIBRA analysis is though [sc-Libra](https://pypi.org/project/sc-libra/) python package. Package [documentation](https://sc-libra.readthedocs.io/en/latest/) is online available using "Read the Docs" platform.

# Getting Started with LIBRA (manuscript codes)

### Basic vignettes:
- Model training and integration/prediction [vignette](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/code_snapshots/vignettes/Jupyter_notebook/LIBRA_main_pipeline_v1.0.1.ipynb) for a quick example. 
- PPJI preservation metric computation [vignette](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/code_snapshots/vignettes/Jupyter_notebook/LIBRA_ppji_metric_v1.0.1.ipynb) for a quick example. 

### Vignettes repository:
- All vignettes can be found [here](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/code_snapshots/vignettes/).
