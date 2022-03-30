LIBRA - Machine Translation between paired <img src="gaf/figures/LIBRA_icon_2.png" width="181px" align="right" />  
Single-Cell Multi-Omics Data 
===========
This repository contains the [LIBRA code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/LIBRA_code/) and [online data](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/data/) used for Single-cell multi-omics integration and prediction analysis employed on [LIBRA manuscript](https://www.biorxiv.org/content/10.1101/2021.01.27.428400v1). [Libra metrics](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/LIBRA_code/) are also available for quantifying outputs quality as well as novel PPJI preservation measurement. [Seurat code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/Seurat_code/) employed to analyze LIBRA input omics as well as for clustering and visualization pipelines are providen.

- [Summary](#summary)
- [Prerequisites](#prerequisites)
- [Datasets](#datasets)
- [Usage](#usage)
- [Getting Started LIBRA](#getting-started-libra)
- [Material of interest](#material-of-interest)

# Summary
LIBRA is a deep learning model that is designed for Single-cell multi-omics integration and prediction. LIBRA performs this by using an unbalance Autoencoder which learns a shared low-dimensional embedding from both experiment omics, combining each sample's uniqueness for generating a enriched representation of integrated data respect to the original experiment independent data. This tool has been develop in [R code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/R/LIBRA_code/). Next, fine-tune LIBRA tool has been develop for paralellize training of LIBRA models using a grid structure for selecting optimal hyperparameters in a automatic way excluding the requirement of doing this by users saving considerable time. This tool is providen in [Python code](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/Python/LIBRA_fine_tune_code/).

For further details, please refer to the [online manuscript](https://www.biorxiv.org/content/10.1101/2021.01.27.428400v1) currently at biorxiv repository (will be updated asap).

# Prerequisites
 
To run LIBRA pipeline or any other metric generated in the manuscript the following environmental settings are required:

- Run [R3_requirements.R](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/files/R3_requirements.R) under R 3.5.2 or higher R 3.X.X for automatically install all dependencies required before using LIBRA.

- Run [R4_requirements.R](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/files/R4_requirements.R) under R 4.0.3 or higher R 4.X.X for automatically install all dependencies required before using LIBRA.

-Following libraries are **not supported under R 3.X.X environment: Seurat_4.0.0, MOFA2_1.0.1**.  

To run LIBRA fine-tune pipeline generated in the manuscript the following environmental settings are required:

-Please install the following **Python libraries** under Python v3.7.1 or higher: scanpy_1.5.0, **scvi_0.8.1 ([for totalVI](https://docs.scvi-tools.org/en/stable/installation.html))**, anndata_0.7.5, pandas_1.3.4, numpy_1.18.1, scipy_1.7.1, keras_2.7.0 and multiprocessing_2_6_2.

# Datasets
| LIBRA name | GSE link | Modalities | Technology | Genomic ref used | Download sparse matrix |
|    :---:    |    :---:    |    :---:    |    :---:    |    :---:    |    :---:    |
| DataSet1 | [GSE126074](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074) | scRNAseq + scATACseq | SNARE-seq | [Mus_musculus.GRCm38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/c9b87f4ac1d1c030e128) and [ATAC](https://figshare.com/s/9ff9ea93a2108478bb36) |
| DataSet2 | [GSE128639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) | scRNAseq + scADT | CITE-seq | [Homo_sapiens.GRCh38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/5f5cfa6fda4ae3512c0d) and [ADT](https://figshare.com/s/5e34cd80455398855ad8) |
| DataSet3 | [GSE130399](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) | scRNAseq + scATACseq | Paired-seq | [Mus_musculus.GRCm38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/a1f4a5ef0735d1b4167d) and [ATAC](https://figshare.com/s/80d9b9d84ada526668a6) |
| DataSet4 | [GSE140203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) | scRNAseq + scATACseq | SHARE-seq | [Mus_musculus.GRCm38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/71312a335649b04972b8) and [ATAC](https://figshare.com/s/0b581450cd6e1f8fb64c) |
| DataSet5 | [10X Genomics](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k) | scRNAseq + scATACseq | 10X multiome | [Homo_sapiens.GRCh38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/90b237227f0cc07d075d) and [ATAC](https://figshare.com/s/4086bce6032f6a206a13) |
| DataSet6 | [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) | scRNAseq + scATACseq | 10X multiome | [Homo_sapiens.GRCh38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/134562c3ec74a3a50c84) and [ATAC](https://figshare.com/s/378a630ec9c6ddadf4f5) |
| DataSet7 | [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) | scRNAseq + scADT | CITE-seq | [Homo_sapiens.GRCh38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/41bdbfe7479e9729c800) and [ADT](https://figshare.com/s/975cd8a5bbc57c8d2c8c) |
| DataSet8 | [GSE109262](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) | scRNAseq + scATACseq | scNMT-seq | [Mus_musculus.GRCm38 Ver: 3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#) | [RNA](https://figshare.com/s/4a158e6d243bcd45a171) and [ATAC](https://figshare.com/s/f13136b52f3b387d1a66) |

# Usage

- **LIBRA pipeline is made easy** to be run especially for any Seurat package user. 
- The code is executed/stored in **Seurat R objects**, this allows the user to **benefit from the long ecosystem of functions and structures present in Seurat**, working under LIBRA modeling. 
- Either **Seurat3** in R 3.X.X environment or **Seurat4** in R 4.X.X environment can be used **by hand of LIBRA**.
- The **valid input for LIBRA** is any pair of omic matrices assigning the cell information in the rows and the feature information in the columns.

# Getting Started LIBRA

### Basic vignettes:
- Model training and integration/prediction [vignette](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/vignettes/Jupyter_notebook/LIBRA_main_pipeline_v1.0.1.ipynb) for a quick example. 
- PPJI preservation metric computation [vignette](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/vignettes/Jupyter_notebook/LIBRA_ppji_metric_v1.0.1.ipynb) for a quick example. 

### Vignettes repository:
- All vignettes can be found [here](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/vignettes/).

# Material of interest

### LIBRA benchmarking comparison:
For validating LIBRA performance we compared it against other:

- Integration performance compared to - non-published/available: [Seurat4](https://github.com/satijalab/seurat).

- Prediction performance compared to - published/available: [Seurat3](https://satijalab.org/seurat/articles/integration_mapping.html), [MOFA+](https://biofam.github.io/MOFA2/index.html), [totalVI](https://github.com/YosefLab/scvi-tools), [BABEL](https://github.com/wukevin/babel).

### LIBRA visual workflow:
![workflow.png](https://github.com/TranslationalBioinformaticsUnit/LIBRA/blob/main/gaf/figures/workflow.png)

### Versatility:
> **LIBRA was develop in R but feel free to use "rpy2" Python library (https://rpy2.github.io/) for running LIBRA on R snippet through Python console otherwise if your preprocessing was performed using Python but you are interested in running libra in its R implementation you are able to move it to R by using "reticulate" package (https://rstudio.github.io/reticulate/).**
