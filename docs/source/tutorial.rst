Tutorial
==========

LIBRA provides four main functions that coverages all required aspects of the proposed pipeline.

Enviroment preparation
------------

First dependencies should be loaded by setting R path and then importing **sc_libra** module (at this moment rest of required dependencies will also be loaded). In order to have an ordered working directory is recommended to set the directory path desired where the tree of folders and files will be saved during the different functions execution. The tree of folders will be generated under working directory path automatically.

.. code-block:: python

    import os
    os.environ['R_HOME'] = "/opt/R/3.5.2/lib64/R/"
    import sc_libra
    os.chdir("/desired_working_directory_path")

Loading input data 
------------------

A universal loading data is provided for easy data loading into Python enviroment. Input file formats will be automatically detected among: **AnnData(".h5ad"), sparse matrix(".mtx",".txt"), coma separated matrix(.csv), 10XGenomics(.mtx folder) and 10XGenomics(.h5 file)**. If format is different data should be loaded by user with corresponding command and harmonized as this function does. 

Datasets should be loaded in order been dataset1 the input modality to LIBRA model and dataset2 the output modality. **We strongly recomend to set RNA as the second modality**. By now only specific actions are performed over RNA omic so it is **important to set "RNA" name to RNA omic input dataset**, in oder cases desired name can by used. If both input datasets are in working directory path, no paths are required by *load_data* function.

Example (input files are in wd):

.. code-block:: python

    par = sc_libra.load_data("ATAC","RNA", dataset1='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2.h5ad', dataset2='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad')
    
If input datasets are in a different location than working directory, *dataset1_path* and *dataset2_path* can be used (dataset2_path will be equal to dataset1_path if not setted).

Example (input files are in other directory than wd):

.. code-block:: python

    par = load_data("ATAC","RNA", dataset1_path='/both_are_in_this_different_location', dataset1='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2.h5ad', dataset2='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad')

If 10X folder contains the input data desired, specify only the path where files are and they will be automatically loaded.

Example (10X folder as input):

.. code-block:: python

    par = load_data("ATAC","RNA", dataset1_path='/location_to_10x_folder_for_input_omic_ATAC', dataset2_path='/location_to_10x_folder_for_output_omic_RNA')

As a result output (par in these examples) will contain a dictionary such as: 
   - {**omic_1_name**: pandas.dataframe.omic1, **omic_2_name**: pandas.dataframe.omic2}.

Training LIBRA model
--------------------
