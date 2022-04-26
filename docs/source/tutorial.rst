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

**Example (input files are in wd)**:

.. code-block:: python

    par = sc_libra.load_data("ATAC","RNA", dataset1='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2.h5ad', dataset2='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad')
    
If input datasets are in a different location than working directory, *dataset1_path* and *dataset2_path* can be used (dataset2_path will be equal to dataset1_path if not setted).

**Example (input files are in other directory than wd)**:

.. code-block:: python

    par = sc_libra.load_data("ATAC","RNA", dataset1_path='/both_are_in_this_different_location', dataset1='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod2.h5ad', dataset2='openproblems_bmmc_multiome_phase2.censor_dataset.output_mod1.h5ad')

If 10X folder contains the input data desired, specify only the path where files are and they will be automatically loaded.

**Example (10X folder as input)**:

.. code-block:: python

    par = sc_libra.load_data("ATAC","RNA", dataset1_path='/location_to_10x_folder_for_input_omic_ATAC', dataset2_path='/location_to_10x_folder_for_output_omic_RNA')

As a result output (par in these examples) will contain a dictionary such as: 
   - {**omic_1_name**: pandas.dataframe.omic1, **omic_2_name**: pandas.dataframe.omic2}.


Training LIBRA model
--------------------

LIBRA can run in many different ways using the *libra* function. This step uses the previously generated variable as input (in this example, par), but any other generated dictionary by user that satisfies the above requirements can be used as input. 

The most basic way is at following example shows. This will train the LIBRA model with default parameters finding a good balance between prediction/integration performance. Will generate integration output file containing latent space for each cell and sotore it in automatically generated tree of directories. Model will also be stored in .hdf5 format.

**Example (default use)**:

.. code-block:: python

    output_data = sc_libra.libra(par)

LIBRA can also be used for training a bunch of models for bosting performance on one of the main tasks over the other (prediction/integration). To this aim a grid of parameters will be used generating hundreds of models and storing the outputs following the same default squema. A custom grid can also be used if desired by user.

**Example (bossting one task over the other)**:

.. code-block:: python

    output_data = sc_libra.libra(par, training_mode = 'fine_tune_prediction') #For prediction best model finding
    output_data = sc_libra.libra(par, training_mode = 'fine_tune_integration') #For prediction best model finding
    output_data = sc_libra.libra(par, training_mode = 'custom') #For custom grid user
 
Extra parameters can by added to the function for example *n_top_genes*. In the case of containing a omic named as "RNA" *libra* function will filter gen space to contain only the most 2000 highly variable genes, this is peformed becose in our experiments RNA has prove to provide better performance over LIBRA model when only using HVG. If a different amount of genes is wanted it can be setted as in following example:

**Example (using other amount of genes than 2000 HVG)**:

.. code-block:: python

    output_data = sc_libra.libra(par, n_top_genes = 3000) #For use 3000 number of HVG
    
For bosting speed (if user hardware is sufficient) and extra parameter can be added, *n_jobs*. This parameter setted as default to 1, can be changed to any amount of cores present in users CPU to perform multiple model trainings in paralel. This is designed specifically for other that the default *libra* option where many models will be trained depending on grid selected. This reduces the time required but also requires more RAM memory.

**Example (parallel training for grid based version)**:

.. code-block:: python

    output_data = sc_libra.libra(par, n_jobs=20) #For training 20 models in parallel (your CPU should have at least 20 cores, and enought RAM to handle them in memmory).

All these parameters can be combined for desired task.


Prediction using LIBRA model
----------------------------

If user what to use LIBRA model generated for a prediction task over same or new input dataset, it can be done through this function, *libra_predict* as following example. Either latent of output spaces can be predicted.

**Example (predict over input dataset)**:

.. code-block:: python
    
    model = load_model('/.../LIBRA_outputs/Models/model_n_layers2_n_nodes512_alpha0.3_dropout0.2_batch_size7000_mid_layer10.hdf5')
    input_data = output_data[0].todense() #For predict over input dataset. A novel one can be used here.
    to_predict = 'integrated_space' #For latent space prediction or 'modality_B' for output prediction.
    
    predicted_data = sc_libra.libra_predict(model, input_data, to_predict)

Metrics computation
-------------------







