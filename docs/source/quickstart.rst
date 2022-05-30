Quickstart
==========

LIBRA is a modular toolbox and hence easy to use.
All outputs from functions and directories tree are generated behind scenes and required user interaction is low.

Follow this steps to a proper installation:

#. Install Python **>=3.7.0**.
#. Install R **>=3.5.2**.
#. Install `sc_libra <https://pypi.org/manage/project/sc-libra/releases/>`_ Python package.
#. **OPTIONAL: Prepare the environment (Only for selecting a specific R version in case many are installed, otherwise avoid this step)**.
#. Importing sc_libra.


Installation
------------

LIBRA is compatible with Python >=3.7, and depends on rpy2, NumPy, SciPy, Pandas and Keras.
All required dependencies will be automatically installed when running:

.. code-block:: bash

    $ pip install sc_libra


(Optional) Environment preparation
----------------------
LIBRA makes use of rpy2 for running some specific R functions. In order to import properly this dependency is mandatory that Python knows where are the R libs for the specific R version used. This is a requirement of rpy2 and should be done, else sc_libra will raise an error when importing it (as it will import all dependencies required as it is imported) This can be done in two steps:

1. Run on console prior to run Python.

.. code-block:: bash

    $ #Typical locations are: 
    $ # export LD_LIBRARY_PATH="/opt/R/3.5.2/lib64/R/lib:$LD_LIBRARY_PATH" (if local installation of R was done) 
    $ # export LD_LIBRARY_PATH="/usr/lib64/R/:$LD_LIBRARY_PATH" (if global installaltion of R was done)
    
    $ export LD_LIBRARY_PATH="/YOUR_PATH_TO_R_LIBS_HERE:$LD_LIBRARY_PATH"

2. Open your Python environment and run:

.. code-block:: bash

    $ import os
    $ # Typical locations are:
    $ # export os.environ['R_HOME'] = "/opt/R/3.5.2/lib64/R/" (if local installation of R was done)
    $ # export os.environ['R_HOME'] = "/usr/lib64/R/" (if global installation of R was done)
    
    $ os.environ['R_HOME'] = "/YOUR_POATH_TO_R_HERE"
