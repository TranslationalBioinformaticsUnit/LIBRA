Quickstart
==========

LIBRA is a modular toolbox and hence easy to use.
All outputs from functions and directorie tree are generated behind scenes and the amount of setup the user needs to do is minimal.

Usually only the following steps are required:

#. Install Python :math:`3.7.0 or above`.
#. Install R :math:`3.5.2 or above`.
#. Install :math:`sc_LIBRA` Python package.
#. Prepare the :math:`enviroment` before importing sc_LIBRA.


Installation
------------

LIBRA is compatible with Python 3.7+, and depends on rpy2, NumPy, SciPy, Pandas and Keras.
All dependencies will by installed when sc_libra is installed through:

.. code-block:: bash

    $ pip install sc_libra


Enviroment preparation
----------------------
LIBRA make use of rpy2 for running some specific R functions. In order to import properly this dependency is mandatory that Python knows where are the R libs for the specific R version used. This is a requirement of rpy2 and should be done, else sc_libra will raise and error when importing it (as it will import all dependecies required as it is imported) This can be done in two steps:

#. Run: :math:`export LD_LIBRARY_PATH="/YOUR_PATH_TO_R_LIBS_HERE:$LD_LIBRARY_PATH"` on console prior to run Python.

Typical locations are: 

- :math:`export LD_LIBRARY_PATH="/opt/R/3.5.2/lib64/R/lib:$LD_LIBRARY_PATH"` (if local installation of R was done) 

- :math:`export LD_LIBRARY_PATH="/usr/lib64/R/:$LD_LIBRARY_PATH"` (if global instalaltion of R was done).


#. Open your Python enviroment and run:

- :math:`import os`

- :math:`os.environ['R_HOME'] = "/opt/R/3.5.2/lib64/R/"`  or "/usr/lib64/R/" or the corresponding one.

Through these two steps possible miss-making errors because incorrect R version pointed can be evade (often if more than one R version is installed).
