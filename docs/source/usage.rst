Usage
=====

.. _installation:

Installation
------------

To use sc-LIBRA, first install it using pip:

.. code-block:: console

   (.venv) $ pip install sc_libra

Using pipeline
----------------

To load input data,
you can use the sc_libra.data_load()`` function:

.. autofunction:: sc_libra.data_load

The omic types parameter should be either ``"omic_1"``, ``"omic_2"``,
or ``"None"``. Otherwise, :py:func:`sc_libra.data_load`
will raise an exception.

.. autoexception:: sc_libra.InvalidKindError

For example:

>>> import sc_libra
>>> sc_libra.load_data("ATAC","RNA")
["ATAC","RNA"]

