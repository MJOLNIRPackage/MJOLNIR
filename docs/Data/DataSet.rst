.. :DataSet:

Data Set
========

.. currentmodule:: Data.DataSet

The DataSet object is the interface between the data files and the data treatment and visualziation. It is both responsible for the conversion of raw '.h5'-files into '.nxs'-files as well as plotting these. Extracting values from this object results in a list of values where the first dimension is determined from the number of data files provided.

.. autosummary::
   :nosignatures:

    DataSet.DataSet
    DataSet.DataSet.convertDataFile
    DataSet.DataSet.cut1D
    DataSet.DataSet.plotCut1D
    DataSet.DataSet.cutQE
    DataSet.DataSet.plotCutQE
    DataSet.DataSet.cutPowder
    DataSet.DataSet.plotCutPowder
    
    DataSet.DataSet.plotQPlane

    DataSet.DataSet.cutQELine
    DataSet.DataSet.plotCutQELine



.. automodule:: Data
   :members:



DataSet Object and Methods
--------------------------

Object to take care of all data conversion and treatment taking it from raw hdf5 files obtained at the instrument into re-binned data sets converted to S(q,omega). 

.. automodule:: DataSet
   :members:

.. _DataSet:

.. autoclass:: DataSet
    :members:




