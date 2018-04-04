.. :DataModule:

Data Module
===========

.. currentmodule:: Data



.. autosummary::
   :nosignatures:

    DataSet.DataSet
    DataSet.binData3D
    DataSet.calculateGrid3D
    Viewer3D.Viewer3D

.. automodule:: Data


Data Set Object
^^^^^^^^^^^^^^^

Object to take care of all data conversion and treatment taking it from raw hdf5 files obtained at the instrument into rebinned data sets converted to S(q,omega). 

.. automodule:: DataSet

.. _DataSet:

.. autoclass:: DataSet
    :members:



Data Set bin data
^^^^^^^^^^^^^^^^^
.. autofunction:: binData3D


Data Set calculated 3D grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: calculateGrid3D



Data viewing tool
^^^^^^^^^^^^^^^^^

Visualization tool designed to deal with the difficulties of handling 3D data. 


.. automodule:: Viewer3D

.. _Viewer3D:

.. autoclass:: Viewer3D
    :members:



