Advanced View3D tutorial
^^^^^^^^^^^^^^^^^^^^^^^^
Assuming that the Quick visualization from the Quick tutorials is understood, this tutorial seeks to introduce more advanced features for the 3D viewer object.

.. code-block:: python
   :linenos:

   from MJOLNIR.Data import DataSet
   # Load and convert data
   fileName = ['/Path/To/Data/camea2018n000136.hdf','/home/lass/Dropbox/PhD/CAMEAData/camea2018n000137.hdf']
   ds = DataSet.DataSet(dataFiles=fileName)
   ds.convertDataFile(saveFile=False)
   

