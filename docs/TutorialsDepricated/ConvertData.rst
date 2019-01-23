Convert data to Q and Omega
^^^^^^^^^^^^^^^^^^^^^^^^^^^
With the above normalization table created, on can easily convert raw data files using the method in the :ref:`DataSet<DataSet>` called ConvertDatafile as 

.. literalinclude:: ../../TutorialsDepricated/Convert_Data.py
    :lines: 4-
    :language: python 
    :linenos:  

The code then converts the scan files, be it either Ei, A4, or A3 scan, and saves it into a new HDF files following the Nexus NXsqom data convention. This is followed in order to facilitate easy interfacing with other software later used. The location of the converted file is the same as the original but the new file has the ending .nxs. Furthermore, in order to store all information and ensure that no data is lost, all information already present in the raw data file is copied into the new. This also include the original raw data counts.



