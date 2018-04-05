Full data treatment
^^^^^^^^^^^^^^^^^^^

Full example of the data treatment starting from the instrument definitions provided in the XML file, through generation of normalization table using 8 software pixels, and to data conversion, rebinning and visualization. The file tree for this example is given below where location of data files for normalization and actual data, together with the generated files are shown.

::

    project
    ├── TestData          
    │   ├── CAMEA_Full_2.xml                         Instrument description
    │   ├── VanNormalization.h5                      Normalization data in raw format
    │   └── cameasim2018n000005.h5                   Phonon scan data in raw format
    └── Tutorials
        ├── EnergyNormalization_8.calib              Generated 8 pixel normalization calibration table
        ├── cameasim2018n000005.nxs                  Converted phonon scan data in NXqom format
        └── Full_example.py                          Python script executed to produce data treatment


+-------------------------------------------------------------------+------------------------------------------------------------------+
|      .. literalinclude:: ../../Tutorials/Full_example.py          |    .. Visualization_E_546:                                       |
|         :lines: 4-                                                |                                                                  |
|         :language: python                                         |    .. figure:: ../../Tutorials/Visualization_E_546.png           |
|         :linenos:                                                 |      :width: 90%                                                 |
|                                                                   |                                                                  |
| Binning of the converted data into Qx and Qy bins of size 0.05AA  |  Cut through data along the energy direction showing Qx and Qy   |
| and energy in 0.2 meV. Intensity is calculated and with the bins  |  for a phonon scan at the energy 5.46 meV.                       |
| passed to the visualizer.                                         |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
|  .. Visualization_Qx_119:                                         |  .. SimpleInstrumentEf:                                          |
|                                                                   |                                                                  |
|  .. figure:: ../../Tutorials/Visualization_Qx_119.png             |  .. figure:: ../../Tutorials/Visualization_Qy_089.png            |
|    :width: 90%                                                    |    :width: 90%                                                   |
|                                                                   |                                                                  |
| Cut of data along the Qx direction.                               | Cut of data along the Qx direction.                              |
|                                                                   |                                                                  |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
