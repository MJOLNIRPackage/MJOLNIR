Full data treatment
^^^^^^^^^^^^^^^^^^^

Full example of the data treatment starting from the instrument definitions provided in the XML file, through generation of normalization table using 8 software pixels, and to data conversion, rebinning and visualization. The file tree for this example is given below where location of data files for normalization and actual data, together with the generated files are shown.

::

    project
    ├── TestData          
    │   └── 1024
    │       └ Magnon_ComponentA3Scan.h5              Magnon scan data in raw format
    └── Tutorials
        └── Full_example.py                          Python script executed to produce data treatment


+-------------------------------------------------------------------+------------------------------------------------------------------+
|      .. literalinclude:: ../../Tutorials/Binning_data.py          |    .. Visualization_E_1_5:                                       |
|         :lines: 4-22,25,28,31-                                    |                                                                  |
|         :language: python                                         |    .. figure:: ../../Tutorials/Visualization_E_1_5.png           |
|         :linenos:                                                 |      :width: 90%                                                 |
|                                                                   |                                                                  |
| Binning of the converted data into Qx and Qy bins of size 0.02 AA |  Cut through data along the energy direction showing Qx and Qy   |
| and energy in 0.1 meV. Intensity is calculated and with the bins  |  for a phonon scan at the energy 1.5 meV.                        |
| passed to the visualizer.                                         |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
|  .. Visualization_Qx_119:                                         |  .. Visualization_Qy_1_84:                                       |
|                                                                   |                                                                  |
|  .. figure:: ../../Tutorials/Visualization_Qx_1_32.png            |  .. figure:: ../../Tutorials/Visualization_Qy_m0_05.png          |
|    :width: 90%                                                    |    :width: 90%                                                   |
|                                                                   |                                                                  |
| Cut of data along the Qx direction.                               | Cut of data along the Qy direction.                              |
|                                                                   |                                                                  |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
