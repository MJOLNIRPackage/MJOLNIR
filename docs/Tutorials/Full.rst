Full data treatment
^^^^^^^^^^^^^^^^^^^

Full example of the data treatment starting from the instrument definitions provided in the XML file, through generation of normalization table using 8 software pixels, and to data conversion, rebinning and visualization. 


+-------------------------------------------------------------------+------------------------------------------------------------------+
|      .. literalinclude:: ../../Tutorials/Full_example.py          |    .. Visualization_E_1_5:                                       |
|         :lines: 4-                                                |                                                                  |
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
