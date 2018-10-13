Bin data and visualize
^^^^^^^^^^^^^^^^^^^^^^
Having converted the data into the Nexus NXsqom format, one wants to both rebin the data and visualize it. As different detector pixels covers different positions in reciprocal space, one needs to rebin the data in order to avoid large areas of no data. This is done using the method in the :ref:`DataSet<DataSet>` called binData3D. As input one needs to provide the step size in the x, y, and z directions, the 3D position and intensity. Furthermore, normalization and monitor count can be specified in order to also bin these. Returned is the rebinned data together with normalization and monitor count, if applicable, and the bins used.

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

The bins and calculated intensity is then passed on to the Viewer3D object, that generates a matplotlib figure. This plot is made interactive through the slider in the bottom, that shows the current position along the axis as well as the value of this. By pressing the up/down arrows (or +/- buttons as well as scrolling) one can change this value and thus investigate the third dimension. By default the energy direction is chosen but it can be changed by pressing the 0, 1, or 2 numerical buttons. This is to be understood as which direction is to be sliced. Further, by clicking on the plot, the x and y value corresponding to the point is printed to the terminal. This is then intented to be used in further data treatment when cuts are needed. 

