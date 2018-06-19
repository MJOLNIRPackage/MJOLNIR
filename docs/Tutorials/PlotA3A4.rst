
Plotting all pixels binned in A3 and A4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Instead of performing rebinning of the measured datapoints in to some sort of regular grid, one could instead try to create a grid that fits the data measured. This is the basis idea behind the plotA3A4 method. It takes a list of files, creates a commong tesselation using the voronoi method in A3-A4 coordinates and maps this into Q-space. 

 .. literalinclude:: ../../Tutorials/PlotA3A4.py
     :lines: 4-  
     :language: python
     :linenos:

As shown above one can provide an axis or a list of axes into which the plot is to be made. This is especially usefull if combined with the RLU axis method calculated in the DataSet object as one then gets the plot in reciprocal lattice units directly.

+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig1:                                                |.. _PlotA3A4_fig2:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/A3A4/000.png                          |.. figure:: ../../Tutorials/A3A4/001.png                          |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 0 through 7 into one plot. This figure is       |Binning of planes 8 through 15 into one plot.                     |
|completely green due to no intensity measured during simulation.  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig3:                                                |.. _PlotA3A4_fig4:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/A3A4/002.png                          |.. figure:: ../../Tutorials/A3A4/003.png                          |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 16 through 23 into one plot.                    |Binning of planes 24 through 31 into one plot.                    |
|                                                                  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig5:                                                |.. _PlotA3A4_fig6:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/A3A4/004.png                          |.. figure:: ../../Tutorials/A3A4/005.png                          |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 32 through 39 into one plot.                    |Binning of planes 40 through 47 into one plot.                    |
|                                                                  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig7:                                                |.. _PlotA3A4_fig8:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/A3A4/006.png                          |.. figure:: ../../Tutorials/A3A4/007.png                          |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 48 through 55 into one plot.                    |Binning of planes 56 through 63 into one plot.                    |
|                                                                  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+

