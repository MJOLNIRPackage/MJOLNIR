.. warning::
    This Page is not up to date!

Plotting all pixels binned in Q space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Instead of performing rebinning of the measured datapoints in to some sort of regular grid, one could instead try to create a grid that fits the data measured. This is the basis idea behind the plotA3A4 method. It takes a list of files, creates a commong tesselation using the voronoi method in Q coordinates. 

 .. literalinclude:: ../../Tutorials/PlotQPatches.py
     :lines: 4-  
     :language: python
     :linenos:

As shown above one can provide an axis or a list of axes into which the plot is to be made. This is especially usefull if combined with the RLU axis method calculated in the DataSet object as one then gets the plot in reciprocal lattice units directly. 

.. warning::
    However, this method is really slow and takes approximately 3.5 minutes when combining 2 files and 8 planes.......

.. warning::
    This method does not work fully yet due to inconsistensies in the calibration file of the data.

+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig1:                                                |.. _PlotA3A4_fig2:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/QPatches/000.png                      |.. figure:: ../../Tutorials/QPatches/001.png                      |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 0 through 7 into one plot. This figure is       |Binning of planes 8 through 15 into one plot.                     |
|completely green due to no intensity measured during simulation.  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig3:                                                |.. _PlotA3A4_fig4:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/QPatches/002.png                      |.. figure:: ../../Tutorials/QPatches/003.png                      |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 16 through 23 into one plot.                    |Binning of planes 24 through 31 into one plot.                    |
|                                                                  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig5:                                                |.. _PlotA3A4_fig6:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/QPatches/004.png                      |.. figure:: ../../Tutorials/QPatches/005.png                      |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 32 through 39 into one plot.                    |Binning of planes 40 through 47 into one plot.                    |
|                                                                  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+
|.. _PlotA3A4_fig7:                                                |.. _PlotA3A4_fig8:                                                |
|                                                                  |                                                                  |
|.. figure:: ../../Tutorials/QPatches/006.png                      |.. figure:: ../../Tutorials/QPatches/007.png                      |
|   :width: 95%                                                    |   :width: 95%                                                    |
|                                                                  |                                                                  |
|Binning of planes 48 through 55 into one plot.                    |Binning of planes 56 through 63 into one plot.                    |
|                                                                  |                                                                  |
+------------------------------------------------------------------+------------------------------------------------------------------+

