Quick plotting of data
^^^^^^^^^^^^^^^^^^^^^^
For the initial overview of what has been measured and for viewing some features in the data, a quick 3D viewer can be initiated. This tutorial shows how to initialize it but not how the interactive parts work; for this see `Viewer3D documentation <../../Data/Gui.html#viewer3d>`_ or the `dedicated tutorial <Viewer3D.html>`_. What is done when using the viewer is that all data from all provided files are binned in 3D pixels (voxels) of equal size independent of density of measured points. This does of course not reflect the actual data quality but gives quick insights into the data.

.. code-block:: python
   :linenos:

   from MJOLNIR.Data import DataSet
   import numpy as np
   
   # Location of raw data file. Can also be given as a single string
   fileName = ['/Path/To/Data/camea2018n000136.hdf','/Path/To/Data/camea2018n000137.hdf']
   # Define the DataSet object and provide file(s)
   ds = DataSet.DataSet(dataFiles=fileName)
   # Run the converter. This automatically generates nxs-file(s) in location
   # of original file. Can be changed with key-word arguments.
   ds.convertDataFile(saveFile=False)
   
   # With saveFile = False, converted files are not saved to disk but stays
   # only in RAM. 
   for df in ds:
   mask = np.zeros_like(df.I.data) # Define mask, see FAQ for explanation
   mask[:,:,:3]=True
   df.mask = np.logical_or(df.mask,mask)
   
   
   # Plotting data quickly in equi-sized voxels can be done by
   Viewer = ds.View3D(0.03,0.03,0.05)
   # Above, all data is binned in voxels of size 0.03/AA, 0.03/AA, and 0.05 meV.
   # Automatically, data is plotted in reciprocal lattice as provided by the
   # UB matrix saved in the data file. Alternatively, one can plot in 'raw'
   # coordinates (orthogonal and in units of 1/AA) by issuing rlu=False above.
   
   Viewer.caxis=(1e-7,2.5e-6)
   # Without any intervention data is usually plotted on a useless colour scale.
   # This is countered by specifying min and max for colour by the call above.
   # Alternatively, one can provide this as input to View3D 
   
   # Set title of plot
   Viewer.ax.set_title('Magnon ComponentA3Scan')
   
   # Change to plane number 44. As default, the viewer cuts for contant energy
   # and with a binning of 8 pixels/segment around 64 planes are available.
   Viewer.setPlane(44)
   
   # Get figure of Viewer
   fig = Viewer.ax.get_figure()
   
   fig.savefig('figure0.png',format='png')
   
   # Change cut to be along first axis (in this case [1,0,0]). This is in general
   # the vector in reciprocal space closest to [1,0,0]. In reciprocal lattice,
   # axis 1 is always orthogonal on axis 0. This is important to remember when
   # dealing with non-orthogonal crystals or non-trivial alignments.
   
   Viewer.setProjection(0)
   
   Viewer.setPlane(31)
   fig.savefig('figure1.png',format='png')
   
   
   Viewer.setProjection(1)
   
   Viewer.setPlane(61)
   fig.savefig('figure2.png',format='png')
   

As indicated in the code there are many further possibilities for making use of the 3Dviewer but for these, see the `Advanced View3D tutorial <../Advanced/View3D.html>`_. 
The plots generated from the above code is shown below: 

.. figure:: EnergyPlaneViewer.png
  :width: 50%
  :align: center

View of constant energy plane for YMnO3 with the dispersion shown

.. figure:: HPlaneViewer.png
  :width: 50%
  :align: center

Change of axis to be along [H,0,0] with a cut through the dispersion. 

.. figure:: OrthogonalPlaneViewer.png
  :width: 50%
  :align: center

Choosing axis to be 1, one gets the QE plane orthogonal to the [H,0,0] vector, which in a hexagonal system is [-H,2K,0] if [H,0,0] and [0,K,0] are in the plane.
