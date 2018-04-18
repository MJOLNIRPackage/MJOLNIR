Full data treatment without Viewer3D
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Full example of a data treatment workflow starting from raw h5 files and ending in a plot. The treatment is done in four steps signified in the code by the comments:

    - Convert raw data from h5 to nxs (NXSqom) data format

    - Open converted files and extract intensities and measurement positions from them

    - Bin the data using polar binning and calculated the intensity for each point

    - Plot an energy slice of the data set in a regular Matplotlib figure



 .. literalinclude:: ../../Tutorials/Full_example_without_Viewer3D.py
    :lines: 1-             
    :language: python      
    :linenos:


In the above example, there are two key points in the treatment; first the binning in polar coordinates which takes advantage of the measurement positions of the instrument setup (equidistant A3 or theta steps), second that the binned data from the binning method is of e.g. shape (20,23,6) while the position arrays contains the bin edges and have thus the shape (21,24,7). This makes it easy to perform 2D cuts of the data along the primal axis (radius, angle and energy) by simply slicing data using regular slicing tools. 

The arguments vmin and vmax controls the colorbar and thus the colors corresponding to different intensities. 

.. note::
    The way that Matplotlib currently has implimented the pcolormesh method for plotting 2D data requires the position matrices and data matrix to be transposed before plotting. Further, applying the 'gouraud' interpolation scheme provided by the method, one needs to give as arguments the centers of the binned data instead of the bin edges. This can simply be given by defining e.g. a new Qx by QX = 0.5*(Qx[:-1,:-1,Eslice]+Qx[1:,1:,Eslice])
