FAQ
===

A collection of questions and answers meant to guide the user(s) of the MJOLNIR package through some of the most common pitfalls.


My RLU-axis has its tick marks positions stupidly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This 'feature' of the RLU axis is unfortunately linked to the current stage of the development of the `GridHelperCurveLinear <https://matplotlib.org/api/_as_gen/mpl_toolkits.axisartist.grid_helper_curvelinear.GridHelperCurveLinear.html#mpl_toolkits.axisartist.grid_helper_curvelinear.GridHelperCurveLinear>`_.
Because of the current lack of customization for this in-progress work some attributes cannot be changed. However, when generating an new RLUaxis 
one can specify the number of tick marks on the x and y axis. This number stays constant  throughout the lifetime of that axis, but a new can be created with a 
different number if needed. Further, specifying only the number of y-axis tick marks the program calculates a suitable number along the x-axis.


Every time I run my script, I am bombarded with 'The file xx exists already. Old file will be renamed to xx.'-warnings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The warning originates from the fact that a data file with the same name has been converted previously in the same directory, 
and the program then tells you that the old conversion will be moved to a new file with the new name. It only serves as to try 
to prevent you from saving a conversion on top of another involuntarily. To  silence this annoyance, convert the data files only 
once or use the 'saveFile=False' flag in the DataSet.convertDataFile() method to not save the converted data.


I do not want to make use of the features in MJOLNIR but only convert and extract my data. How do I do this?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is definitely possible to convert and extract the data from MJOLNIR without many troubles. Below is a code example:

.. code-block:: python
   :linenos:

    from MJOLNIR.Data import DataSet
    fileName = ['/Path/To/Data/camea2018n000136.hdf','/Path/To/Data/camea2018n000137.hdf']
    ds = DataSet.DataSet(dataFiles=fileName)
    ds.convertDataFile(binning=3,saveFile=False)

    Int = np.divide(ds.I,ds.Monitor*ds.Norm) # Calcualte floating point intensity

    data = np.array([ds.h.extractData(),ds.k.extractData(),ds.l.extractData(),Int.flatten()])

    saveFile = '/Path/To/Data/Data.csv'
    np.savetxt(saveFile,data,delimiter=',')

However, a small warning. The data file created from this can be of several hundred MB.



My data has annoying lines for constant energy that looks like background
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When measuring data on CAMEA and utilizing the prismatic concept one can split the nominal eight final energies into a maximum of 64 energies. This is performed in the conversion and is encoded in the normalization
files saved in each data file. Due to difference in sensitivity of these prismatic pixels any uniform background on the detector tube will be inflated at the points of low sensitivity resulting in lines of 
higher scattering intensity for constant energies. These then show up as lines, both in 2D plots but more importantly especially in energy cuts. Currently. there are two solutions to be employed, one for during experiments 
and one for post experiment data treatment if first possibility is not done:

- When performing a measurement, interlace all A3/sample rotation scans with additional scans with an Ei slightly changed.

- Perform a powder background subtraction as explained in the `Background Subtraction Tutorial <Tutorials/Advanced/backgroundSubtraction.html>`_.

The amount of change needed in :math:`E_i` to counter the stripes can be difficult to figure out and depends on the question at hand. However, as a standard approach 
one has to shift with half the mean difference in final energy. This ensures that the parts of the detectors with low sensitivity will be covered by parts of the detector tube 
with better sensitivity. Taking this half average, one finds the value: 0.127 meV. However, if one is applying a high resolution measurement where only a subset
of the final energies are useful (on the energy-loss side) and one would move :math:`E_i` by half the average of the final energies in use. Most often one uses
the (approximate) 3.6, 3.4, and 3.2 meV final energies where a shift of ~0.094 meV. Lastly, one can also add multiple shifts in :math:`E_i` whith sizes of 0.064, 0.127, and 0.191 meV, 
corresponding to 0.5, 1.0, and 1.5 times the energy shift.




.. table:: Mean Final Energies (`Commissioning Report <https://aip.scitation.org/doi/10.1063/5.0128226>`_):
   :widths: auto
   :align: center

   +-----------+------------------+
   | Ef [meV]  | Difference [meV] |
   +===========+==================+
   |  3.200    |  N/A             |
   +-----------+------------------+
   |  3.382    |  0.182           |
   +-----------+------------------+
   |  3.574    |  0.192           |
   +-----------+------------------+
   |  3.787    |  0.213           |
   +-----------+------------------+
   |  4.035    |  0.248           |
   +-----------+------------------+
   |  4.312    |  0.277           |
   +-----------+------------------+
   |  4.631    |  0.319           |
   +-----------+------------------+
   |  4.987    |  0.356           |
   +-----------+------------------+



I aligned my crystal slightly off when setting up the experiment. Can I correct for an offset in A3?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is definitely possible through the code below, where 10 degrees is added to the offset.

.. code-block:: python

    files = '276-279'

    # Generate file locations as normal
    dataFiles = _tools.fileListGenerator(files,'/path/to/data/')
    
    ds = DataSet.DataSet(dataFiles) # load the data files

    for d in ds: # Loop through all files 
        d.A3Off += 10

    ds.convertDataFile()
    # Continue as normal

Because my crystal has a 120 degrees symmetry, we only scanned 120 degrees in our scans. How do I symmetrize it?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When performing the Q scans most often it is not needed at all to scan 360 degrees as the scattering planes has a rotation symmetry of e.g. 90, 120, or the like. There is no method of the DataSet object to do this but
one can, in a script, load in the measured data twice where the second set of data files is rotated by the symmetry. For the example in this question, I will assume that one has measured the files '276-279' of 120 degrees and the system has 120 degrees symmetry, thus datafiles need to replicated twice and ritated 120 and 240 degrees:

.. code-block:: python

    files = '276-279'

    # Generate file locations as normal
    dataFiles = _tools.fileListGenerator(files,'/path/to/data/')
    numFiles = len(dataFiles) # The number of data files actually measured 
    dataFiles = dataFiles+dataFiles+dataFiles # Repeat the file list to contain 3 copies

    ds = DataSet.DataSet(dataFiles) # load the 12 data files

    for i,d in enumerate(ds): # Loop through all files and keep a counter 
        if i>=numFiles: # If the file number is bigger than 4 add 120 degrees offset to A3
            d.A3Off += 120.0
        if i>=2*numFiles: # If the file number is bigger than 8 add another 120 degrees offset to A3
            d.A3Off += 120.0    

One then has a DataSet object with 12 data files covering the full scattering plane, the so-called donut plot. If a crystal has a symmetry smaller than what is measured (e.g. 110 degrees are measured but symmetry is 90 degrees), there is no problem with this procedure. It merely results in double data coverage at the positions in Q where there is an overlap. 


When I save my nice figures as .eps white lines appear around all pixels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This has to do with how Matplotlib shows figures in contrast to how the eps engine saves the file. To remove the lines add the key word argument 

.. code-block:: python

    linewidth=0,rasterized=True

to the plotting function. This then forces the eps back-end to draw the surrounding edges as the face colour instead of white/transparent.


I cannot change the colour scale when running MJOLNIR through the terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to be able to interact with the figures generated by MJOLNIR when using the terminal make sure to use the interactive flag of matplotlib. That is

.. code-block:: python

    import matplotlib.pyplot as plot
    plt.ion() # Set interactive matplotlib windows

This will make the generation of plots and the plt.show() non-blocking and thus allow the change of axes and other aspects of the plots.

.. _MaskFAQ:

I want to created a mask for my data to exclude specific points in Q
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When dealing with multiple data files at once in a DataSet object, one needs to keep in mind that MJOLNIR does not expect these to be of equal size. That is, there is no requirement for the number of step points to be equal
 (and further, that the instrument is the same for all files). This then results in the mask of at DataSet being a list of np arrays with the same size as the data. Most easily, this is taken care of by looping through the DataFiles in the DataSet as:

.. code-block:: python

    ds = DataSet(dataFiles)
    
    mask = []
    # Loop through all data files in the DataSet
    for d in ds:
        mask.append(d.h>0.0)

    ds.mask = mask

The above code ensures that all points with a value of H larger than 0.0 are masked out. Usually multiple conditions are required for the mask to be correctly created. As an example below is code that creates masks removing all data points within a radius of 0.1 1/A from the provided QPoints.


.. code-block:: python

    QPoints = [[1,1,0],[-1,1,0],[1,0,0],[0,1,0]]
    
    for d in ds:
    # calculate position in qx,qy for QPoints (may differ from file to file)
        localMask = []

        for h,k,l in QPoints:
            qx,qy = d.sample.calculateHKLToQxQy(h,k,l)
            m = np.sqrt((d.qx-qx)**2+(d.qy-qy)**2)<radius
            localMask.append(m)
        trueMask = localMask[0]
        for m in localMask[1:]:
            trueMask = np.logical_or(trueMask,m)
        
        mask.append(trueMask)
    

or in one line

.. code-block:: python

    QPoints = [[1,1,0],[-1,1,0],[1,0,0],[0,1,0]]
    radius = 0.1
    mask = [reduce(np.logical_or,[np.sqrt((d.qx-qx)**2+(d.qy-qy)**2)<radius for qx,qy in [d.sample.calculateHKLToQxQy(*HKL) for HKL in QPoints]]) for d in ds]

    ds.mask = mask



After I have performed a cut, the data I receive as a data type of object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, when working with pandas DataFrames, the data type of the different columns might go from e.g. int64 to object. In order to change this 
when wanting to perform numeric calculations, on can format the DataFrame to be numeric. That is:

.. code-block:: python

    ...
    >>> Data = ds.cut1D(...)
    >>> Data['Intensity'].dtypes
    Name: Intensity, Length: 4244, dtype: object

    >>> Intensity = pd.to_numeric(Data['Intensity'])
    >>> Intensity.dtypes
    dtype('int64')


Is it possible to calculate the resolution of the instrument at a given position in reciprocal space?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Indeed it is through the tools in the Geometry.Instrument settings, but if you want to overplot the resolution ellipsoid directly on your data
this can be done as explain in `Interactivity <InDepthDocumentation/Interactivity.html>`_.


I want to make a cut directly on the 2D data that I have plotted. How can I do this?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Just as for the question above, this is possible as explained in `Interactivity <InDepthDocumentation/Interactivity.html>`_.
