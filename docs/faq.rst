FAQ
===

A collection of questions and answers ment to guide the user(s) of the MJOLNIR package through some of the most common pitfalls.


My RLU-axis has its tick marks positions stupidly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This 'feature' of the RLU axis is unfortunately linked to the current stage of the development of the `GridHelperCurveLinear <https://matplotlib.org/api/_as_gen/mpl_toolkits.axisartist.grid_helper_curvelinear.GridHelperCurveLinear.html#mpl_toolkits.axisartist.grid_helper_curvelinear.GridHelperCurveLinear>`_.
Because of the current lack of customization for this in-progress work some attributes cannot be changed. However, when generating an new RLUaxis 
one can specify the number of tick marks on the x and y axis. This number stays constant  throughout the lifetime of that axis, but a new can be created with a 
different number if needed. Further, specifying only the number of y-axis tick marks the program calculates a suitable number along the x-axis.


Everytime I run my script, I am bombarded with 'The file xx exists alread. Old file will be renamed to xx.'-warnings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The warning originates from the fact that a data file with the same name has been converted priviously in the sampe directory, 
and the program then tells you that the old conversion will be moved to a new file with the new name. It only serves as to try 
to prevent you from saving a conversion on top of another involuntarily. To  silence this annoyance, convert the data files only 
once or use the 'saveFile=False' flag in the DataSet.convertDataFile() method to not save the converted data.


I do not want to make use of the features in MJOLNIR but only convert and extract my data. How do I do this?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is definitely possible to convert and extract the data from MJOLNIR without many trouples. Below is a code example:

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

The reason for this can be that the normalziation file has not converged for all prismatic pixels. Due to a spurion, the lowest energy pixels across all 
detector tubes are wrongly fitted resulting in them being moved and normalized with a wrong value. This is, however only a temporary problem and does also 
only apply to data sets that are converted with binning = 8. To counter this, apply a mask to the data masking out 
the bottom 2 prismatic pixels as:

.. code-block:: python
   :linenos:

    from MJOLNIR.Data import DataSet
    fileName = ['/Path/To/Data/camea2018n000136.hdf']
    ds = DataSet.DataSet(dataFiles=fileName)
    ds.convertDataFile(binning=8,saveFile=False)

    # Create a masks that has the same shape like e.g. the intensity
    mask = np.zeros_like(DS.I.data,dtype=bool) # Make it boolean as well (not explicitly necessary)
    mask[:,:,:2] = True # Remember the shape of I being (#ScanSteps,#Detectors,#Binning*#Analysers)

    # The #Detectors is 104 and there are #Analysers = 8 analysers

    # Apply the mask to the DataSet object.
    ds.mask = mask


Because my crystal has a 120 degrees symmetry, we only scanned 120 degrees in our scans. How do I symmetrize it?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When performing the Q scans most often it is not needed at all to scan 360 degress as the scattering planes has a rotation symmetry of e.g. 90, 120, or the like. There is no method of the DataSet object to do this but
one can, in a script, load in the measured data twice where the second set of data files is rotated by the symmetry. For the example in this question, I will assume that one has measured the files '276-279' of 120 degrees and the system has 120 degrees symmetry, thus datafiles need to replicated twice and ritated 120 and 240 degrees:

.. code-block:: python

	files = '276-279'

	# Generate file locations as normal
	dataFiles = _tools.fileListGenerator(files,'/paht/to/data/')
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

    edgecolors='face'

to the plotting function. This then forces the eps back-end to draw the surrounding edges as the face colour instead of white/transparent.


I cannot change the colour scale when running MJOLNIR through the terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to be able to interact with the figures generated by MJOLNIR when using the terminal make sure to use the interactive flag of matplotlib. That is

.. code-block:: python

    import matplotlib.pyplot as plot
    plt.ion() # Set interactive matplotlib windows

This will make the generation of plots and the plt.show() non-blocking and thus allow the change of axes and other aspects of the plots.


I want to created a mask for my data to exclude specific points in Q
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When dealing with multiple data files at once in a DataSet object, one needs to keep in mind that MJOLNIR does not expect these to be of equal size. That is, there is no requirement for the number of step points to be equal
 (and further, that the instrument is the same for all files). This then results in the mask of at DataSet being a list of np arrays with the same size as the data. Most easily, this is takne care of by looping through the DataFiles in the DataSet as:

.. code-block:: python

	ds = DataSet(dataFiles)
	
	mask = []
	# Loop through all data files in the DataSet
	for d in ds:
		mask.append(d.h>0.0)

	ds.mask = mask

The above code ensures that all points with a value of H larger than 0.0 are masked out. Usually multiple conditions are required for the mask to be corretly created. As an example below is code that creates masks removing all data points within a radius of 0.1 1/A from the provided QPoints.


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


