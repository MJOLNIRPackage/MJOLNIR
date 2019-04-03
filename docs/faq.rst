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


When I save my nice figures as .eps white lines appear around all pixels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This has to do with how Matplotlib shows figures in contrast to how the eps engine saves the file. To remove the lines add the key word argument 

.. code-block:: python

    edgecolors='face'

to the plotting function. This then forces the eps backend to draw the surrounding edges as the face color instead of white/transparent.
