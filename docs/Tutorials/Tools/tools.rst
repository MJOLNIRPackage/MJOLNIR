Overview of _tools
^^^^^^^^^^^^^^^^^^
Some different tools have been developed taking care of different aspects of the MJOLNIR package. Some are used as decorators and does thus not make sense to utilize outside of the package itself. Below, the binEdges and fileListGenerator will be shown. 

.. code-block:: python
   :linenos:

   import numpy as np
   import matplotlib.pyplot as plt
   from MJOLNIR import _tools
   from MJOLNIR.Data import DataSet
   # To generate a list of data files from nubers
   
   numbers = '457-458,460,461'
   files = _tools.fileListGenerator(numbers,'/Path/To/Data/',year=2018)
   
   ## Creating a non-linear distribution of points
   points = np.exp(np.linspace(-1,np.log(10),21))
   
   minimumBinSize = 1
   
   Bins = _tools.binEdges(points,minimumBinSize)
   
   fig,ax = plt.subplots()
   
   ax.scatter(points,np.ones_like(points),c='r',label='Data Points',zorder=2)
   [ax.plot([Bins[i],Bins[i]],[0.5,1.5],'k',zorder=1) for i in range(len(Bins))]
   ax.plot(np.concatenate([Bins,np.flip(Bins)]),np.concatenate([np.ones_like(Bins)*0.5,np.ones_like(Bins)*1.5]),c='k',label='Bins',zorder=1)    
   ax.set_xticks(np.linspace(0,10,11))
   ax.set_yticks([])
   ax.set_ylim(-1,3)
   ax.grid(True,c='k',zorder=0)
   fig.legend()
   fig.savefig('figure0.png',format='png')
   

Generate Data Path
------------------
Mostely, one needs to combine files with different run numbers. These are often sequential and thus making a loop to generate the file path is straight forward. However, sometimes one is to combine files differently. This is taken care of with the function _tools.fileListGenerator, which takes a number string, a folder, and the year of the data files. From this data files are generated. The number string works by creating paths to all files between and including numbers provided like "xx-yy". Adding numbers by comma separation simply appends it. Further, a combination of these two methods is possible. The resulting files from above is: 

 - /home/lass/Dropbox/PhD/CameaData/camea2018n000457.hdf
 - /home/lass/Dropbox/PhD/CameaData/camea2018n000458.hdf
 - /home/lass/Dropbox/PhD/CameaData/camea2018n000460.hdf
 - /home/lass/Dropbox/PhD/CameaData/camea2018n000461.hdf


Binning
-------
When dealing with data that is non-uniformly distributed it is sometimes necessary to have a function that creates suitable bins all containing  data. This is exactly what the _tools.binEdges does. It starts the first bin starts 0.1 times the tolerance away from the first unique entry and itterates through the list of points. For each step, as long as the distance between the current value and the next is less than the tolerance the bin is expanded. If the distance is greater, a bin is created in the middle between the last accepted and the rejected point. This insures that all points are within a bin and that no bins are empty. For the last point of the list, it is checked if the difference between last bin edge and point is smaller than tolerance*1.1 and if so, a bin is created of size tolerance. Otherwise, the las bin edge is 0.1*tolerance away from last point. An example is shown in the above code generating the figure below:
 .. figure:: Binning.png
  :width: 30%
  :align: center

