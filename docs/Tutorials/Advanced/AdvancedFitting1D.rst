Fitting of sequental 2D plots in as 1D Gaussians
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After a spetrum has been measured, it is sensible to extract some of the physical quantities of excitations. In this tutorial, one possible method of doing so by the use of QELine as well as LMFIT is presented. The data is the MnF2 data set where only as subset of the data is used. Secifically from 4.8 meV to 6.9 meV and the parameters wanted are the positions and widths of the magnon along two cuts in Q. More specifically, there will be performed a cut from (-1,0,-67) RLU to (-1,0,1) RLU to (0,0,1) RLU. For each of the two segments, the data is split into horizontal (energy) slices of same energy width and Gaussians are used to model the spin wave line shapes. Below is the data as seen before the fitting routine has been performed.
 .. figure:: RawQEPlot.png
  :width: 80%
  :align: center

From this it is clear that the number of horizontal Gaussians changes for the second segment. This is all done in the code. Further, the start guess for the parameters are found from the estimators build in to the Gaussian objects of the LMFIT module.

.. code-block:: python
   :linenos:

   
   from MJOLNIR.Data import DataSet
   from MJOLNIR import _tools
   import numpy as np
   from lmfit.models import GaussianModel
   
   numbers = '499,500'
   files = _tools.fileListGenerator(numbers,'/Path/To/Data/',year=2018)
   ds = DataSet.DataSet(files)
   ds.convertDataFile(binning=8,saveFile=False)
   
   # utility function to return points in array before a given value
   def index_of(arrval, value):
   """Return index of array *at or below* value."""
   if value < min(arrval):
       return 0
   if np.sum(np.diff(arrval))>0: # Positive change
       return max(np.where(arrval <= value)[0])
   else:
       return max(np.where(arrval >= value)[0])
   
   
   # Define qut in Q space
   q1 = [-1,0,-1.2]
   q2 = [-1,0,1.0]
   q3 = [0,0,1.0]
   
   # Create two different energy ranges for the cut from q1 to q2 and q2 to q3
   EBins = [np.linspace(ds.energy.min(),ds.energy.max(),31),np.linspace(ds.energy.min(),ds.energy.max(),25)]
   
   # Perform the cut and plot it. Returns a cutObject
   cutObject = ds.plotCutQELine(QPoints=[q1,q2,q3],width=0.05,minPixel=0.01,EnergyBins=EBins)
   
   # Extract axis and save figure
   ax = cutObject._ax
   ax.get_figure('figure0.png',format='png').savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Advanced/RawQEPlot.png',format='png',dpi=600)
   
   # Define plotting parameters for error bar plot to be created
   ErrorBarKwargs = {'markersize':4, 'capsize':2,'elinewidth':1,'markeredgewidth':2,'mfc':'white','fmt':'o'}
   
   out = [] # Holder for fitting results
   
   # Utilize the fitting itterator of the cutObject to loop through data
   for ID,x,data,xplot in cutObject.fittingIterator():
   # ID: which segment, either 0: q1->q2, or 1: q2->q3
   # x: 4D position of center of bins (h,k,l,E)
   # data: Data list for current 1D cutout in the form [Intensity, Normalization, Monitor, Count]
   
   # Transpose as to easier perform fit
   x = x.T
   
   # Calculate intensity and remove points that has not been measured
   y = data[0]*data[-1]/(data[1]*data[2])
   y = y.flatten()
   NoNNans = (data[2]>0).flatten()
   
   y = y[NoNNans]
   x = x[:,NoNNans]
   
   
   # Fit depends on which cut is in question    
   if ID==0:
       # Fit consists of three 1D Gaussians
       gauss1 = GaussianModel(prefix='g1_')
       gauss2 = GaussianModel(prefix='g2_') # Give the Gaussians different prefixes
       gauss3 = GaussianModel(prefix='g3_')
       
       # Model is simply the sum of these
       mod = gauss1 + gauss2 + gauss3
       
       # Direction along which to fit
       fittingX = x[2]
       
       # Cut out portions to make use of the automatic parameter estimation
       ix1 = index_of(fittingX, -0.1) # All points befor -0.1 is for first Gaussian
       ix2 = index_of(fittingX, 0.5) # Between ix1 and ix2 is second Gaussian
       # Third Gaussian is from 0.5 and the rest.
       
       # Use lmfit's parameter starting guess
       pars = gauss1.guess(y[:ix1],x=fittingX[:ix1])
       pars.update(gauss2.guess(y[ix1:ix2],x=fittingX[ix1:ix2]))
       pars.update(gauss3.guess(y[ix2:],x=fittingX[ix2:]))
   
   elif ID == 1:
       # Repeat procedure from above for second segment with only 2 Gaussians or 1 depending on energy
       
       fittingX = x[0]
       gauss1 = GaussianModel(prefix='g1_')
       E = x[-1][0]
       if E<6.1:
   
           gauss2 = GaussianModel(prefix='g2_')
           
           mod = gauss1 + gauss2
           
           ix1 = index_of(fittingX, -0.5)
   
           pars = gauss1.guess(y[:ix1],x=fittingX[:ix1])
           pars.update(gauss2.guess(y[ix1:],x=fittingX[ix1:]))
       elif E<6.6:
           mod = gauss1
           pars = gauss1.guess(y,x=fittingX)
       else:
           continue
   
   # Perform fit and save results in the 'out' array
   result = mod.fit(y, pars, x=fittingX)
   out.append(result)
   
   # For plotting Centres and widths are needed (Their errors are extracted as well, but unused)
   centres = []
   centres_err = []
   widths = []
   widths_err = []
   for parname, par in result.params.items():
       if 'center' in parname:
           centres.append(par.value)
           centres_err.append(par.stderr)
       if 'sigma' in parname:
           widths.append(par.value)
           widths_err.append(par.stderr)
   
   # Depending on segment, the centres of the Gaussians are eihter [-1,0,c] or [c,0,-1]
   # 'Errors' referes to the errorbars and are in this case plotted as the widths
   if ID == 0:
       Center3D    = [[-1.0,0,c] for c in centres]
       Center3DErr = [[-1.0,0,c+err] for c,err in zip(centres,widths)]
   
   elif ID == 1:
       Center3D    = [[c,0,-1] for c in centres]
       Center3DErr = [[c+err,0,-1] for c,err in zip(centres,widths)]
   
   # Calculate the position along the plot for the HKL points.
   XPosition = [ax.converterFunction(c,ID=ID) for c in Center3D]
   XPositionErr = [ax.converterFunction(CE,ID=ID)-XC for CE,XC in zip(Center3DErr,XPosition)]
   # The above is needed as an axis only has one x-axis and plotting multiple segments require some trickery, 
   # resulting in rescaling and offsets. It is all taken care of in the converterFunction of the axis.
   
   
   ax.errorbar(XPosition,[x[-1][0]]*len(XPosition),xerr=XPositionErr,c = 'b', **ErrorBarKwargs)
   # plot the errorbar on top of intensity data
   
   ax.get_figure('figure1.png',format='png').savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Advanced/FittedQEPlot.png',format='png',dpi=600)
   

This is just a simple example code generating a figure denoted "figure0.png". However, in the actual code this figure is called "MJOLNIR/docs/Tutorials/Tools/Here.png" and can be shown with the use of rst commands
 .. figure:: FittedQEPlot.png
  :width: 80%
  :align: center
