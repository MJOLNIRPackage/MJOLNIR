import sys
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec

import warnings

from matplotlib import rc

from distutils.spawn import find_executable
if find_executable('latex'): # pragma: no cover
    
    rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    USETEX = True
else:
    rc('text', usetex=False)
    USETEX = False
import scipy.optimize
import pyperclip
from MJOLNIR.Statistics.FittingFunction import Gaussian, Lorentz
from MJOLNIR import _tools

class State:#pragma: no cover
    def __init__(self,parent):
        self.parent = parent

    def button(self, event):
        assert 0, "button not implemented"
    def mouse(self, event):
        assert 0, "mouse not implemented"


class Initial(State):#pragma: no cover
    def __init__(self,parent):
        super(Initial,self).__init__(parent)
        self.__name__ = 'Initial'
    def button(self,event):
        if event.key in ['ctrl+i','i']:
            return FitInitialization(self.parent)
        return self
    def mouse(self,event):
        return self
        
            
class FitInitialization(State):#pragma: no cover
    def __init__(self,parent):
        super(FitInitialization,self).__init__(parent)
        
        self.__name__ = 'Fit Initialization - Click to choose bold parameter(s) or number for other models.'
        if USETEX:
            self.__name__ = '\ '.join([x for x in self.__name__.split(' ')])
        self.fitParameterIndex = 0 # Index of current fit parameter
        
        self.ps=' '.join(['{}: {}'.format(i,self.parent.fitObjects[i].__name__) for i in range(len(self.parent.fitObjects))])
        self.parent.updateText(self.fitParameterIndex,title=self.__name__,ps=self.ps)

        
    def button(self,event): # Change fitting function
        if event.key in [str(x) for x in np.arange(len(self.parent.fitObjects))]:
            
            self.parent.fitFunction = self.parent.fitObjects[int(event.key)]()
            self.fitParameterIndex = 0
            self.parent.updateText(self.fitParameterIndex,title=self.__name__,ps=self.ps)
        elif event.key in ['ctrl+i','i']:
            print('Current fit function is: {}({}) with values: {}'.format(self.parent.fitFunction.__name__,\
                  ', '.join([str(x) for x in self.parent.fitFunction.parameterNames]),', '.join([str(x) for x in self.parent.fitFunction.parameters])))
        elif event.key in ['r']:
            self.fitParameterIndex = np.mod(self.fitParameterIndex-1,self.parent.fitFunction.parameterLength)
            self.parent.fitFunction.parameters[self.fitParameterIndex] = np.NaN
            
        elif event.key in ['e']: # Execute fitting
            return Execute(self.parent)
            
        return self            
    def mouse(self,event):
        if not event.inaxes is None: # Only if inside the axis
        
            self.fitParameterIndex = self.parent.fitFunction.setParameter(event,self.fitParameterIndex)
            self.parent.plotFit() # Plot if possible
            
            self.parent.updateText(self.fitParameterIndex,title=self.__name__,ps=self.ps)
        return self

                
class Execute(State):#pragma: no cover
    def __init__(self,parent):
        super(Execute,self).__init__(parent)
        self.__name__='Fit Executed - Press "i" for new fit or "ctrl+c" to copy parameters' 
        if USETEX:
            self.__name__ = '\ '.join([x for x in self.__name__.split(' ')])
        ## Perform fit
        f = lambda *args: self.parent.fitFunction.func.__func__(None,*args)
        xdata = self.parent.xData
        ydata = self.parent.yData
        yerr = self.parent.yErr
        guess = self.parent.fitFunction.parameters

        # Set all nan-values to zero
        guess[np.isnan(guess)]=0.0

        yerr[yerr==0]==1 # Set all zero error values to 1
        fit = scipy.optimize.curve_fit(f,xdata,ydata,p0=guess,sigma=yerr,absolute_sigma=True)
        self.parent.fitFunction.parameters = fit[0]
        self.parent.fitFunction.fitError = fit[1]
        self.parent.plotFit()
        self.parent.updateText(title=self.__name__)
        
    def button(self,event):
        if event.key in ['ctrl+i','i']: # New fitting
            return FitInitialization(self.parent)
        else:
            return self
    
    def mouse(self,event):
        return self
        

    
class Viewer1D: # pragma: no cover
    initialText = 'Press "i" to initialize fitting procedure.'    
    fitObjects = [Gaussian,Lorentz]

    @_tools.KwargChecker
    def __init__(self, XData,YData,YErr,fitFunction=fitObjects[0](),xLabel='',dataLabel='',xID = 0, yID = 0, plotAll = False,**kwargs):
        """Interactive visualization of 1D data with fitting capabilities. Currently only inteded for 1 scan file.
        
        Args:
            
            - XData (list): List of x-valies in shape (n) or (n,m) for n scan points and m scan parameters.
            
            - YData (list): List of y-values in shape (n) or (n,m) or (n,m,l) for n scan points, m and l free.
            
            - YErr (list): List of y errors in same shape as YData.
            
        Kwargs:
            
            - fitFunction (FittingFunction): Custumized object to perform fitting (default Gaussian).
            
            - xLabel (list): X label text in shape (m) for m scan parameters (default '', nothing plotted).
            
            - dataLabel (list): Label to be shown in legend in shape (m) or (m,l), m and l free (default '', nothing plotted)
            
            - xID (int): Index of x axis to be used (default 0)
            
            - yID (int): Index of y axis to be plotted first (default 0)
            
            - plotAll (bool): Boolean deciding whether or not to plot all data (default False)
            
        Raises:
            
            - AttributeError
                
        Example:

        >>> from Mjolnir.Data import DataSet,Viewer1D
        >>> file = 'TestData/ManuallyChangedData/A3.nxs'
        >>> ds = DataSet.DataSet(convertedFiles = file)
        >>> data = ds.extractDetectorData(ID=10)
        >>>
        >>> Y = np.concatenate(data)[:,:5] # Only first 5 energies
        >>> Y_err = np.sqrt(Y) # Calculate errors
        >>> X = np.arange(Y.shape[0])
        >>> 
        >>> xlabel = ['Arbitrary [arb]']
        >>> dataLabel = np.array(['Detector 0: pos 0', 'Detector 0: pos 1', 'Detector 0: pos 2','Detector 0: pos 3', 'Detector 0: pos 4'])
        >>> 
        >>> # Initialize the viewer
        >>> Viewer = Viewer1D.Viewer1D(XData=X,YData=Y,\
        >>> YErr=Y_err,xLabel=xlabel,dataLabel = dataLabel,plotAll=True)

            
        For a walkthrough of the interface see :ref:`Raw plotting and fitting<Raw-plotting-and-fitting>`. 
        """
        ## Make x data into shape (N,M) for N: num of variables, M: scanpoints
        if len(XData.shape) >= 3 or len(XData.shape)==0:
            raise AttributeError('Expected size of xData is 1 or 2 dimensional')
        self.xDataAll = XData
        if len(XData.shape)==1:
            self.xDataAll = np.array(self.xDataAll).reshape(-1,1)
        
        
        if len(YData.shape) >= 4 or len(YData.shape)==0:
            raise AttributeError('Expected size of xData is 1, 2, or 3 dimensional')
        self.yDataAll = YData
        self.yErrAll = YErr
        if len(YData.shape)==1:
            self.yDataAll = np.array(self.yDataAll).reshape(-1,1)
            self.yErrAll = np.array(self.yErrAll).reshape(-1,1)
        elif len(YData.shape)==3: # yDataAll into shape (N,M) as above

            self.yDataAll = np.array(self.yDataAll).reshape(YData.shape[0],-1)
            self.yErrAll = np.array(self.yErrAll).reshape(YData.shape[0],-1)
            
            if not dataLabel is '':
                dataLabel = np.array(dataLabel).reshape(-1)
                dataLabel.dtype = object
                    
        
        if len(xLabel)!=self.xDataAll.shape[1] and xLabel!='':
            raise AttributeError('Provided x-labels do not match number of x-values. ({} labels and {} x-values)'.format(len(xLabel),len(self.xDataAll)))
        if len(dataLabel)!=self.yDataAll.shape[1] and dataLabel!='':
            raise AttributeError('Provided data labels do not match number of y values. ({} labels and {} x-values)'.format(len(dataLabel),len(self.yDataAll)))
        
        
        if xID>self.xDataAll.shape[1]-1:
            raise AttributeError('Provided xId ({}) is outside of x values provided ({})'.format(xID,self.xDataAll.shape[1]))
        self.xID = xID
        if yID>self.yDataAll.shape[1]-1:
            raise AttributeError('Provided yId ({}) is outside of y values provided ({})'.format(yID,self.yDataAll.shape[1]))
        self.yID = yID
        
        if dataLabel == '': # If no label is given, generate them as [0,1,2,3,...]
            dataLabel = [str(x) for x in np.arange(len(self.yDataAll))]
        
        self.xLabel = xLabel
        self.dataLabel = dataLabel
        gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 6]) 
        self.figure = plt.figure()
        self.textAx = plt.subplot(gs[0])
        self.textAx.remove()
        self.ax = plt.subplot(gs[1])
        
        if self.yDataAll.shape[1]>10 and plotAll==True:
            warnings.warn('More than 10 plots are being plotted. This may take some time...',category=RuntimeWarning,stacklevel=2)
        self.plotAll = plotAll
        
        self.text = self.initialText
        
        
        self.figure.canvas.mpl_connect('key_press_event',lambda event: self.button(event) )
        self.figure.canvas.mpl_connect('button_press_event',self.mouse )
        
        self.initData()
        
        
        self.fitFunction = fitFunction
        self.currentState = Initial(self)
        self.fitPlot = None # No initial plot
        
    @property
    def text(self):
        return self._text

    @text.getter
    def text(self):
        return self._text

    
    @text.setter
    def text(self,text):
        self._text = text
        try:
            self.textObject.remove()
        except:
            pass
        BBox = self.textAx.get_position()
        p0 = BBox.p0
        p1 = BBox.p1
        x = p0[0]-0.05
        y = p1[1]+0.1
        self.textObject = self.figure.text(x,y, self._text, fontsize=10, transform=self.figure.transFigure,horizontalalignment='left',\
                                           verticalalignment='top')
        
    def plotData(self):
        """Plot current data. First destroy previous plot if possible"""
        try:
            self.dataPlot.remove()
        except:
            pass
        self.ax.clear()
        self.ax.set_ylabel('Int [arb]')
        if self.plotAll==False:
            self.dataPlot = self.ax.errorbar(self.xData,self.yData,yerr=self.yErr,fmt='.')
            
            
            if not self.dataLabel is '':
                self.ax.legend([self.dataLabel[self.yID]])
                self.ax.legend_.draggable(True)
        else:
            for i in range(self.yDataAll.shape[1]):
                self.dataPlot = self.ax.errorbar(self.xData,self.yDataAll[:,i],yerr=self.yErrAll[:,i],fmt='.')

            if not self.dataLabel is '':
                labels = self.dataLabel.copy()
                    
                labels[self.yID] = '*'+labels[self.yID]
                self.ax.legend(labels)
                self.ax.legend_.draggable(True)
                
        if not self.xLabel is '':
            self.ax.set_xlabel(self.xLabel[self.xID])
        
    def plotFit(self):
        """Plot current guess or fit"""
        if self.fitFunction.executable:
            self.removeFitPlot()
            self.y = self.fitFunction(self.x)
            self.fitPlot = self.ax.plot(self.x,self.y)

    def initData(self):
        """Update with new indices both X and Y (+Yerr)"""
        self.xData = self.xDataAll[:,self.xID]
        self.yData = self.yDataAll[:,self.yID]
        self.yErr = self.yErrAll[:,self.yID]
        self.x = np.linspace(np.min(self.xData),np.max(self.xData),201)
        
        
        self.plotData()
            
    def removeFitPlot(self):
        """Try to remove previous fitPlot if it exists"""
        if not self.fitPlot is None:
            try:
                self.fitPlot[0].remove()
            except:
                pass
                
    def button(self,event):
        if event.key in ['up','down','left','right']: # Change xID or yID
            if event.key == 'left':
                if self.xDataAll.shape[1]>1:
                    self.xID = np.mod(self.xID-1,self.xDataAll.shape[1])
                else:
                    return
            elif event.key == 'right':
                if self.xDataAll.shape[1]>1:
                    self.xID = np.mod(self.xID+1,self.xDataAll.shape[1])
                else:
                    return
            elif event.key == 'down':
                self.yID = np.mod(self.yID-1,self.yDataAll.shape[1])
            elif event.key == 'up':
                self.yID = np.mod(self.yID+1,self.yDataAll.shape[1])
            self.initData()

        elif event.key in ['ctrl+c']:
            pyperclip.copy(', '.join([str(x) for x in self.fitFunction.parameters]))
        else:
            self.currentState = self.currentState.button(event)
    def mouse(self,event):
        self.currentState = self.currentState.mouse(event)    
        
    def updateText(self,highlight=None,title=None,ps=None):
        if USETEX:
            text = '$\mathbf{'+'{}'.format(title)+'}$'+'\nCurrent fitting model: {}\n'.format(self.fitFunction.latex(highlight=highlight))
        else:
            text = '{}'.format(title)+'}'+'\nCurrent fitting model: {}\n'.format(self.fitFunction.latex(highlight=highlight))
        for i in range(self.fitFunction.parameterLength):
            if np.isnan(self.fitFunction.parameters[i]):
                val = '   '
            else:
                val = '{:.3e}'.format(self.fitFunction.parameters[i])
            if USETEX:
                text+='${}$:    {}     '.format(self.fitFunction.variableNames[i],val)
            else:
                text+='{}:    {}     '.format(self.fitFunction.variableNames[i],val)
        if not self.fitFunction.fitError is False:
            text+='\n'
            for i in range(self.fitFunction.parameterLength):
                if np.isnan(self.fitFunction.parameters[i]):
                    val = '   '
                else:
                    val = '{:.3e}'.format(np.sqrt(self.fitFunction.fitError[i,i]))
                if USETEX:
                    text+='${}{}{}$:  {}    '.format('\sigma_{',self.fitFunction.variableNames[i],'}',val)
                else:
                    text+='{}{}{}:  {}    '.format('sigma_',self.fitFunction.variableNames[i],'',val)
        if not ps is None:
            text+='\n'+ps
        self.text = text