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
from MJOLNIR.Statistics.FittingFunction import *
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
        
        self.__name__ = r'Fit Initialization - Click to choose "\{"parameter(s)"\}" or number for other models.'
        #if USETEX:
        #    self.__name__ = '\ '.join([x for x in self.__name__.split(' ')])
        self.fitParameterIndex = 0 # Index of current fit parameter
        
        self.ps=' '.join(['{}: {}'.format(i,self.parent.fitObjects[i].__name__) for i in range(len(self.parent.fitObjects))])
        self.ps+='\nPress e for execution.'
        self.parent.updateText(self.fitParameterIndex,title=self.__name__,ps=self.ps)
        plt.draw()
        
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
        plt.draw()    
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
        self.__name__=r'Fit Executed - Press "i" for new fit or "ctrl+c" to copy parameters' 
        #if USETEX:
        #    self.__name__ = '\ '.join([x for x in self.__name__.split(' ')])
        ## Perform fit
        f = lambda *args: self.parent.fitFunction.func.__func__(None,*args)
        xdata = self.parent._xData
        ydata = self.parent._yData
        yerr = self.parent._yErr
        guess = self.parent.fitFunction.parameters

        # Set all nan-values to zero
        guess[np.isnan(guess)]=0.0

        #yerr[ydata==0]=1 # Set all zero error values to 1
        fit = scipy.optimize.curve_fit(f,xdata,ydata,p0=guess,sigma=yerr,absolute_sigma=True)
        self.parent.fitFunction.parameters = fit[0]
        self.parent.fitFunction.fitError = fit[1]
        self.parent.plotFit()
        self.parent.updateText(title=self.__name__)
        plt.draw()
        
    def button(self,event):
        if event.key in ['ctrl+i','i']: # New fitting
            return FitInitialization(self.parent)
        else:
            return self
    
    def mouse(self,event):
        return self
        

    
class Viewer1D: # pragma: no cover
    initialText = r'Press "i" to initialize fitting procedure.'    
    fitObjects = FittingFunction.__subclasses__()

    @_tools.KwargChecker()
    def __init__(self, XData,YData,YErr,fitFunction=fitObjects[0](),xLabel='',dataLabel='',xID = 0, plotAll = False,**kwargs):
        """Interactive visualization of 1D data with fitting capabilities. Currently only inteded for 1 scan file.
        
        Args:
            
            - XData (list): List of x-valies in shape (m,n) for m data series and n scan points.
            
            - YData (list): List of y-values in shape (m,n) for m data series and n scan points.
            
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
                
        Example: # TODO: REDO!!

        >>> from MJOLNIR.Data import DataSet,Viewer1D
        >>> file = 'TestData/ManuallyChangedData/A3.h5'
        >>> ds = DataSet.DataSet(dataFiles = file)
        >>> ds.convertDataFile(binning=1)
        >>> data = ds.extractData(A4Id=30)
        >>>
        >>> Y = data[:,:5] # Only first 5 energies
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
            raise AttributeError('Expected size of xData is 2 dimensions')
        
        if np.all([len(YData[i]) != len(XData[i]) for i in range(len(YData))]):
        #if len(YData.shape) >= 4 or len(YData.shape)==0:
            raise AttributeError('Shape of YData({}) does not match XData({}).'.format(YData.shape,XData.shape))
        if YErr.shape != YData.shape:
            raise AttributeError('Shape of YErr{}) does not match YData({}).'.format(YErr.shape,YData.shape))
                
        if not dataLabel is '':
            dataLabel = np.array(dataLabel).reshape(-1)
            #dataLabel.dtype = object
        else:
            dataLabel = [str(x) for x in np.arange(len(YData))]    
        
        if len(xLabel)!=len(XData) and xLabel!='':
            raise AttributeError('Provided x-labels do not match number of x-values. ({} labels and {} x-values)'.format(len(xLabel),len(self.XData)))
        if len(dataLabel)!=len(YData) and dataLabel!='':
            raise AttributeError('Provided data labels do not match number of y values. ({} labels and {} x-values)'.format(len(dataLabel),len(YData)))
        
        
        if xID>len(XData)-1:
            raise AttributeError('Provided xId ({}) is outside of x values provided ({})'.format(xID,len(XData)))
        
        
        # Save all data to self
        self.xData = XData
        self.yData = YData
        self.yErr = YErr
        

        self.xLabel = xLabel
        self.dataLabel = dataLabel
        self.xID = xID
        
        gs = matplotlib.gridspec.GridSpec(2, 1, height_ratios=[1, 6]) 
        self.figure = plt.figure()
        self.textAx = plt.subplot(gs[0])
        self.textAx.remove()
        self.ax = plt.subplot(gs[1])
        
        if len(self.yData)>10 and plotAll==True:
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
            self.dataPlot = self.ax.errorbar(self._xData,self._yData,yerr=self._yErr,fmt='.')
            
            
            if not self.dataLabel is '':
                self.ax.legend([self.dataLabel[self.xID]])
                self.ax.legend_.draggable(True)
        else:
            for i in range(len(self.yData)):
                self.dataPlot = self.ax.errorbar(self.xData[i],self.yData[i],yerr=self.yErr[i],fmt='.')

            if not self.dataLabel is '':
                labels = self.dataLabel.copy()
                    
                labels[self.xID] = '*'+labels[self.xID]
                self.ax.legend(labels)
                self.ax.legend_.draggable(True)
                
        if not self.xLabel is '':
            self.ax.set_xlabel(self.xLabel[self.xID])
        plt.draw()
        
    def plotFit(self):
        """Plot current guess or fit"""
        if self.fitFunction.executable:
            self.removeFitPlot()
            self.y = self.fitFunction(self.x)
            self.fitPlot = self.ax.plot(self.x,self.y)
            plt.draw()

    def initData(self):
        """Update with new indices both X and Y (+Yerr)"""
        self._xData = self.xData[self.xID]
        self._yData = self.yData[self.xID]
        self._yErr = self.yErr[self.xID]
        self.x = np.linspace(np.min(self._xData),np.max(self._xData),len(self._xData)*10)
        self.plotData()
            
    def removeFitPlot(self):
        """Try to remove previous fitPlot if it exists"""
        if not self.fitPlot is None:
            try:
                self.fitPlot[0].remove()
                plt.draw()
            except:
                pass
                
    def button(self,event):
        if event.key in ['up','down','left','right']: # Change xID or yID
            if event.key == 'left':
                if len(self.xData)>1:
                    self.xID = np.mod(self.xID-1,len(self.xData))
                else:
                    return
            elif event.key == 'right':
                if len(self.xData)>1:
                    self.xID = np.mod(self.xID+1,len(self.xData))
                else:
                    return
            #elif event.key == 'down':
            #    self.yID = np.mod(self.yID-1,self.yDataAll.shape[1])
            #elif event.key == 'up':
            #    self.yID = np.mod(self.yID+1,self.yDataAll.shape[1])
            self.initData()

        elif event.key in ['ctrl+c']:
            pyperclip.copy(', '.join([str(x) for x in self.fitFunction.parameters]))
        else:
            self.currentState = self.currentState.button(event)
    def mouse(self,event):
        self.currentState = self.currentState.mouse(event)    
        
    def updateText(self,highlight=None,title=None,ps=None):
        #if USETEX:
        text = '{}'.format(title)+'\nCurrent fitting model: {}\n'.format(self.fitFunction.latex(highlight=highlight))
        #else:
        #    text = '{}'.format(title)+'}'+'\nCurrent fitting model: {}\n'.format(self.fitFunction.latex(highlight=highlight))
        for i in range(self.fitFunction.parameterLength):
            if np.isnan(self.fitFunction.parameters[i]):
                val = r'   '
            else:
                val = r'{:.3e}'.format(self.fitFunction.parameters[i])
            #if USETEX:
            text+=r'${}$:    {}     '.format(self.fitFunction.variableNames[i],val)
            #else:
            #    text+=r'{}:    {}     '.format(self.fitFunction.variableNames[i],val)
        if not self.fitFunction.fitError is False:
            text+='\n'
            for i in range(self.fitFunction.parameterLength):
                if np.isnan(self.fitFunction.parameters[i]):
                    val = r'   '
                else:
                    val = r'{:.3e}'.format(np.sqrt(self.fitFunction.fitError[i,i]))
                #if USETEX:
                text+=r'${}{}{}$:  {}    '.format('\sigma_{',self.fitFunction.variableNames[i],'}',val)
                #else:
                #    text+=r'{}{}{}:  {}    '.format('sigma_',self.fitFunction.variableNames[i],'',val)
        if not ps is None:
            text+='\n'+ps
        self.text = text