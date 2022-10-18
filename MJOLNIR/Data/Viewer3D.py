import os
import sys
from typing import Type

sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')

import warnings
import copy
import matplotlib.cm as cm
import matplotlib.gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
from MJOLNIR import _tools
from MJOLNIR._interactiveSettings import Viewer3DSettings, States, cut1DHolder
import functools

pythonVersion = sys.version_info[0]
pythonSubVersion = sys.version_info[1]

class Viewer3D(object):  
    @_tools.KwargChecker(include=[_tools.MPLKwargs])
    def __init__(self,Data,bins,axis=2, log=False ,ax = None, grid = False, adjustable=True, outputFunction=print, 
                 cmap=None, CurratAxeBraggList=None, plotCurratAxe=False,Ei=None,EfLimits=None, dataset = None, cut1DFunctionRectangle=None,\
                    cut1DFunctionCircle = None, cut1DFunctionRectanglePerp=None,cut1DFunctionRectangleHorizontal=None,cut1DFunctionRectangleVertical=None, **kwargs):#pragma: no cover
        """3 dimensional viewing object generating interactive Matplotlib figure. 
        Keeps track of all the different plotting functions and variables in order to allow the user to change between different slicing modes and to scroll through the data in an interactive way.

        Args:

            - Data (3D array): Intensity array in three dimensions. Assumed to have Qx, Qy, and E along the first, second, and third directions respectively.

            - bins (List of 1D arrays): Coordinates of the three directions as returned by the BinData3D functionality of DataSet.

        Kwargs:

            - axis (int): Axis along which the interactive plot slices the data (default 2).

            - log (bool): If true, the log 10 of the intensity is plotted (default False).
            - ax (matplotlib axis): Matplotlib axis into which one plots data (Default None).

            - grid (bool/int): If int, grid will be plotted with zorder=int, if True, grid is plotted at zorder=-10 (Default False).

            - adjustable (bool): If set true, 2 sliders will be present allowing to fine tune the c-axis (Default True)

            - outputFunction (function): Function called on output string (default print)

            - CurratAxeBraggList (list): List of Bragg peaks for which Currat-Axe spurions are to be calculated, [[H1,K1,L1],[H2,K2,L2],..] (default None)

            - plotCurratAxe (bool): Flag to determine whether or not to plot the Currat-Axe spurions
            
            - Ei (list): List of incoming energies (default None)
            
            - EfLimits (float,float): Minimal and maximal Ef for current instrument (default None)

            - dataset (DataSet): Reference to data set used - needed for interactive cutting with (default None)

            - cut1DFunctionRectangle (function): Function to be called when performing an interactive rectangle (default None)

            - cut1DFunctionCircle (function): Function to be called when performing an interactive circle (default None)

        For an example, see the `quick plotting tutorial <../Tutorials/Quick/QuickView3D.html>`_ under scripting tutorials.

        """
        if cmap is None:
            cmap = 'viridis'
        
        self.Ei = Ei
        self.EfLimits = EfLimits

        self.outputFunction = outputFunction


        self.currentData = None
        self.ds = dataset

        self.plotCurratAxe = plotCurratAxe # Set to false but will change to True when correct list is created
        self._CurratAxeBraggList = None
        if len(Data)==4: # If data is provided as I, norm, mon, normcount
            with warnings.catch_warnings() as w:
                warnings.simplefilter("ignore")
                self.Data = np.divide(Data[0]*Data[3],Data[1]*Data[2])
            
            self.Counts,self.Monitor,self.Normalization,self.NormCounts = Data
            self.allData = True
        elif len(Data) == 2: # From bin3D with no normalization or monitor
            with warnings.catch_warnings() as w:
                warnings.simplefilter("ignore")
                self.Data = np.divide(Data[0],Data[1])
            self.allData = False
        else:
            self.Data = Data
            self.allData = False
        if log:
            self.Data = np.log10(self.Data+1e-20)
        self.bins = bins
        self.dataLimits = [np.nanmin(Data),np.nanmax(Data)]

        gs = matplotlib.gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 
        
        if not grid == False:
            gridArg = grid
            if isinstance(gridArg,(bool)):
                self.grid = True
                self.gridZOrder=-10
            elif isinstance(gridArg,(int,float)):
                self.grid = True
                self.gridZOrder=gridArg
        else:
            self.grid = False
            self.gridZOrder = 0

        if ax is None: # TODO: REDO with actual correct axes!
            self.figure = plt.figure()
            self.ax = plt.subplot(gs[0])#self.figure.add_subplot(111)
            self.ax.isActive = lambda : True
            self.ax.drawState = States.INACTIVE
            self.xlabel = r'Qx [$A^{-1}$]'
            self.ylabel = r'Qy [$A^{-1}$]'
            self.zlabel = r'E [meV]'
            self.rlu = False
            self._axes = [self.ax]
        else:
            if isinstance(ax,plt.Axes): # Assuming only RLU - energy plot is provided
                self.axRLU = ax
                self.figure = ax.get_figure() # Get the correct figure
                self.axNorm,ax2  = self.figure.subplots(1,2,gridspec_kw={'width_ratios':[4, 1]}) # Create figure on top of the other
                ax2.remove() # Remove the excess figure
                
                self.axRLU.set_position(self.axNorm.get_position()) # Update RLU to correct position

                self._axes = [self.axNorm,self.axNorm,self.axRLU]
                self._axes[0].set_xlabel(r'Qx [$A^{-1}$]')
                self._axes[0].set_ylabel(r'E [meV]')
                self._axes[1].set_xlabel(r'Qy [$A^{-1}$]')
                self._axes[1].set_ylabel(r'E [meV]')
                self.ax = self.axNorm
                self.xlabel = r'Qx [$A^{-1}$]'
                self.ylabel = r'Qy [$A^{-1}$]'
                self.zlabel = 'E [meV]'
                self.rlu = True
            elif len(ax)==3: # All axes provided in order QxE,QyE,QxQy
                self.axQxE = ax[0]
                self.axQyE = ax[1]
                self.axRLU = ax[2]
                self.figure = self.axQyE.get_figure() # Get the correct figure
                self.axNorm,ax2  = self.figure.subplots(1,2,gridspec_kw={'width_ratios':[4, 1]}) # Create figure on top of the other
                ax2.remove() # Remove the excess figure
                
                self.axQxE.set_position(self.axNorm.get_position()) # Update axQxE to correct position
                self.axQyE.set_position(self.axNorm.get_position()) # Update axQyE to correct position
                self.axRLU.set_position(self.axNorm.get_position()) # Update axRLU to correct position
                self._axes = [self.axQxE,self.axQyE,self.axRLU]
                
                self.ax = self.axNorm
                hkl = ['H','K','L']
                xlabelSplit = self.axQyE.get_xlabel().replace(' [RLU]','').split(',')
                ylabelSplit = self.axQxE.get_xlabel().replace(' [RLU]','').split(',')
                self.xlabel = '\n'.join(['{}: '.format(hkl[i])+'{:+.3f}'.format(float(xlabelSplit[i])) for i in range(len(xlabelSplit))])
                self.ylabel = '\n'.join(['{}: '.format(hkl[i])+'{:+.3f}'.format(float(ylabelSplit[i])) for i in range(len(ylabelSplit))])
                self.zlabel = 'E [meV]'
                self.rlu = True
                
                self.EnergySliderTransform=[self.axQxE._length,self.axQyE._length,1.0] # Factor to divide the Energy slider value with (only applicable for QE axes)

            else:
                raise AttributeError('Number of provided axes is {} but only 1 or 3 is accepted.'.format(len(ax)))

        self.figure.set_size_inches(11,7)
        self.value = 0
        self.figure.subplots_adjust(bottom=0.25)
        self.cmap = cmap # Update to accommodate deprecation warning
        self.value = 0
        self.setAxis(2)
        # Set up interactive generation of 1DCuts


        viewAxis = axis
        axis_color='white'
        
        self.figure.canvas.mpl_connect('key_press_event',lambda event: onkeypress(event, self) )
        self.figure.canvas.mpl_connect('scroll_event',lambda event: onscroll(event, self))
        
        zeroPoint = np.argmin(np.abs(0.5*(self.Z[0,0][1:]+self.Z[0,0][:-1])))
        
    
        self.Energy_slider_ax = self.figure.add_axes([0.15, 0.1, 0.65, 0.03])
        

        self.Energy_slider = Slider(self.Energy_slider_ax, label=self.label, valmin=self.lowerLim, valmax=self.upperLim, valinit=zeroPoint)
        self.Energy_slider.valtext.set_visible(False)
        
        self.Energy_slider.on_changed(lambda val: sliders_on_changed(self,val))
            
        if self.rlu:
            self.units = ['','',' meV']
        else:
            self.units = [' 1/AA',' 1/AA',' meV']
            
        textposition = [self.Energy_slider_ax.get_position().p1[0]+0.005,self.Energy_slider_ax.get_position().p0[1]+0.005]
        self.text = self.figure.text(textposition[0], textposition[1],s=self.stringValue())
        self.shading = 'flat'
        #self.imcbaxes = self.figure.add_axes([0.0, 0.2, 0.2, 0.7])
        #self.im = self.ax.imshow(self.masked_array[:,:,self.value].T,cmap=self.cmap,extent=[self.X[0],self.X[-1],self.Y[0],self.Y[-1]],origin='lower')
        if self.shading=='flat':
            self.ax.grid(False)
            self.im = self.ax.pcolormesh(self.X[:,:,0].T,self.Y[:,:,0].T,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=cmap)
        elif self.shading=='gouraud':  # pragma: no cover
            XX = 0.5*(self.X[:-1,:-1,self.value]+self.X[1:,1:,self.value]).T
            YY = 0.5*(self.Y[:-1,:-1,self.value]+self.Y[1:,1:,self.value]).T
            self.ax.grid(False)
            self.im = self.ax.pcolormesh(XX,YY,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=cmap) # ,vmin=1e-6,vmax=6e-6
        else:
            raise AttributeError('Did not understand shading {}.'.format(self.shading))
        if self.grid:
            self.ax.grid(self.grid,zorder=self.gridZOrder)
        else:
            self.ax.grid(self.grid)

        self._caxis = self.im.get_clim()
        self.figpos = [0.125,0.25,0.63,0.63]#self.ax.get_position()
        
        self.cbaxes = self.figure.add_axes([0.8, 0.2, 0.03, 0.7])
        self.colorbar = self.figure.colorbar(self.im,cax = self.cbaxes)
        warnings.simplefilter("ignore")
        #self.figure.tight_layout(rect=[0,0.1,0.9,0.9])
        warnings.simplefilter("once")

        self.text.set_text(self.stringValue())
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        #if self.axis == 2:
        #    self.ax.set_xlim(np.min([xlim[0],ylim[0]]),np.max([xlim[1],ylim[1]]))


        dQx = np.diff(self.bins[0][:2,0,0])[0]
        dQy = np.diff(self.bins[1][0,:2,0])[0]
        ddE = np.diff(self.bins[2][0,0,:2])[0]
        self.dQE = [dQx,dQy,ddE]

        self.resolution = [0.05,0.05]
        

        self.Energy_slider.set_val(self.value)
        if self.rlu:

            self.axRLU.onClick = lambda event: eventdecorator(onclick,self,self.axRLU,event,outputFunction=outputFunction)#,extra=' axRLU')
            self.axRLU._button_press_event = self.figure.canvas.mpl_connect('button_press_event', self.axRLU.onClick)

            def updateLimitsAxRLU(self,lower,upper):
                self.EMin = lower
                self.EMax = upper
            self.axRLU.updateLimits = lambda lower,upper: updateLimitsAxRLU(self.axRLU,lower,upper)
            
            self.axQxE.onClick = lambda event: eventdecorator(onclick,self,self.axQxE,event,outputFunction=outputFunction)#,extra=' axQxE')
            self.axQxE._button_press_event = self.figure.canvas.mpl_connect('button_press_event', self.axQxE.onClick)

            def updateLimitsAxQ(self,lower,upper):
                self.QPoints = [self.sample.calculateHKLToQxQy(*x) for x in np.array([self.calculateRLU(-1,0)[0],self.calculateRLU(1,0)[0]])]
                
            self.axQxE.updateLimits = lambda lower,upper: updateLimitsAxQ(self.axQxE,lower,upper)
            self.axQyE.updateLimits = lambda lower,upper: updateLimitsAxQ(self.axQyE,lower,upper)
            self.axQyE.onClick = lambda event: eventdecorator(onclick,self,self.axQyE,event,outputFunction=outputFunction)#,extra=' axQyE')
            self.axQyE._button_press_event = self.figure.canvas.mpl_connect('button_press_event', self.axQyE.onClick)

            self.axRLU.rlu =  self.axQxE.rlu = self.axQyE.rlu = True
            
            # Set up initial active axes
            self.axRLU.isActive = True
            self.axQxE.isActive = False
            self.axQyE.isActive = False

            if not cut1DFunctionRectangle is None:
                self.axRLU.cut1DFunctionRectangle = lambda dr: cut1DFunctionRectangle(viewer=self,dr=dr)

            if not cut1DFunctionCircle is None:
                self.axRLU.cut1DFunctionCircle = lambda dr: cut1DFunctionCircle(viewer=self,dr=dr)

            
            if not cut1DFunctionRectanglePerp is None:
                self.axQxE.cut1DFunctionRectanglePerpendicular = lambda dr: cut1DFunctionRectanglePerp(viewer=self,dr=dr)
                self.axQyE.cut1DFunctionRectanglePerpendicular = lambda dr: cut1DFunctionRectanglePerp(viewer=self,dr=dr)

            if not cut1DFunctionRectangleHorizontal is None:
                self.axQxE.cut1DFunctionRectangleHorizontal = lambda dr: cut1DFunctionRectangleHorizontal(viewer=self,dr=dr)
                self.axQyE.cut1DFunctionRectangleHorizontal = lambda dr: cut1DFunctionRectangleHorizontal(viewer=self,dr=dr)

            if not cut1DFunctionRectangleVertical is None:
                self.axQxE.cut1DFunctionRectangleVertical = lambda dr: cut1DFunctionRectangleVertical(viewer=self,dr=dr)
                self.axQyE.cut1DFunctionRectangleVertical = lambda dr: cut1DFunctionRectangleVertical(viewer=self,dr=dr)

        else:
            self.ax._button_press_event = self.figure.canvas.mpl_connect('key_press_event',lambda event: onkeypress(event, self) )
            self.onClick = lambda event: eventdecorator(onclick,self,self.ax,event,outputFunction=outputFunction)
            self.button_press_event = self.figure.canvas.mpl_connect('button_press_event', self.onClick)

            self.ax.rlu =  False
            
        
        
        try:
            maxVal = np.nanmax(self.masked_array[np.isfinite(self.masked_array)])
        except ValueError:
            maxVal = 1
        self.caxis = [np.nanmin(self.masked_array),maxVal]
        if self.grid:
            self.ax.grid(self.grid,zorder=self.gridZOrder)
        else:
            self.ax.grid(self.grid)

        self.setAxis(viewAxis) # Set view plane to correct
        ## Hack for this to look nice as just changing direction does not render correctly
        self.setPlane(1)
        self.setPlane(0)

        if not self.EfLimits is None:
            self.Ef = np.arange(*self.EfLimits,self.dQE[2])
            
            self.CurratAxeBraggList = CurratAxeBraggList

        if adjustable and pythonVersion>2 and pythonSubVersion>5:

            ax_cmin = plt.axes([0.87, 0.1, 0.05, 0.7])
            ax_cmax  = plt.axes([0.93, 0.1, 0.05, 0.7])

            addColorbarSliders(self,c_min=self.caxis[0],c_max=self.caxis[1],c_minval=self.caxis[0],\
                c_maxval=self.caxis[1],ax_cmin=ax_cmin,ax_cmax=ax_cmax,log=False)

    @property 
    def caxis(self):
        return self._caxis

    @caxis.getter
    def caxis(self):
        return self._caxis

    @caxis.setter
    def caxis(self,caxis):
        ErrMsg = 'Provided caxis is not of correct format. Expected 2 values but recieved "{}" of type {}'
        if not isinstance(caxis,(list,np.ndarray,tuple)):
            raise AttributeError(ErrMsg.format(caxis,type(caxis)))
        if len(list(caxis))!=2:
            raise AttributeError(ErrMsg.format(caxis,type(caxis)))
        self._caxis = caxis
        self.im.set_clim(caxis)
        cmin,cmax = caxis
        fig = self.ax.get_figure()
        #if not _internal:
        if hasattr(fig,'s_cmin'): # The method addColorbarSliders has been called
            for s in [fig.s_cmin,fig.s_cmax]:
                if s._log:
                    s.valmin = np.log10(cmin)
                    s.valmax = np.log10(cmax)
                else:
                    s.valmin = cmin
                    s.valmax = cmax
                    s.ax.set_ylim(s.valmin,s.valmax)
            fig.s_cmin.set_val(cmin)
            fig.s_cmax.set_val(cmax)

        self.colorbar.update_normal(self.im)

    def set_clim(self,cmin,cmax=None):
        if cmax is None:
            self.caxis = cmin
        else:
            self.caxis = (cmin,cmax)

    def setAxis(self,axis):
        if hasattr(self,'im'): # this function is also called before any plot has been performed
            self.im.set_array(self.emptyData) # Set data in current im to a fully masked data set
        if axis==2:
            if self.rlu:
                if hasattr(self.ax,'_button_press_event'):
                    self.ax.get_figure().canvas.mpl_disconnect(self.ax._button_press_event)
                self.ax.isActive = False
                self.figure.delaxes(self.ax)
                self.ax = self.figure.add_axes(self._axes[axis])
                self.ax.isActive = True
                if hasattr(self.ax,'onClick'):
                    self.ax._button_press_event = self.ax.get_figure().canvas.mpl_connect('button_press_event', self.ax.onClick)
            else:
                self.ax.set_xlabel(self.xlabel)
                self.ax.set_ylabel(self.ylabel)
            axes = (0,1,2)
            label = self.zlabel#self.ax.get_ylabel
        elif axis==1:  # pragma: no cover
            if self.rlu:
                if hasattr(self.ax,'_button_press_event'):
                    self.ax.get_figure().canvas.mpl_disconnect(self.ax._button_press_event)
                self.ax.isActive = False
                self.figure.delaxes(self.ax)
                self.ax = self.figure.add_axes(self._axes[axis])
                self.ax.isActive = True
                if hasattr(self.ax,'onClick'):
                    self.ax._button_press_event = self.ax.get_figure().canvas.mpl_connect('button_press_event', self.ax.onClick)
            else:
                self.ax.set_xlabel(self.xlabel)
                self.ax.set_ylabel(self.zlabel)
            axes = (0,2,1)
            label =  self.ylabel#self.ax.get_ylabel
        elif axis==0:  # pragma: no cover
            if self.rlu:
                if hasattr(self.ax,'_button_press_event'):
                    self.ax.get_figure().canvas.mpl_disconnect(self.ax._button_press_event)
                self.ax.isActive = False
                self.figure.delaxes(self.ax)
                self.ax = self.figure.add_axes(self._axes[axis])
                self.ax.isActive = True
                if hasattr(self.ax,'onClick'):
                    self.ax._button_press_event = self.ax.get_figure().canvas.mpl_connect('button_press_event', self.ax.onClick)
            else:
                self.ax.set_xlabel(self.ylabel)
                self.ax.set_ylabel(self.zlabel)
            axes = (1,2,0)
            label =  self.xlabel#self.ax.get_xlabel()
            
        else:
            raise AttributeError('Axis provided not recognized. Should be 0, 1, or 2 but got {}'.format(axis))
        #self.ax.format_coord = self._axes[axis].format_coord
        if hasattr(self.ax,'_step'):
            self.ax._step=self.calculateValue()
        X=self.bins[axes[0]].transpose(axes)
        Y=self.bins[axes[1]].transpose(axes)
        Z=self.bins[axes[2]].transpose(axes)

        
        masked_array = np.ma.array (self.Data, mask=np.isnan(self.Data)).transpose(axes)
        self.emptyData = masked_array[:,:,0].T.flatten().copy()
        self.emptyData.mask = np.ones_like(self.emptyData,dtype=bool)
        self._axesChanged = True
        upperLim = self.Data.shape[axis]-1
        self.label = label
        self.X = X
        self.Y = Y
        self.Z = Z
        self.masked_array = masked_array
        self.axes = axes
        self.upperLim = upperLim
        self.lowerLim = 0
        self.axis = axis

    def calculateValue(self,index=None):
        if index is None:
            index = int(self.value)
        try:
            val = 0.5*(self.Z[0,0,index+1]+self.Z[0,0,index])
        except:
            val = 0.5*(2*self.Z[0,0,index]-self.Z[0,0,index-1])
        if hasattr(self,'EnergySliderTransform'):
            val/=self.EnergySliderTransform[self.axis]
        return val
        
    def stringValue(self):
        unit = self.units[self.axis]
        val = self.calculateValue()
        return str(np.round(val,2))+unit
    
    def setProjection(self,value):
        """Change projection between Qx,Qy, and E, or along principal, orthogonal Q direction, or E if plotting in RLU."""
        self.figure.canvas.key_press_event(str(value))

    def setPlane(self,value):
        """Change plotting plane to new along same axis"""
        self.Energy_slider.set_val(value)
        
    @property
    def CurratAxeBraggList(self):
        return self._CurratAxeBraggList
    
    @CurratAxeBraggList.setter
    def CurratAxeBraggList(self,peaks):
        self._CurratAxeBraggList = peaks

        if peaks is None:
            self.plotCurratAxe = False
        elif np.array(self.plotCurratAxe).size ==0:
            self.plotCurratAxe = False
        elif self.rlu :
            self.monoPoints = self.axRLU.sample.CurratAxe(Ei=self.Ei,Ef=self.Ef,Bragg=self.CurratAxeBraggList,spurionType='Monochromator',HKL=False)[:,:,:,:2]
            self.anaPoints = self.axRLU.sample.CurratAxe(Ei=self.Ei,Ef=self.Ef,Bragg=self.CurratAxeBraggList,spurionType='Analyser',HKL=False)[:,:,:,:2]

            self.dE = np.repeat((self.Ei[:,np.newaxis]-self.Ef[np.newaxis])[np.newaxis],len(self.CurratAxeBraggList),axis=0)
            dEIndex = np.round((self.dE-self.bins[2][0,0,0])/self.dQE[2])
            dQxIndexMono = np.round((self.monoPoints[:,:,:,0]-self.bins[0][0,0,0])/self.dQE[0])
            dQyIndexMono = np.round((self.monoPoints[:,:,:,1]-self.bins[1][0,0,0])/self.dQE[1])
            dQxIndexAna = np.round((self.anaPoints[:,:,:,0]-self.bins[0][0,0,0])/self.dQE[0])
            dQyIndexAna = np.round((self.anaPoints[:,:,:,1]-self.bins[1][0,0,0])/self.dQE[1])

            self.CurratAxeIndicesMono = np.array([dQxIndexMono,dQyIndexMono,dEIndex])
            self.CurratAxeIndicesAna = np.array([dQxIndexAna,dQyIndexAna,dEIndex])
            try:
                self.plot()
                self.plotCurratAxe = True
            except:
                pass
            
        else:
            self.plotCurratAxe = False

    @CurratAxeBraggList.getter
    def CurratAxeBraggList(self):
        return self._CurratAxeBraggList

    @property
    def plotCurratAxe(self):
        return self._plotCurratAxe
    
    @plotCurratAxe.setter
    def plotCurratAxe(self,value):
        self._plotCurratAxe = value
        if hasattr(self,'_axes'):
            # for ax in self._axes:
            #     for plotter in ['BraggScatterMono','BraggScatterAna']:
            #         if hasattr(ax,plotter):
            #             try:
            #                 getattr(ax,plotter).remove()
            #             except ValueError:
            #                 pass
            if self._plotCurratAxe:
                self.plot()
    
    def plot(self):

        self.text.set_text(self.stringValue())
        try:
            self.im.set_array(self.emptyData)
        except TypeError:
            pass
        if self._axesChanged:
            tempData = np.ma.array(self.im.get_array().T)
            tempData.mask = np.ones_like(tempData,dtype=bool)
            self.im.set_array(tempData)
            self._axesChanged = False
        else:
            self.im.set_array(self.masked_array[:,:,int(self.value)].T.flatten())
        
        if not self.CurratAxeBraggList is None and not self.Ei is None and self.plotCurratAxe is True and self.rlu is True:

            insidePointsMono = self.monoPoints[np.isclose(self.CurratAxeIndicesMono[self.axis],self.Energy_slider.val)]
            insidePointsAna = self.anaPoints[np.isclose(self.CurratAxeIndicesAna[self.axis],self.Energy_slider.val)]
            if not hasattr(self.ax,'BraggScatterMono'):
                self.ax.BraggScatterMono = self.ax.plot([],[],zorder=20, linestyle="", marker="o", mfc='none', color='r',label='Currat-Axe Monochromator')[0]
            if not hasattr(self.ax,'BraggScatterAna'):
                self.ax.BraggScatterAna = self.ax.plot([],[],zorder=20, linestyle="", marker="o", mfc='none', color='k',label='Currat-Axe Analyzer')[0]
            
            
            self.ax.BraggScatterMono.set_data([],[])
            self.ax.BraggScatterAna.set_data([],[])
            
            if len(insidePointsMono) > 0: # Only if there are any points to plot, do it
                insidePointsMono = insidePointsMono.T
                
                                
                
                if self.axis == 2:
                    insidePointsAna = insidePointsAna.T
                    self.ax.BraggScatterMono.set_data(*insidePointsMono)
                    self.ax.BraggScatterAna.set_data(*insidePointsAna)
                else:
                    energy = self.dE[np.isclose(self.CurratAxeIndicesMono[self.axis],self.Energy_slider.val)]
                    
                    plotPointsMoni = np.array([insidePointsMono[1-self.axis],energy])
                    self.ax.BraggScatterMono.set_data(*plotPointsMoni)
                    
            if len(insidePointsAna) > 0 and self.axis != 2: # Only if there are any points to plot, do it
                insidePointsAna = insidePointsAna.T
                
                energy = self.dE[np.isclose(self.CurratAxeIndicesAna[self.axis],self.Energy_slider.val)]

                plotPointsAna = np.array([insidePointsAna[1-self.axis],energy])
                self.ax.BraggScatterAna.set_data(*plotPointsAna)
            #self.ax.

        self.im.set_clim(self.caxis)
        self.ax.set_position(self.figpos)
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.axis == 2:
            pass
        if self.grid:
            self.ax.grid(self.grid,zorder=self.gridZOrder)
        else:
            self.ax.grid(self.grid)
    def set_title(self,title):
        self.ax.set_title(title)

    def saveToFile(self,folder,extension,gui=None):
        """Save all planes to files
        
        Args:

            - folder (path): Folder in which to save files

            - extension (str): name of wanted extension. Must be accepted by matplotlib

        Kwargs:

            - gui (PyQt5 MainApplication): Used for MJOLNIR Gui
        """
        XLength,YLength,ZLength = self.Data.shape

        startProjection = self.axis
        startPosition = self.Energy_slider.val

        totalFigures = XLength+YLength+ZLength
        if not gui is None:
            gui.setProgressBarMaximum(totalFigures)
            progress = 0
            gui.setProgressBarValue(progress) 
        self.setProjection(2)
        for zidx in range(ZLength):
            self.setPlane(zidx)
            Energy = self.bins[2][0,0,zidx]
            self.figure.savefig(os.path.join(folder,'E{}.{}'.format('{:.2f}'.format(Energy).replace('.','p').replace('-','m'),extension)))
            if not gui is None:
                progress += 1
                gui.setProgressBarValue(progress)

        self.setProjection(0)

        for xidx in range(XLength):
            self.setPlane(xidx)
            position = self.bins[0][xidx,0,0]
            self.figure.savefig(os.path.join(folder,'X{}.{}'.format('{:.2f}'.format(position).replace('.','p').replace('-','m'),extension)))
            if not gui is None:
                progress += 1
                gui.setProgressBarValue(progress)
            
        self.setProjection(1)

        for yidx in range(YLength):
            self.setPlane(yidx)
            position = self.bins[1][0,yidx,0]
            self.figure.savefig(os.path.join(folder,'Y{}.{}'.format('{:.2f}'.format(position).replace('.','p').replace('-','m'),extension)))
            if not gui is None:
                progress += 1
                gui.setProgressBarValue(progress)
            
        self.setProjection(startProjection)
        self.setPlane(startPosition)


def eventdecorator(function,self,ax,event,*args,**kwargs):# pragma: no cover
    if event.xdata is not None and ax.in_axes(event) and ax.isActive:
        try:
            C = ax.get_figure().canvas.cursor().shape() # Only works for pyQt5 backend
        except:
            pass
        else:
            if C != 0 or ax.drawState != States.INACTIVE:
                return
        return function(self,event.xdata,event.ydata,*args,**kwargs)


def onclick(self,x,y,returnText=False, outputFunction=print,extra=None): # pragma: no cover
    idz = self.value
    axis = self.axis

    XX,YY = self.X[:,:,idz],self.Y[:,:,idz]
    XX = 0.25*(XX[:-1,:-1]+XX[1:,:-1]+XX[:-1,1:]+XX[1:,1:])
    YY = 0.25*(YY[:-1,:-1]+YY[1:,:-1]+YY[:-1,1:]+YY[1:,1:])
    idx = np.unravel_index(np.argmin(np.abs(XX-x)),XX.shape)[0]
    idy = np.unravel_index(np.argmin(np.abs(YY-y)),YY.shape)[1]
    I = self.masked_array[idx,idy,idz]

    masked = np.ma.is_masked(I)
    printString = ''
    printString+=self.ax.format_coord(x, y)+', '

    if masked:
        I = np.nan

    printString+='I = {:.4E}'.format(I)
    
    if self.allData is True and not masked:
        if self.axis == 0:
            flipper = [2,0,1]
        elif self.axis == 1:
            flipper = [0,2,1]
        else:
            flipper = [0,1,2]
        ID = np.array([idx,idy,idz])[flipper]
        cts = self.Counts[ID[0],ID[1],ID[2]]
        Norm = self.Normalization[ID[0],ID[1],ID[2]]
        Mon = self.Monitor[ID[0],ID[1],ID[2]]
        NC = self.NormCounts[ID[0],ID[1],ID[2]]
        printString+=', Cts = {:}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
    if not extra is None:
        printString+=extra
    if returnText:
        return printString
    else:
        if hasattr(self.ax,'suppressPrint'):
            if not self.ax.suppressPrint:
                outputFunction(printString)
        else:
            outputFunction(printString)
        
def onkeypress(event,self): # pragma: no cover
    if event.key in Viewer3DSettings['upwards']:
        increaseAxis(event,self)
    elif event.key in Viewer3DSettings['downwards']:
        decreaseAxis(event,self)
    elif event.key in Viewer3DSettings['home']:
        self.Energy_slider.set_val(self.Energy_slider.valmin)
    elif event.key in Viewer3DSettings['end']:
        self.Energy_slider.set_val(self.Energy_slider.valmax)

    elif event.key in Viewer3DSettings['QxE']:
        if self.axis!=0:
            reloadslider(self,0)
            del self.im
            #self.currentData = None
            self.ax.grid(False)
            if self.shading=='flat':
                self.im = self.ax.pcolormesh(self.X[:,:,0].T,self.Y[:,:,0].T,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=self.cmap)
            elif self.shading=='gouraud':
                self.im = self.ax.pcolormesh(0.5*(self.X[:-1,:-1,0]+self.X[1:,1:,0]).T,0.5*(self.Y[:-1,:-1,0]+self.Y[1:,1:,0]).T,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=self.cmap) # ,vmin=1e-6,vmax=6e-6
            else:
                raise AttributeError('Did not understand shading {}.'.format(self.shading))
            if self.grid:
                self.ax.grid(self.grid,zorder=self.gridZOrder)
            else:
                self.ax.grid(self.grid)
            self.im.set_clim(self.caxis)
            self.Energy_slider.set_val(0)
            self.plot()
            self.ax.set_xlim([np.min(self.X),np.max(self.X)])
            self.ax.set_ylim([np.min(self.Y),np.max(self.Y)])
    elif event.key in Viewer3DSettings['QyE']:
        if self.axis!=1:
            reloadslider(self,1)
            del self.im
            #self.currentData = None
            self.ax.grid(False)
            if self.shading=='flat':
                self.im = self.ax.pcolormesh(self.X[:,:,0].T,self.Y[:,:,0].T,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=self.cmap)
            elif self.shading=='gouraud':
                self.im = self.ax.pcolormesh(0.5*(self.X[:-1,:-1]+self.X[1:,:1:]).T,0.5*(self.Y[:-1,-1]+self.Y[1:,1:]).T,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=self.cmap) # ,vmin=1e-6,vmax=6e-6
            else:
                raise AttributeError('Did not understand shading {}.'.format(self.shading))
            if self.grid:
                self.ax.grid(self.grid,zorder=self.gridZOrder)
            else:
                self.ax.grid(self.grid)
            self.im.set_clim(self.caxis)
            self.Energy_slider.set_val(0)
            self.plot()
            self.ax.set_xlim([np.min(self.X),np.max(self.X)])
            self.ax.set_ylim([np.min(self.Y),np.max(self.Y)])
    elif event.key in Viewer3DSettings['QxQy']:
        if self.axis!=2:
            reloadslider(self,2)
            del self.im
            #self.currentData = None
            self.ax.grid(False)
            if self.shading=='flat':
                self.im = self.ax.pcolormesh(self.X[:,:,0].T,self.Y[:,:,0].T,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=self.cmap)
            elif self.shading=='gouraud':
                XX = 0.5*(self.X[:-1,:-1,self.value]+self.X[1:,1:,self.value]).T
                YY = 0.5*(self.Y[:-1,:-1,self.value]+self.Y[1:,1:,self.value]).T
                self.im = self.ax.pcolormesh(XX,YY,self.masked_array[:,:,self.value].T,zorder=10,shading=self.shading,cmap=self.cmap) # ,vmin=1e-6,vmax=6e-6
            else:
                raise AttributeError('Did not understand shading {}.'.format(self.shading))
            if self.grid:
                self.ax.grid(self.grid,zorder=self.gridZOrder)
            else:
                self.ax.grid(self.grid)
            self.im.set_clim(self.caxis)
            self.Energy_slider.set_val(0)
            self.plot()
            self.ax.set_xlim([np.min(self.X),np.max(self.X)])
            self.ax.set_ylim([np.min(self.Y),np.max(self.Y)])


def reloadslider(self,axis):
    self.setAxis(axis)
    self.Energy_slider.set_val(0)
    self.Energy_slider.label.remove()
    self.Energy_slider.disconnect_events()#(self.Energy_slider.cids[0])
    self.Energy_slider.vline.set_visible(False)
    
    del self.Energy_slider
    
    zeroPoint = np.argmin(np.abs(0.5*(self.Z[0,0][1:]+self.Z[0,0][:-1])))
    self.Energy_slider = Slider(self.Energy_slider_ax, label=self.label, valmin=self.lowerLim, valmax=self.upperLim, valinit=zeroPoint)
    self.Energy_slider.valtext.set_visible(False)
    self.Energy_slider.on_changed(lambda val: sliders_on_changed(self,val))
    self.value=0
    self.im.remove()
    
        
def onscroll(event,self): # pragma: no cover
    if(event.button in Viewer3DSettings['upwardsScroll']):
        increaseAxis(event,self)
    elif event.button in Viewer3DSettings['downwardsScroll']:
        decreaseAxis(event,self)


def increaseAxis(event,self): # pragma: no cover
    self.Energy_slider.set_val(self.Energy_slider.val+1)
    
    
def decreaseAxis(event,self): # pragma: no cover
    self.Energy_slider.set_val(self.Energy_slider.val-1)   
    

def sliders_on_changed(self,val): # pragma: no cover
    value = int(np.round(val))
    
    if value>self.Energy_slider.valmax:
        self.Energy_slider.set_val(self.Energy_slider.valmax)
        
    elif value<self.Energy_slider.valmin:
        self.Energy_slider.set_val(self.Energy_slider.valmin)
    
    if value<=self.Energy_slider.valmax and value>=self.Energy_slider.valmin:
       
        if hasattr(self.ax,'updateLimits'):
            self.ax.updateLimits(self.Z[0,0,value],self.Z[0,0,value+1])

        if value!=val:
            self.value = val
            self.Energy_slider.set_val(value)
            self.plot()
        else:
            self.value = val
            #self.Energy_slider.set_val(value)
            self.plot()
    if hasattr(self.ax,'_step'):
        self.ax._step=self.calculateValue()


def addColorbarSliders(self,c_min,c_max,c_minval,c_maxval,ax_cmin,ax_cmax,log=True):
    """Add two colorbars controling the colour axis

    args:

        self (Viewer3D object): Current object

        c_min (float): Minimal color value

        c_max (float): Maximal color value

        c_minval (float): Starting value of lower bound

        c_maxval (float): Starting value of upper  bound

        ax_cmin (mpl axes): Axis in which lower bound slider is to be shown
        
        ax_cmax (mpl axes): Axis in which upper bound slider is to be shown
        
    Kwargs:
    

        log (bool): If true, sliders are logarithmic

    """

    fig = self.ax.get_figure()
    if log==True:
        fig.s_cmin = Slider(ax_cmin, 'min', np.log10(c_min+1e-20), np.log10(c_max), valinit=np.log10(c_min+1e-20),orientation='vertical',valfmt='%2.1f')
        fig.s_cmax = Slider(ax_cmax, 'max', np.log10(c_min+1e-20), np.log10(c_max), valinit=np.log10(c_max),orientation='vertical',valfmt='%2.1f')
    else:
        fig.s_cmin = Slider(ax_cmin, 'min', c_min, c_max, valinit=c_min,orientation='vertical',valfmt='%2.1e')
        fig.s_cmax = Slider(ax_cmax, 'max', c_min, c_max, valinit=c_max,orientation='vertical',valfmt='%2.1e')

    fig.s_cmin._log = log
    fig.s_cmax._log = log
    for s in [fig.s_cmin,fig.s_cmax]:
        s.label.set_fontsize(12)
        s.valtext.set_fontsize(12)
    
    def update(fig,val, bar=None, s=None,log=log):
        _cmin =fig.s_cmin.val# np.log10(fig.s_cmin.val)
        _cmax =fig.s_cmax.val# np.log10(fig.s_cmax.val)
        if _cmin>_cmax:
            if bar == 'min': # If lower bar is change push max upwards
                _cmax = _cmin
                if log:
                    fig.s_cmax.set_val(np.power(10,_cmax))
                else:
                    fig.s_cmax.set_val(_cmax)
            else:
                _cmin = _cmax
                if log:
                    fig.s_cmin.set_val(np.power(10,_cmin))
                else:
                    fig.s_cmin.set_val(_cmin)
        
        if log:
            
            #self.caxis = (np.power(10,_cmin), np.power(10,_cmax))
            self.im.set_clim([np.power(10,_cmin), np.power(10,_cmax)])
            fig.s_cmax.valtext.set_text(_cmin)
            fig.s_cmax.valtext.set_text(_cmax)
            self._caxis = (np.power(10,_cmin), np.power(10,_cmax))
        else:
            
            #self.caxis = (_cmin,_cmax)
            self.im.set_clim([_cmin,_cmax])
            self._caxis = (_cmin,_cmax)
        plt.draw()
        self.colorbar.update_normal(self.im)
    
    
    fig.s_cmin.on_changed(lambda val,*arg,**kwargs: update(fig,val,*arg,bar='min',**kwargs))
    fig.s_cmax.on_changed(lambda val,*arg,**kwargs: update(fig,val,*arg,bar='max',**kwargs))
    
     
    
   
    fig._savefig = fig.savefig
    

    @functools.wraps(fig.savefig)
    def savefig(fname,hide = [ax_cmin,ax_cmax],fig=fig,**kwargs):
    
        unhide = []
        for obj in hide:
            if obj.get_visible():
                obj.set_visible(False)
                unhide.append(obj)
            
        fig._savefig(fname,**kwargs)
        for obj in unhide:
            obj.set_visible(True)
            
        fig.canvas.draw()
    
    fig.savefig = savefig

def cut1DFunctionDefault(self,dr):# pragma: no cover
    global cut1DHolder
    parameters = extractCut1DProperties(dr.rect,self.ax.sample)
    step = self.dQE[self.axis]
    pos = np.zeros(3,dtype=int)
    pos[self.axis]=self.Energy_slider.val
    
    EMin = self.bins[self.axis][pos[0],pos[1],pos[2]]
    EMax = EMin+step
    cut1DHolder.append([self.ds.plotCut1D(**parameters,Emin=EMin,Emax=EMax)])

