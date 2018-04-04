import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from matplotlib.widgets import Slider

import warnings

class Viewer3D(object):
    def __init__(self,Data,bins,axis=2):
        """3 dimensional viewing object generating interactive Matplotlib figure. Keeps track of all the different plotting functions and variables in order to allow
        the user to change between different slicing modes and to scroll through the data in an interactive way.
        """
        self.Data = Data
        self.bins = bins
        
        self.figure = plt.figure()
        self.ax = self.figure.add_subplot(111)
        
       
        self.value = 0
        self.figure.subplots_adjust(bottom=0.25)
        
        self.cmap = cm.jet
        self.cmap.set_bad('white',1.)
        self.value = 0
        
        axis_color='white'
        self.setAxis(axis)
        
        self.figure.canvas.mpl_connect('key_press_event',lambda event: onkeypress(event, self) )
        self.figure.canvas.mpl_connect('scroll_event',lambda event: onscroll(event, self))
        
        zeroPoint = np.argmin(np.abs(self.Z))
        
    
        self.Energy_slider_ax = self.figure.add_axes([0.15, 0.1, 0.7, 0.03], facecolor=axis_color)
        self.Energy_slider = Slider(self.Energy_slider_ax, label=self.label, valmin=self.lowerLim, valmax=self.upperLim, valinit=zeroPoint,valfmt='%0f')
        self.Energy_slider.valtext.set_visible(False)
        
        self.Energy_slider.on_changed(lambda val: sliders_on_changed(self,val))
            
            
        textposition = [self.Energy_slider_ax.get_position().p1[0]+0.005,self.Energy_slider_ax.get_position().p0[1]+0.005]
        self.text = self.figure.text(textposition[0], textposition[1],s=self.stringValue())
        
        
        self.im = self.ax.imshow(self.masked_array[:,:,self.value].T,cmap=self.cmap,extent=[self.X[0],self.X[-1],self.Y[0],self.Y[-1]],origin='lower')
        self.cbaxes = self.figure.add_axes([0.8, 0.2, 0.03, 0.7])
        self.colorbar = self.figure.colorbar(self.im,cax = self.cbaxes)
        warnings.simplefilter("ignore")
        self.figure.tight_layout(rect=[0,0.1,0.9,0.9])
        warnings.simplefilter("once")

        self.text.set_text(self.stringValue())
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.axis == 2:
            self.ax.set_xlim(np.min([xlim[0],ylim[0]]),np.max([xlim[1],ylim[1]]))
        self.Energy_slider.set_val(self.value)

        cid = self.figure.canvas.mpl_connect('button_press_event', onclick)

    def setAxis(self,axis):
        if axis==2:
            axes = (0,1,2)
            self.ax.set_xlabel('Qx [1/AA]')
            self.ax.set_ylabel('Qy [1/AA]')
            label = 'Energy'
        elif axis==1:
            axes = (0,2,1)
            self.ax.set_xlabel('Qx [1/AA]')
            self.ax.set_ylabel('E [meV]')
            label = 'Qy'
        elif axis==0:
            axes = (1,2,0)
            self.ax.set_xlabel('Qy [1/AA]')
            self.ax.set_ylabel('E [meV]')
            label = 'Qx'
        else:
            raise AttributeError('Axis provided not recognized. Should be 0, 1, or 2 but got {}'.format(axis))
        
        X=self.bins[axes[0]]
        Y=self.bins[axes[1]]
        Z=self.bins[axes[2]]

        masked_array = np.ma.array (self.Data, mask=np.isnan(self.Data)).transpose(axes)
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

        
    def stringValue(self):
        if self.axis==2:
            unit = ' meV'
        else:
            unit = ' 1/AA'
        
        return str(np.round(self.Z[self.value],2))+unit
    
    
    def plot(self):
        self.im = self.ax.imshow(self.masked_array[:,:,self.value].T,cmap=self.cmap,extent=[self.X[0],self.X[-1],self.Y[0],self.Y[-1]],origin='lower')
        self.colorbar.update_bruteforce(self.im)
        self.text.set_text(self.stringValue())
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.axis == 2:
            self.ax.set_xlim(np.min([xlim[0],ylim[0]]),np.max([xlim[1],ylim[1]]))


def onclick(event):
    if event.xdata is not None:
        print('x={}, y={}, xdata={}, ydata={}'.format(event.x, event.y, event.xdata, event.ydata))



def onkeypress(event,self):
        
        if event.key in ['+','up']:
            increaseAxis(event,self)
        elif event.key in ['-','down']:
            decreaseAxis(event,self)
        elif event.key in ['0']:
            if self.axis!=0:
                reloadslider(self,0)
                del self.im
                self.plot()
                warnings.simplefilter("ignore")
                self.figure.tight_layout(rect=[0,0.1,0.9,0.9])
                warnings.simplefilter("once")
        elif event.key in ['1']:
            if self.axis!=1:
                reloadslider(self,1)
                del self.im
                self.plot()
                warnings.simplefilter("ignore")
                self.figure.tight_layout(rect=[0,0.1,0.9,0.9])
                warnings.simplefilter("once")
        elif event.key in ['2']:
            if self.axis!=2:
                reloadslider(self,2)
                del self.im
                self.plot()
                warnings.simplefilter("ignore")
                self.figure.tight_layout(rect=[0,0.1,0.9,0.9])
                warnings.simplefilter("once")


def reloadslider(self,axis):
    self.setAxis(axis)
    self.Energy_slider.label.remove()#self.Energy_slider.label.text('')
    self.Energy_slider.disconnect(self.Energy_slider.cids[0])
    self.Energy_slider.vline.set_visible(False)
    self.Energy_slider.set_val(0)
    del self.Energy_slider
    zeroPoint = np.argmin(np.abs(self.Z))
    self.Energy_slider = Slider(self.Energy_slider_ax, label=self.label, valmin=self.lowerLim, valmax=self.upperLim, valinit=zeroPoint)
    self.Energy_slider.valtext.set_visible(False)
    self.Energy_slider.on_changed(lambda val: sliders_on_changed(self,val))
    self.value=0
    self.Energy_slider.set_val(self.value)
    
        
def onscroll(event,self):
    if(event.button=='up'):
        increaseAxis(event,self)
    elif event.button=='down':
        decreaseAxis(event,self)





def increaseAxis(event,self):
    self.Energy_slider.set_val(self.Energy_slider.val+1)
    
    
def decreaseAxis(event,self):
    self.Energy_slider.set_val(self.Energy_slider.val-1)   
    

def sliders_on_changed(self,val):
        value = int(np.round(val))
        
        if value>self.Energy_slider.valmax:
            self.Energy_slider.set_val(self.Energy_slider.valmax)
            
        elif value<self.Energy_slider.valmin:
            self.Energy_slider.set_val(self.Energy_slider.valmin)
        
        if value<=self.Energy_slider.valmax and value>=self.Energy_slider.valmin:
            if value!=val:
                self.Energy_slider.set_val(value)
                
            else:
                self.value = val
                self.plot()

