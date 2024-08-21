# -*- coding: utf-8 -*-

import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')


import numpy as np
import pickle as pickle
import matplotlib.pyplot as plt
from MJOLNIR import _tools
from MJOLNIR.Data import Mask

import functools
import warnings



# If an update has happened to the self.int, and the originalDataSet has a backgroundIntensities pointer array
def recalculateOnChange(func):
    @functools.wraps(func)
    def recalculateFunction(self,*args,**kwargs):
        if self.isChanged and hasattr(self.originalDataSet,'_backgroundIntensities'):
            #print('RECALCULATING!!')
            self.generateFullBackground()
            #print('Done',f'{self.isChanged=}', flush=True)
        #else:
        #    print('No recalculation needed')
        return func(self,*args,**kwargs)
    return recalculateFunction
        

class BackgroundObject(object):
    
    def __init__(self,originalDataSet,dQ=0.025, dE = 0.05, backgroundMask=None):
        """Object to stream line the background generation
        
        Args:
            
            - originalDataSet (DataSet): DataSet containing the data set used for the background
            
        Kwargs:
            
            - dQ (float): Step size along the q direction in the powder average bin in 1/AA (default 0.02)
            
            - dE (float): Step size along the energy direction in the powder average bin in meV (default 0.1)
            
            - backgroundMask (Mask): Masking object to utilize for background generation. If none use DataSet mask (default None)
        
        """
        
        self.originalDataSet = originalDataSet
        
        self.dQ = dQ
        self.dE = dE
        self.backgroundMask = backgroundMask
        
        # Flag to keep track of updates in model not yet applied
        self.isChanged = False
        
        
    def generatePowderBackground(self,dQ=None, dE = None):
        if not dQ is None:
            self.dQ = dQ
        if not dE is None:
            self.dE = dE
        
        ## Needed to reset data set after calculating background
        foregroundMask = self.originalDataSet.mask
        
        if self.backgroundMask is None:
            self.backgroundMask = foregroundMask
        else:
            if isinstance(self.backgroundMask,(Mask.MaskingObject)): # Is mask
                self.backgroundMask = self.backgroundMask(self.originalDataSet)
            
            self.originalDataSet.mask = self.backgroundMask
        # Find Q and E ranges
        
        qLength = np.linalg.norm([self.originalDataSet.qx.extractData(),self.originalDataSet.qy.extractData()],axis=0)
        
        QMax = np.nanmax(qLength) # QMin is always set to 0
        #print(QMax)
        Energies = self.originalDataSet.energy.extractData()
        
        EMin = np.nanmin(Energies)
        EMax = np.nanmax(Energies)#self.originalDataSet.energy.extractData().max()
                
        #print(EMin,EMax)
        I = self.originalDataSet.I.extractData()
        Monitor = self.originalDataSet.Monitor.extractData()
        Norm = self.originalDataSet.Norm.extractData()
        
        # Position in the powder average (QLength, Energy)
        positions2D = np.array([np.linalg.norm([self.originalDataSet.qx.extractData(),
                                                self.originalDataSet.qy.extractData()],axis=0),
                                self.originalDataSet.energy.extractData()])
    
        ## Reset mask 
        self.originalDataSet.mask = foregroundMask
        
        
        # Generate the bins with suitable extensions
        self.QBins = np.arange(0,QMax+self.dQ*1.1,self.dQ)
        self.EnergyBins = np.arange(EMin-self.dE*1.1,EMax+self.dE*1.1,self.dE)
    
        
        # Perform 2D histogram
        self._intensity,self._norm,self._monitor,self._counts = _tools.histogramdd(positions2D.T,bins=(self.QBins,self.EnergyBins),weights=[I,Norm,Monitor],returnCounts=True)
        
        
        
        ## Calcualte the intensities
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self._int = np.divide(self.intensity*self.counts,self.monitor*self.norm)
        self.isChanged = True
        
    @recalculateOnChange
    def sample(self,positions,norm=None,monitor=None):#,returnAll=False):
        """Sample the background at the given positions (|Q|,E)
        
        Args:
            
            - positions (2D list): List of positions (|Q|, E) with shape N x 2
            
        Kwargs:
            
            - norm (list): List of normalizations of shape N (default None)
            
            - monitor (list): List of monitor of shape N (default None)
            
            - returnAll (bool): If True, returns intensity, monitor, norm, counts
            
            
        """
        BGBins = [self.QBins,self.EnergyBins]
            
        
        xQ = _tools.findFlattenedIndex(positions[0].reshape(-1,1), [BGBins[0]])-1
        yE = _tools.findFlattenedIndex(positions[1].reshape(-1,1), [BGBins[1]])-1

        inside = np.logical_and(xQ<len(self.QBins)-1,yE<len(self.EnergyBins)-1)
        newI = np.zeros_like(inside,dtype=float)
        
        
        
        newI[inside] = self.int[xQ[inside],yE[inside]]
        if not norm is None:
            newI[inside]*=norm.reshape(newI.shape)[inside]
        if not monitor is None:
            newI[inside]*=monitor.reshape(newI.shape)[inside]
        return newI
    
    def generateFullBackground(self):
        """Generate the full fetched background data set from the original one"""
        
        self.isChanged = False
        for df in self.originalDataSet:
            positions = np.asarray([np.linalg.norm([df.qx.flatten(),df.qy.flatten()],axis=0),df.energy.flatten()])
            df.backgroundIntensities = self.sample(positions,norm=df.Norm.flatten(),monitor=df.Monitor.flatten()).reshape(df.I.shape)
        
        files = self.originalDataSet.convertedFiles if len(self.originalDataSet.convertedFiles) > 0 else self.originalDataSet.dataFiles
        self.originalDataSet.backgroundIntensities = _tools.PointerArray('backgroundIntensities',files)
        
        
        
        
    @property
    def dQ(self):
        return self._dQ
    
    @dQ.getter
    def dQ(self):
        return self._dQ
    
    @dQ.setter
    def dQ(self,value):
        if value > 0:
            self._dQ = value
        else:
            raise AttributeError('Provided dQ is negative or null. Received {:}'.format(value))
    
    @property
    def dE(self):
        return self._dE
    
    @dE.getter
    def dE(self):
        return self._dE
    
    @dE.setter
    def dE(self,value):
        if value > 0:
            self._dE = value
        else:
            raise AttributeError('Provided dE is negative or null. Received {:}'.format(value))
            
    @property
    def int(self):
        if not hasattr(self,'_int'):
            warnings.warn('No background was found. Generating using standard parameters...')
            self.generatePowderBackground()
        return self._int
    
    
    
    @property
    def intensity(self):
        if not hasattr(self,'_intensity'):
            warnings.warn('No background was found. Generating using standard parameters...')
            self.generatePowderBackground()
        return self._intensity
    
    @property
    def monitor(self):
        if not hasattr(self,'_monitor'):
            warnings.warn('No background was found. Generating using standard parameters...')
            self.generatePowderBackground()
        return self._monitor
    
    @property
    def norm(self):
        if not hasattr(self,'_norm'):
            warnings.warn('No background was found. Generating using standard parameters...')
            self.generatePowderBackground()
        return self._norm
    
    @property
    def counts(self):
        if not hasattr(self,'_counts'):
            warnings.warn('No background was found. Generating using standard parameters...')
            self.generatePowderBackground()
        return self._counts
    
    
    
    
    def plotPowderAverage(self,ax=None):
        if ax is None:
            fig,ax = plt.subplots()
            
        if not hasattr(self,'Int'):
            self.generatePowderBackground()
            
        
        X,Y = np.meshgrid(self.QBins,self.EnergyBins)
        ax.p = ax.pcolormesh(X,Y,self.int.T)
        ax.set_xlabel('|Q| [1/A]')
        ax.set_ylabel('E [meV]')
        ax.get_figure().colorbar(ax.p)
        
        def set_clim(ax,vmin,vmax):
            ax.p.set_clim(vmin,vmax)
            
        ax.set_clim = lambda vmin,vmax: set_clim(ax,vmin,vmax)
        
        return ax