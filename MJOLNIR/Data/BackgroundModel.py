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

try: 
    from AMBER import background
    AMBER_located = True
except ModuleNotFoundError:
    AMBER_located = True


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
            
        #if not hasattr(self,'_int'):
        #    self.generatePowderBackground()
            
        
        X,Y = np.meshgrid(self.QBins,self.EnergyBins)
        ax.p = ax.pcolormesh(X,Y,self.int.T)
        ax.set_xlabel('|Q| [1/A]')
        ax.set_ylabel('E [meV]')
        ax.get_figure().colorbar(ax.p)
        
        def set_clim(ax,vmin,vmax):
            ax.p.set_clim(vmin,vmax)
            
        ax.set_clim = lambda vmin,vmax: set_clim(ax,vmin,vmax)
        
        return ax

class AMBERBackground(BackgroundObject):
    def __init__(self, dataset, dQx = 0.03, dQy = 0.03, dE = 0.05, l = None, mu = None, beta = 1.0, backgroundMask=None):
        """
        Create an AMBER background model

        Args: 
            - dataset: MJOLNIR.DataSet
                DataSet to be utilized for the background

        Kwargs: 
            - dQx: float
                Voxel binning along the Qx direction in Ångstrom (default 0.03)
            - dQy: float
                Voxel binning along the Qy direction in Ångstrom (default 0.03)
            - dE: float
                Voxel binning along the energy direction in meV (default 0.05)
            - l: float
                Cost function parameter concerning signal sparcity. If None, found 
                through the Median Absolute Deviation (MAD) (default None)
            - mu: float
                Cost function parameter concerning signal changes along the energy
                direction. If None, mu is estiamted through the magnitude of the signal 
                variance along energy (default None)
            - beta: float
                Cost function parameter concerning background smoothness. Can be estimated through,
                a cross validation method. (default 1.0)
            - backgroundMask: MJOLNIR.Mask
                Mask to be applied to the data set before AMBER is run (default None)
        """
        if not AMBER_located:
            raise ModuleNotFoundError("The AMBER package was not found. This can be install by\npython -m pypi AMBER-ds4ms")
        
        dQ = np.mean([dQx, dQy])
        super(AMBERBackground,self).__init__(originalDataSet = dataset, dQ=dQ, dE = dE)
        
        self.originalDataSet = dataset

        self.dQx = dQx
        self.dQy = dQy

        self.data, self.bins = self.originalDataSet.binData3D(self.dQx, self.dQy, self.dE)

        # Compute the normalized intensity (I) using binned data
        self.I = np.divide(self.data[0] * self.data[-1], self.data[1] * self.data[2])

        # Extract coordinates in (Qx, Qy, E) dimensions
        self.Qx = self.bins[0][:-1, 0, 0]
        self.Qy = self.bins[1][0, :-1, 0]
        self.E = self.bins[2][0, 0, :-1]

        self.l = l
        self.mu = mu
        self.beta = beta

        self.dtype = self.I.dtype


    def set_radial_bins(self, max_radius=None, n_bins=None):
        """
        Define ranges for radial bins.

        Args:
            - max_radius: float
                Maximum radial distance for the bins. If None found from DataSet (default None)
            - n_bins: int
                Number of bins to create. If None found from max_radius and largest of dQx and dQy (default None)
        """
        if max_radius is None:
            local_max = 0
            for df in self.originalDataSet:
                local_max = np.max([local_max,np.max(np.linalg.norm([df.qx,df.qy],axis=0))])
            max_radius = local_max
        # Set the maximum radial distance and the number of bins
        self.max_radius = max_radius

        if n_bins is None:
            dQ = np.max([self.dQx, self.dQy])
            n_bins = int(np.ceil(self.max_radius/dQ))
        self.n_bins = n_bins

        # Generate radial bins using torch linspace
        # The range starts from 0 and goes up to the square root of max_radius, creating n_bins+1 points
        self.r_range = np.linspace(0, self.max_radius, self.n_bins, dtype=self.dtype)
        self.dQ = np.mean(np.diff(self.r_range))
        self.QBins = self.r_range
        self.EnergyBins = self.E

    def generateAMBER(self, n_epochs=20, verbose=True):
        """
        Run AMBER.

        Kwargs: 
            - n_epocs: int
                Maximum number of epochs AMBER runs (default 20)
            - verbose: bool
                If set True, Display convergence logs (default True)
        """
        self.AMBER = background.background(dtype=self.dtype)
        self.AMBER.set_gridcell_size(dqx=self.dQx, dqy=self.dQy, dE=self.dE)
        self.AMBER.set_binned_data(Qx=self.Qx, Qy=self.Qy, E=self.E, Int=self.I)
        self.AMBER.set_radial_bins(max_radius=self.max_radius, n_bins=self.n_bins)
        self.AMBER.set_variables()

        if self.l is None:
            self.l = self.AMBER.MAD_lambda()

        if self.mu is None:
            self.mu = self.AMBER.mu_estimator()

        self.AMBER.denoising(self.AMBER.Ygrid, lambda_=self.l, beta_=self.beta, mu_=self.mu, \
                             n_epochs=n_epochs, verbose=verbose)

    def applyAMBER(self):
        self._int = self.AMBER.b.T
        self.generateFullBackground()
        
        self.originalDataSet.backgroundModel = self

    def cross_validation(self, q=0.75, beta=None, l=None, mu=None, n_epochs=15, verbose=True):
        """
        Run a cross-validation for a specified set of parameters

        Kwargs: 
            - q: floay
                quantile level (default 0.75)
            beta: list of floats
                Parameter value beta, if None use previously given beta (default None)
            l: float
                Parameter value for lambda, None use previously given beta (default None)
            mu: float
                Parameter value for mu, None use previously given beta (default None)
            - n_epocs: int
                Maximum number of epochs AMBER runs (default 20)
            - verbose: bool
                If set True, Display convergence logs (default True)
        """

        if not hasattr(self,'AMBER'):
            self.AMBER = background.background(dtype=self.dtype)
            self.AMBER.set_gridcell_size(dqx=self.dQx, dqy=self.dQy, dE=self.dE)
            self.AMBER.set_binned_data(Qx=self.Qx, Qy=self.Qy, E=self.E, Int=self.I)
            self.AMBER.set_radial_bins(max_radius=self.max_radius, n_bins=self.n_bins)
            self.AMBER.set_variables()

        if l is None:
            self.l = self.AMBER.MAD_lambda()
            l = self.l
        if mu is None:
            self.mu = self.AMBER.mu_estimator()
            mu = self.mu
        if beta is None:
            beta = self.beta
        beta = np.asarray(beta)

        return self.AMBER.cross_validation(q=q, beta_range=beta, lambda_=l, mu_=mu, n_epochs=n_epochs, verbose=verbose)
