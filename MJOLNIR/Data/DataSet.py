# -*- coding: utf-8 -*-
import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')

import datetime
import h5py as hdf
import numpy as np
import pickle as pickle
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection,PolyCollection
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon
from MJOLNIR.Data import Viewer3D,RLUAxes
from MJOLNIR.Data import Mask
import MJOLNIR.Data.DataFile
import MJOLNIR.Data.Sample
from MJOLNIR import _tools
from mpl_toolkits.axisartist.grid_helper_curvelinear import \
    GridHelperCurveLinear
from mpl_toolkits.axisartist import SubplotHost
import pandas as pd
import pytest
from scipy.ndimage import filters
import scipy.optimize
from scipy.spatial import Voronoi,ConvexHull,KDTree
# from shapely.geometry import Polygon as PolygonS
# from shapely.geometry import Point as PointS
# from shapely.vectorized import contains
import time
from . import Viewer3DPyQtGraph
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import warnings
from ufit import Dataset

pythonVersion = sys.version_info[0]

_cache = []

class DataSet(object):
    @_tools.KwargChecker(include=['Author']) 
    def __init__(self, dataFiles=None, normalizationfiles=None, 
                 calibrationfiles=None, convertedFiles=None, **kwargs):
        """DataSet object to hold all informations about data.
        
        Kwargs:
            
            - dataFiles (string, DataFile or list of strings or DataFiles): List of datafiles or DataFile objects to be used in conversion (default None).

            - normalizationfiles (string or list of strings): Location of Vanadium normalization file(s) (default None).

            - calibrationfiles (string or list of strings): Location of calibration normalization file(s) (default None).

            - convertedFiles (string, DataFile or list of strings): Location of converted data files (default None).

        Raises:

            - ValueError
            
            - NotImplementedError
        
        """
        
        self._dataFiles = []
        self._normalizationfiles = []
        self._convertedFiles = []
        self._calibrationfiles = []
        self._mask = False
        self.index = 0
        self.absolutNormalized = 0 # Float to keep track of normalization has taken place (0 means false)


        if dataFiles is not None:
            self.dataFiles = dataFiles
            self._getData()

        if normalizationfiles is not None:
            self.normalizationfiles = normalizationfiles
        
        if convertedFiles is not None:
            self.convertedFiles = convertedFiles
            self._getData()

        if calibrationfiles is not None:
            self.calibrationfiles = calibrationfiles


        self._settings = {}
            
        
        if len(self.convertedFiles)!=0:
            self.sample = [d.sample for d in self]
        elif len(self.dataFiles)!=0:
            self.sample = [d.sample for d in self]

        # Add all other kwargs to settings
        for key in kwargs:
            if hasattr(kwargs,'shape'): # They are an array
                if kwargs[key].shape in [(),(1,)]:
                    self.settings[key] = float(kwargs[key])
            else:
                self.settings[key]=kwargs[key]
        
    @property
    def dataFiles(self):
        return self._dataFiles

    @dataFiles.getter
    def dataFiles(self):
        return self._dataFiles

    @dataFiles.setter
    def dataFiles(self,dataFiles):
        try:
            correctDataFiles = isListOfDataFiles(dataFiles)
            [self._dataFiles.append(file) for file in correctDataFiles if file.type in MJOLNIR.Data.DataFile.supportedRawFormats]
            [self._convertedFiles.append(file) for file in correctDataFiles if file.type in MJOLNIR.Data.DataFile.supportedConvertedFormats]
        except Exception as e:
            raise


    @property
    def normalizationfiles(self):
        return self._normalizationfiles

    @normalizationfiles.getter
    def normalizationfiles(self):
        return self._normalizationfiles

    @normalizationfiles.setter
    def normalizationfiles(self,normalizationfiles):
        try:
            self._normalizationfiles = isListOfStrings(normalizationfiles)
        except Exception as e:
            raise(e)


    @property
    def convertedFiles(self):
        return self._convertedFiles

    @convertedFiles.getter
    def convertedFiles(self):
        return self._convertedFiles

    @convertedFiles.setter
    def convertedFiles(self,convertedFiles):
        try:
            correctDataFiles = isListOfDataFiles(convertedFiles)
            [self._dataFiles.append(file) for file in correctDataFiles if file.type=='hdf']
            [self._convertedFiles.append(file) for file in correctDataFiles if file.type=='nxs']
        except Exception as e:
            raise(e)
        self._getData()


    @property
    def calibrationfiles(self):
        return self._calibrationfiles

    @calibrationfiles.getter
    def calibrationfiles(self):
        return self._calibrationfiles

    @calibrationfiles.setter
    def calibrationfiles(self,calibrationfiles):
        try:
            self._calibrationfiles = isListOfStrings(calibrationfiles)
        except Exception as e:
            raise(e)

    @property
    def mask(self):
        return self._mask

    @mask.getter
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self,mask):
        if isinstance(mask,list):
            if isinstance(mask[0],Mask.MaskingObject): # If list of maskobjects is provided
                if not len(mask) == len(self): #not the same number of masks as data files
                    raise AttributeError('Provided number of masks ({}) does not match number of data files ({})'.format(len(mask),len(self)))
                self._maskingObject = mask
                for m,df in zip(mask,self):
                    df.mask = m
                masksum = -1
            elif hasattr(self,'_maskingObject'): # not masking object but self has one, delete it
                del self._maskingObject
            masksum = np.sum([np.sum(x) for x  in mask])
            self.maskIndices = np.cumsum([np.sum(1-M) for M in mask])[:-1]

        elif isinstance(mask,Mask.MaskingObject):
            self._maskingObject = mask
            m = []
            for df in self:
                df.mask = mask
                m.append(df._mask)
            mask = m
            masksum = -1
        else:
            if hasattr(self,'_maskingObject'): # not masking object but self has one, delete it
                del self._maskingObject
            masksum = np.sum(mask)
            self.maskIndices = np.cumsum([np.sum(1-M) for M in mask])[:-1]
        if masksum==0:
            pass#warnings.warn('Provided mask has no masked elements!')
        elif masksum==self.I.size:
            warnings.warn('Provided mask masks all elements!')
        self._mask = mask
        for _,val in self.__dict__.items():
            if hasattr(val,'extractData'):
                val.mask = mask

        self.maskIndices = np.cumsum([np.sum(1-M) for M in mask])[:-1]

    @property
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')    

    def save(self, filename):
        try:                                # Opening the given file with an error catch
            fileObject = open(filename, 'wb')
        except IOError as e:                        # Catch all IO-errors
            print("Error in opening file:\n{}".format(e))
        else:
            pickle.dump(self, fileObject)
            fileObject.close()  


    def __eq__(self, other): 
        return np.logical_and(set(self.__dict__.keys()) == set(other.__dict__.keys()),self.__class__ == other.__class__)

    def __str__(self):
        string = '{} with settings:\n'.format(self.__class__)
        for attrib in self.settings:
            string+='{}:\t{}\n'.format(attrib,self.settings[attrib])
        for attrib in self.__dict__:
            if attrib=='_settings':
                continue
            string+='{}:\t{}\n'.format(attrib,self.__dict__[attrib])
        string+='\n'    
        return string

    def __getitem__(self,index):
        if len(self.convertedFiles) > 0:
            data = self.convertedFiles
        else:
            data = self.dataFiles
        try:
            return data[index]
        except IndexError:
            raise IndexError('Provided index {} is out of bounds for DataSet with length {}.'.format(index,len(self)))

    def __len__(self):
        if len(self.convertedFiles) > 0:
            data = self.convertedFiles
        else:
            data = self.dataFiles
        return len(data)
    
    def __iter__(self):
        self._index=0
        return self
    
    def __next__(self):
        if self._index >= len(self):
            raise StopIteration
        if len(self.convertedFiles) > 0:
            data = self.convertedFiles
        else:
            data = self.dataFiles
        result = data[self._index]
        self._index += 1
        return result

    def next(self):
        return self.__next__()
    
    def append(self,item):
        try:
            correctDataFiles = isListOfDataFiles(item)
            [self._dataFiles.append(file) for file in correctDataFiles if file.type=='hdf' or file.type=='MultiFLEXX' or file.type=='FlatCone']
            [self._convertedFiles.append(file) for file in correctDataFiles if file.type=='nxs']
        except Exception as e:
            raise(e)
        self._getData
    
    def __reversed__(self):
        return self[::-1]
    
    def __delitem__(self,index):
        if len(self.convertedFiles) > 0:
            if index < len(self.convertedFiles):
                del self.convertedFiles[index]
            else:
                raise IndexError('Provided index {} is out of bounds for DataSet with length {}.'.format(index,len(self.convertedFiles)))
        else:
            if index < len(self.dataFiles):
                del self.dataFiles[index]
            else:
                raise IndexError('Provided index {} is out of bounds for DataSet with length {}.'.format(index,len(self.dataFiles)))
        self._getData




    @_tools.KwargChecker()
    def convertDataFile(self,dataFiles=None,binning=None,saveLocation=None,saveFile=False):
        """Conversion method for converting scan file(s) to hkl file. Converts the given hdf file into NXsqom format and saves in a file with same name, but of type .nxs.
        Copies all of the old data file into the new to ensure complete redundancy. Determines the binning wanted from the file name of normalization file.

        Kwargs:

            - dataFiles (DataFile, string or list of): File path(s), file must be of hdf format (default self.dataFiles).

            - binning (int): Binning to be used when converting files (default 8).

            - saveLocation (string): File path to save location of data file(s) (defaults to same as raw file).

            - saveFile (bool): If true, the file(s) will be saved as nxs-files. Otherwise they will only persis in memory.

        Raises:

            - IOError

            - AttributeError
            
        """


        if dataFiles is None:
            if len(self.dataFiles)==0:
                raise AttributeError('No data files file provided either through input of in the DataSet object.')
        else:
            dataFiles = isListOfDataFiles(dataFiles)
        

        
        dataFiles = self.dataFiles
        convertedFiles = []
        for rawfile in dataFiles:
            convFile = rawfile.convert(binning)

            if saveFile: # TODO:
                if not saveLocation is None:
                    directory,file = os.path.split(saveLocation)
                    directory = os.path.abspath(directory)
                    if file == '':
                        file = os.path.split(rawfile.fileLocation)[1]
                    fileName = os.path.splitext(file)[0]
                    saveloc = os.path.join(directory,fileName+'.nxs')
                else:
                    saveloc = rawfile.fileLocation.replace('.hdf','.nxs')
                
                convFile.saveNXsqom(saveloc)
            
            #file.close()
            
            convertedFiles.append(convFile)
        self._convertedFiles = []
        self.convertedFiles = convertedFiles    
        self._getData()

        # check if any masking is needed due to defect detectors (If from CAMEA)
            
    def _getData(self): # Internal method to populate I,qx,qy,energy,Norm and Monitor
        
        if len(self.convertedFiles)!=0:
            self.I,self.qx,self.qy,self.energy,self.Norm,self.Monitor,self.a3,self.a3Off,self.a4,self.a4Off,self.instrumentCalibrationEf, \
            self.instrumentCalibrationA4,self.instrumentCalibrationEdges,self.Ei,self.scanParameters,\
            self.scanParameterValues,self.scanParameterUnits,self.h,self.k,self.l,self.mask = MJOLNIR.Data.DataFile.extractData(self.convertedFiles)
        else:
            self.I,self.Monitor,self.a3,self.a3Off,self.a4,self.a4Off,self.instrumentCalibrationEf, \
            self.instrumentCalibrationA4,self.instrumentCalibrationEdges,self.Ei,self.scanParameters,\
            self.scanParameterValues,self.scanParameterUnits,self.mask = MJOLNIR.Data.DataFile.extractData(self.dataFiles)

        if len(self.convertedFiles)!=0:
            self.sample = [d.sample for d in self]
        elif len(self.dataFiles)!=0:
            self.sample = [d.sample for d in self]

    @_tools.KwargChecker()
    def binData3D(self,dx,dy,dz,rlu=True,dataFiles=None):
        """Bin a converted data file into voxels with sizes dx*dy*dz. Wrapper for the binData3D functionality.

        Args:

            - dx (float): step sizes along the x direction (required).

            - dy (float): step sizes along the y direction (required).

            - dz (float): step sizes along the z direction (required).

        Kwargs:

            - rlu (bool): If True, the rotate QX,QY is used for binning (default True)

            - datafile (string or list of strings): Location(s) of data file to be binned (default converted file in DataSet).

        Raises:

            - AttributeError

        Returns:

            - Datalist: List of converted data files having 4 sub arrays: Intensity(counts), Monitor, Normalization, Normalization count

            - bins: 3 arrays containing edge positions in x, y, and z directions.
        """
        
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                I = self.I.extractData()
                qx = self.qx.extractData()
                qy = self.qy.extractData()
                energy = self.energy.extractData()
                Norm = self.Norm.extractData()
                Monitor = self.Monitor.extractData()
                samples = self.sample
                maskIndices = self.maskIndices

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor,samples,maskIndices = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData(),DS.sample,DS.maskIndices
        if rlu: # Rotate data
            Q = [[QX,QY] for QX,QY in zip(np.split(qx,maskIndices),np.split(qy,maskIndices))]
            qx,qy = np.concatenate([np.einsum('ij,j...->i...',s.RotMat,q) for s,q in zip(samples,Q)],axis=1)
        pos=[qx,qy,energy]
        returnData,bins = binData3D(dx=dx,dy=dy,dz=dz,pos=pos,data=I,norm=Norm,mon=Monitor)

        return returnData,bins

    @_tools.KwargChecker()
    def cut1D(self,q1,q2,width,minPixel,Emin,Emax,rlu=True,plotCoverage=False,extend=True,dataFiles=None,constantBins=False,positions=None,I=None,Norm=None,Monitor=None,ufit=False):
        """Wrapper for 1D cut through constant energy plane from q1 to q2 function returning binned intensity, monitor, normalization and normcount. The full width of the line is width while height is given by Emin and Emax. 
        the minimum step sizes is given by minPixel.
        
        .. note::
            Can only perform cuts for a constant energy plane of definable width.
        
        Args:
            
            - q1 (3D or 2D array): Start position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
            
            - q2 (3D or 2D array): End position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
            
            - width (float): Full width of cut in q-plane in 1/AA.
            
            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.
            
            - Emin (float): Minimal energy to include in cut.
            
            - Emax (float): Maximal energy to include in cut
            
        Kwargs:
            
            - rlu (bool): If True, coordinates given are interpreted as (h,k,l) otherwise as (qx,qy)

            - plotCoverage (bool): If True, generates plot of all points in the cutting plane and adds bounding box of cut (default False).

            - extend (bool): Whether or not the cut from q1 to q2 is to be extended throughout the data (default true)

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.

            - constantBins (bool): If True only bins of size minPixel is used (default False)

            - ufit (bool): If True a uFit Dataset object is returned in stead of pandas data frame
        
        
        Returns:
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int for 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
        """
        if np.all([positions is None,I is None,Norm is None,Monitor is None]):
            
            if dataFiles is None:
                if len(self.convertedFiles)==0:
                    raise AttributeError('No data file to be binned provided in either input or DataSet object.')
                else:
                    I = self.I.extractData()
                    qx = self.qx.extractData()
                    qy = self.qy.extractData()
                    energy = self.energy.extractData()
                    Norm = self.Norm.extractData()
                    Monitor = self.Monitor.extractData()
                    samples = self.sample
                    maskIndices = self.maskIndices

            else: 
                DS = DataSet(convertedFiles = dataFiles)
                I,qx,qy,energy,Norm,Monitor,samples,maskIndices = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData(),DS.sample,DS.maskIndices
                

            if rlu==True: # Recalculate H,K,L to qx
                q1,q2 = self.convertToQxQy([q1,q2])
                # Rotate all data files to fit with first data file
                #thetaDifference = [s.theta-samples[0].theta for s in samples]
                rotationMatrices = [np.dot(samples[0].RotMat.T,s.RotMat) for s in samples]#[_tools.Rot(theta,deg=False) for theta in thetaDifference]
                Q = [[QX,QY] for QX,QY in zip(np.split(qx,maskIndices),np.split(qy,maskIndices))]
                qx,qy = np.concatenate([np.einsum('ij,j...->i...',rot,q) for rot,q in zip(rotationMatrices,Q)],axis=1)

                positions = np.array([qx,qy,energy])
                
                
            else:
                positions = np.array([qx,qy,energy])
            
        if np.all(np.isclose(q1,q2)):
            if rlu:
                q1,q2 = self.convertToHKL([q1,q2])
            raise AttributeError('Provided Q points are equal. Got ({}) and ({}).'.format(', '.join([str(x) for x in q1]),', '.join([str(x) for x in q2])))

        Data,[binpositionsTotal,orthopos,EArray] = cut1D(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,
                                                                minPixel=minPixel,Emin=Emin,Emax=Emax,plotCoverage=plotCoverage,
                                                                extend=extend,constantBins=constantBins)

        if len(binpositionsTotal) == 0:
            return pd.DataFrame([],columns=['Qx','Qy','H','K','L','Energy','Intensity','Monitor','Normalization','BinCount','Int']),[binpositionsTotal,orthopos,EArray]    
        QxBin,QyBin = binpositionsTotal[:,:2].T

        if rlu:
                binpositionsTotal = np.concatenate([self.convertToHKL(binpositionsTotal[:,:2]),binpositionsTotal[:,-1].reshape(-1,1)],axis=1)
                orthopos = self.convertToHKL(orthopos)
                HBin,KBin,LBin,EnergyBin = binpositionsTotal.T
        else:
            HBin,KBin,LBin,EnergyBin = np.concatenate([self.convertToHKL(binpositionsTotal[:,:2]),binpositionsTotal[:,-1].reshape(-1,1)],axis=1).T

        def meaning(X):
            return 0.5*(X[:-1]+X[1:])

        Qx = meaning(QxBin)
        Qy = meaning(QyBin)
        H = meaning(HBin)
        K = meaning(KBin)
        L = meaning(LBin)
        Energy = meaning(EnergyBin)


        I,Mon,Norm,BinC = Data
        DataValues = [Qx,Qy,H,K,L,Energy,I,Mon,Norm,BinC]
        columns = ['Qx','Qy','H','K','L','Energy','Intensity','Monitor','Normalization','BinCount']
        dtypes = [np.float]*6+[np.int]*2+[np.float]+[np.int]

        pdData = pd.DataFrame()
        if not len(I) == 0:
            for dat,col,typ in zip(DataValues,columns,dtypes):
                pdData[col] = dat.astype(typ)
            pdData['Int'] = pdData['Intensity']*pdData['BinCount']/(pdData['Normalization']*pdData['Monitor'])

        
        if not ufit:
            return pdData,[binpositionsTotal,orthopos,EArray]
        
        if rlu:
            q1,q2 = self.convertToHKL([q1,q2])
        ufitData = self.generateUFitDataset(pdData,q1,q2,rlu,width=width,Emin=Emin,Emax=Emax,minPixel=minPixel)
        
        return ufitData

        
    
    @_tools.KwargChecker(function=plt.errorbar,include=np.concatenate([_tools.MPLKwargs,['ticks','tickRound','mfc','markeredgewidth','markersize']])) #Advanced KWargs checker for figures
    def plotCut1D(self,q1,q2,width,minPixel,Emin,Emax,rlu=True,ax=None,plotCoverage=False,extend=True,Data=None,dataFiles=None,constantBins=False,ufit=False,outputFunction=print,**kwargs):  
        """plot new or already performed cut.
        
        Args:
                
                - q1 (3D or 2D array): Start position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
                
                - q2 (3D or 2D array): End position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
                
                - width (float): Full width of cut in q-plane in 1/AA, only needed for new cut (default None)
                
                - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.

                - Emin (float): Minimal energy to include in cut.
                
                - Emax (float): Maximal energy to include in cut.
                
            Kwargs:
                
                - Data (pandas Dataframe): Data, if previously made cut is to be plotted (default None)
                
                - rlu (bool): If True, coordinates given are interpreted as (h,k,l) otherwise as (qx,qy)
        
                - extend (bool): Whether or not the cut from q1 to q2 is to be extended throughout the data (default true)
                
                - plotCoverage (bool): If True, generates plot of all points in the cutting plane and adds bounding box of cut (default False).
        
                - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.
        
                - constantBins (bool): If True only bins of size minPixel is used (default False)

                - outputFunction (function): Function called on output string (default print)
        
        """
        if Data is None:
            Data, bins = self.cut1D(q1=q1,q2=q2,width=width,minPixel=minPixel,Emin=Emin,Emax=Emax,rlu=rlu,plotCoverage=plotCoverage,constantBins=constantBins)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            INT = np.divide(Data['Intensity']*Data['BinCount'],Data['Monitor']*Data['Normalization'])
            INT_err = np.divide(np.sqrt(Data['Intensity'])*Data['BinCount'],Data['Monitor']*Data['Normalization'])
        
        
        
        num=len(Data)
        q1 = np.array(q1,dtype=float)
        q2 = np.array(q2,dtype=float)
        if rlu:
            variables = ['H','K','L']
        else:
            variables = ['Qx','Qy']
    
        variables = variables+['Energy']
        if not 'ticks' in kwargs:
            ticks = 5
        else:
            ticks = kwargs['ticks']
            kwargs = _tools.without_keys(dictionary=kwargs, keys='ticks')
        
        if not 'tickRound' in kwargs:
            tickRound = 3
        else:
            tickRound = kwargs['tickRound']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='tickRound')
            
        if not 'fmt' in kwargs:
            kwargs['fmt'] = '.'
        
        xvalues = np.round(np.linspace(0,num-1,ticks)).astype(int)
        my_xticks=[]
        for i in xvalues:
            my_xticks.append('\n'.join(map(lambda x:('{:.'+str(tickRound)+'f}').format(x),[np.round(Data[var][i],tickRound) for var in variables])))
        
        
        Data['binDistance'] = np.linalg.norm(Data[variables]-np.array(Data[variables].iloc[1]),axis=1)
        
        if ax is None:
            plt.figure()
            ax = plt.gca()
        
        if not 'label' in kwargs:
            kwargs['label'] = 'Data'
            
        ax.errorbar(Data['binDistance'],INT,yerr=INT_err,**kwargs)
        
        ax.set_xticks(Data['binDistance'].iloc[xvalues])
        ax.set_xticklabels(my_xticks, multialignment="center",ha="center")
        
        def calculateIndex(binDistance,x):
            idx = np.argmin(np.abs(binDistance-x))
            return idx
        
        ax.calculateIndex = lambda x: calculateIndex(Data['binDistance'],x)
        
        if rlu==False:
            ax.set_xlabel(r'$Q_x [\AA^{-1}]$'+'\n'+r'$Q_y [\AA^{-1}]$'+'\n'+r'E [meV]')
            def format_coord(x,y,ax,binCenter):# pragma: no cover
                index = ax.calculateIndex(x)
                qx,qy,E = binCenter[index]
                return  "qx = {0:.3e}, qy = {1:.3e}, E = {2:.3f}, I = {3:0.4e}".format(qx,qy,E,y)
        else:
            def format_coord(x,y,ax,binCenter):# pragma: no cover
                index = ax.calculateIndex(x)
                h,k,l,E = binCenter[index]
                return  "H = {0:.3e}, K = {1:.3e}, L = {2:.3e}, E = {3:.3f}, I = {4:0.4e}".format(h,k,l,E,y)
            ax.set_xlabel('$Q_h$ [RLU]\n$Q_k$ [RLU]\n$Q_l$ [RLU]\nE [meV]')
        
        
        def onclick(event,ax,Data,outputFunction):# pragma: no cover
            if ax.in_axes(event):
                try:
                    C = ax.get_figure().canvas.cursor().shape() # Only works for pyQt5 backend
                except:
                    pass
                else:
                    if C != 0: # Cursor corresponds to arrow
                        return
        
                x = event.xdata
                y = event.ydata
                printString = ax.format_coord(x,y)
                index = ax.calculateIndex(x)
            
                cts = int(Data['Intensity'][index])
                Mon = int(Data['Monitor'][0])
                Norm = float(Data['Normalization'][0])
                NC = int(Data['BinCount'][0])
                printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
                outputFunction(printString)
        
        
        ax.xaxis.set_label_coords(1.15, -0.025)
        ax.set_ylabel('$I$ [arb.u.]')
        plt.tight_layout()
        
        
        ax.format_coord = lambda x,y: format_coord(x,y,ax,np.array(Data[variables]))
        ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,Data,outputFunction=outputFunction))
        
        if ufit==True:
            ufitdata = self.generateUFitDataset(pdData=Data,q1=q1,q2=q2,rlu=rlu,width=width,Emin=Emin,Emax=Emax,minPixel=minPixel)
            return ax,ufitdata
        
        return ax,Data,bins


    @_tools.KwargChecker()
    def cutQE(self,q1,q2,width,minPixel,EnergyBins,rlu=True,extend=True,dataFiles=None,constantBins=False):
        """Wrapper for cut data into maps of q and intensity between two q points and given energies. This is performed by doing consecutive constant energy planes.

        Args:

            - q1 (3D or 2D array): Start position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
            
            - q2 (3D or 2D array): End position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
            
            - width (float): Full width of cut in q-plane.
            
            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.

            - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

        Kwargs:

            - rlu (bool): If True, coordinates given are interpreted as (h,k,l) otherwise as (qx,qy)

            - extend (bool): If True, cut is extended to edge of measured area instead of only between provided points.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.
    
            - constantBins (bool): If True only bins of size minPixel is used (default False)


        Returns:
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for 2D cut.
            
            - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (n * 3D arrays): n instances of center positions for the bins.

            - binDistance (n arrays): n instances of arrays holding the distance in q to q1.


        """
    
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                I = self.I.extractData()
                qx = self.qx.extractData()
                qy = self.qy.extractData()
                energy = self.energy.extractData()
                Norm = self.Norm.extractData()
                Monitor = self.Monitor.extractData()
                samples = self.sample
                maskIndices = self.maskIndices
    
        else: 
            #dataFiles = isListOfDataFiles(dataFiles)
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor,samples,maskIndices = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData(),DS.sample,DS.maskIndices
        
        
        if rlu==True: # Recalculate H,K,L to qx
            q1,q2 = self.convertToQxQy([q1,q2])
            # Rotate all data files to fit with first data file
            #thetaDifference = [s.theta-samples[0].theta for s in samples]
            rotationMatrices = [np.dot(samples[0].RotMat.T,s.RotMat) for s in samples]#[_tools.Rot(theta,deg=False) for theta in thetaDifference]
            Q = [[QX,QY] for QX,QY in zip(np.split(qx,maskIndices),np.split(qy,maskIndices))]
            qx,qy = np.concatenate([np.einsum('ij,j...->i...',rot,q) for rot,q in zip(rotationMatrices,Q)],axis=1)
            positions = np.array([qx,qy,energy])
                
        else:
            positions = np.array([qx,qy,energy])

        intensityArray = []
        monitorArray = []
        normalizationArray = []
        normcountArray = []
        centerPos = []
        returnpositions = []
        binDistance = []

        dirvec = (np.array(q2) - np.array(q1)).astype(float)
        dirvec /= np.linalg.norm(dirvec)

        dataFrame = []
        for i in np.arange(len(EnergyBins)-1):

            _local,position = self.cut1D(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,
                                    width=width,minPixel=minPixel,Emin=EnergyBins[i],Emax=EnergyBins[i+1],
                                    plotCoverage=False,extend=extend,constantBins=constantBins,dataFiles=dataFiles,rlu=False)                                      

            _local['energyCut'] = i
            dataFrame.append(_local)

            if len(_local)==0:
                continue
            returnpositions.append(position)


            thisCenterPos = 0.5*(position[0][:-1]+position[0][1:])
            centerPos.append(thisCenterPos)
            thisBinDistance = np.dot(thisCenterPos[:,:len(q1)] - q1, dirvec)
            binDistance.append(thisBinDistance)
        if len(dataFrame)>1:
            dataFrame = pd.concat(dataFrame)
        elif len(dataFrame) == 1:
            dataFrame = dataFrame[0]
        

        return dataFrame,returnpositions,centerPos,binDistance

 
    
    @_tools.KwargChecker(function=plt.errorbar,include=_tools.MPLKwargs)
    def plotCutQE(self,q1,q2,width,minPixel,EnergyBins,rlu=True,ax=None,dataFiles=None,constantBins=False,**kwargs): 
        """Plotting wrapper for the cutQE method. Generates a 2D intensity map with the data cut by cutQE. 
    
        .. warning::
           Deprecated! Instead use the plotCutQELine tool with only two q points

        .. note::
            Positions shown in tool tip reflect the closes bin center and are thus limited to the area where data is present.
        
        Args:

            - q1 (3D or 2D array): Start position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
            
            - q2 (3D or 2D array): End position of cut in format (h,k,l) or (qx,qy) depending on rlu flag.
            
            - width (float): Full width of cut in q-plane.
            
            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.

            - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

        Kwargs:

            - rlu (bool): If True, coordinates given are interpreted as (h,k,l) otherwise as (qx,qy)
            
            - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.

            - constantBins (bool): If True only bins of size minPixel is used (default False)
        
            - kwargs: All other keywords will be passed on to the ax.errorbar method.
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for 1D cut.
            
            - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (n * 3D arrays): n instances of center positions for the bins.

            - binDistance (n arrays): n instances of arrays holding the distance in q to q1.
        """
        
        
        #if dataFiles is None:
        #    if len(self.convertedFiles)==0:
        #        raise AttributeError('No data file to be binned provided in either input or DataSet object.')
        #    else:
        #        I = self.I.extractData()
        #        qx = self.qx.extractData()
        #        qy = self.qy.extractData()
        #        energy = self.energy.extractData()
        #        Norm = self.Norm.extractData()
        #        Monitor = self.Monitor.extractData()
        #        samples = self.sample
        #        maskIndices = self.maskIndices

        #else: 
        #    DS = DataSet(convertedFiles = dataFiles)
        #    I,qx,qy,energy,Norm,Monitor,samples,maskIndices = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData(),DS.sample,DS.maskIndices
            
        #if rlu==True: # Recalculate H,K,L to qx
        #    q1,q2 = self.convertToQxQy([q1,q2])
        #    # Rotate all data files to fit with first data file
        #    thetaDifference = [s.theta-samples[0].theta for s in samples]
        #    rotationMatrices = [np.dot(samples[0].RotMat.T,s.RotMat) for s in samples]#[_tools.Rot(theta,deg=False) for theta in thetaDifference]
        #    Q = [[QX,QY] for QX,QY in zip(np.split(qx,maskIndices),np.split(qy,maskIndices))]
        #    qx,qy = np.concatenate([np.einsum('ij,j...->i...',rot,q) for rot,q in zip(rotationMatrices,Q)],axis=1)


        #positions = np.array([qx,qy,energy])
        #return plotCutQELine(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,
        #                minPixel=minPixel,EnergyBins=EnergyBins,rlu=rlu,ax = ax,constantBins=constantBins,**kwargs)
        return self.plotCutQELine(QPoints=[q1,q2],width=width,minPixel=minPixel,EnergyBins=EnergyBins,rlu=rlu,ax=ax,dataFiles=dataFiles,constantBins=constantBins,**kwargs)
    
    @_tools.KwargChecker()
    def cutPowder(self,EBinEdges,qMinBin=0.01,dataFiles=None,constantBins=False):
        """Cut data powder map with intensity as function of the length of q and energy. 

        Args:
            
            - EBinEdges (list): Bin edges between which the cuts are performed.

        Kwargs:

            - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.

            - constantBins (bool): If True only bins of size minPixel is used (default False)


        Returns:
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for powder cut.
            
            - qbins (n arrays): n arrays holding the bin edges along the length of q

        """
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                I = self.I.extractData()
                qx = self.qx.extractData()
                qy = self.qy.extractData()
                energy = self.energy.extractData()
                Norm = self.Norm.extractData()
                Monitor = self.Monitor.extractData()

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData()
        
        positions = np.array([qx,qy,energy])

        return cutPowder(positions=positions,I=I,Norm=Norm,Monitor=Monitor,
                        EBinEdges=EBinEdges,qMinBin=qMinBin,constantBins=constantBins)

    
    @_tools.KwargChecker(function=plt.pcolormesh,include=np.concatenate([_tools.MPLKwargs,['vmin','vmax','edgecolors']]))
    def plotCutPowder(self,EBinEdges,qMinBin=0.01,ax=None,dataFiles=None,constantBins=False,log=False,colorbar=True,outputFunction=print,**kwargs):
        """Plotting wrapper for the cutPowder method. Generates a 2D plot of powder map with intensity as function of the length of q and energy.  
        
        .. note::
            Can only perform cuts for a constant energy plane of definable width.
        
        Args:

            - EBinEdges (list): Bin edges between which the cuts are performed.

        Kwargs:
            
            - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.
            
            - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
            
            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.

            - constantBins (bool): If True only bins of size minPixel is used (default False)

            - log (bool): If true, logarithm to intensity is plotted (default False)

            - colorbar (bool): If True a colorbar is added to the figure (default True)

            - outputFunction (function): Function called on output strung (default print)

            - kwargs: All other keywords will be passed on to the ax.pcolormesh method.
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for poweder cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).

        """

        Data,qbins = self.cutPowder(EBinEdges=EBinEdges,qMinBin=qMinBin,dataFiles=dataFiles,constantBins=constantBins)
        
        if ax is None:
            plt.figure()
            ax = plt.gca()
        pmeshs = []
        
        eMean = 0.5*(EBinEdges[:-1]+EBinEdges[1:])
        for i,dat in Data[['Int','EnergyCut']].groupby('EnergyCut'):
            if len(dat)==0:
                continue
            if log:
                ints = np.log10(dat['Int']+1e-20).values.reshape((len(qbins[i])-1,1)).T
                
            else:
                ints = dat['Int'].values.reshape((len(qbins[i])-1,1)).T
            pmeshs.append(ax.pcolormesh(qbins[i],[EBinEdges[i],EBinEdges[i+1]],ints,**kwargs))
        
        def calculateIndex(x,y,eMean,qBin,EBinEdges):# pragma: no cover
            if y>EBinEdges[-1] or y<EBinEdges[0]:
                return -1,-1
            EIndex = np.argmin(np.abs(y-eMean))
            if x>qBin[EIndex][-1] or x<qBin[EIndex][0]:
                return EIndex,-1
            qIndex = np.argmin(np.abs(x-0.5*(qBin[EIndex][:-1]+qBin[EIndex][1:])))
            return EIndex,qIndex
        ax.calculateIndex = lambda x,y: calculateIndex(x,y,eMean,qbins,EBinEdges)

        def format_coord(x,y,ax,dat,qbins):# pragma: no cover
            EIndex,qIndex = ax.calculateIndex(x,y)
            if EIndex == -1 or qIndex == -1:
                return "outside range"
            Intensity = dat[dat['EnergyCut']==EIndex]['Int'][qIndex]
            return  "|q| = {0:.3f}, E = {1:.3f}, I = {2:0.4e}".format(qbins[EIndex][qIndex],eMean[EIndex],Intensity)
            
        ax.format_coord = lambda x,y: format_coord(x,y,ax,Data[['Int','EnergyCut']],qbins)
        ax.set_xlabel(r'$|q| [\AA^{-1}]$')
        ax.set_ylabel('E [meV]')
        
        def set_clim(VMin,VMax,pmesh):
            for pm in pmeshs:
                pm.set_clim(VMin,VMax)

        ax.set_clim = lambda VMin,VMax: set_clim(VMin,VMax,pmeshs)
        
        def onclick(event,ax,dat,outputFunction):# pragma: no cover
            if ax.in_axes(event):
                try: 
                    c = ax.get_figure().canvas.cursor().shape()
                except:
                    pass
                else:
                    if c!=0:
                        return
                x = event.xdata
                y = event.ydata
                printString = ax.format_coord(x,y)
                EIndex,index = ax.calculateIndex(x,y)
                if EIndex == -1 or index == -1:
                    return
                _local = dat[dat['EnergyCut']==EIndex]
                cts = _local['Intensity'][index]
                Mon = _local['Monitor'][index]
                Norm = _local['Normalization'][index]
                NC = _local['BinCount'][index]
                printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
                outputFunction(printString)

        if not 'vmin' in kwargs or not 'vmax' in kwargs:
            if log:
                minVal = np.log10(Data['Int'].min()+1e-20)
                maxVal = np.log10(Data['Int'].max()+1e-20)
            else:
                minVal = Data['Int'].min()
                maxVal = Data['Int'].max()
            
            ax.set_clim(minVal,maxVal)
        ax.pmeshs = pmeshs
        ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,Data,outputFunction=outputFunction))
        if colorbar:
            ax.colorbar = ax.get_figure().colorbar(ax.pmeshs[0],pad=0.1)

        return ax,Data,qbins

    @_tools.KwargChecker()
    @_tools.overWritingFunctionDecorator(RLUAxes.createRLUAxes)
    def createRLUAxes(*args,**kwargs): # pragma: no cover
        raise RuntimeError('This code is not meant to be run but rather is to be overwritten by decorator. Something is wrong!! Should run {}'.format(RLUAxes.createRLUAxes))

    @_tools.KwargChecker()
    @_tools.overWritingFunctionDecorator(RLUAxes.createQEAxes)
    def createQEAxes(*args,**kwargs): # pragma: no cover
        raise RuntimeError('This code is not meant to be run but rather is to be overwritten by decorator. Something is wrong!! Should run {}'.format(RLUAxes.createQEAxes))
    
    #@_tools.KwargChecker(function=plt.pcolormesh,include=['vmin','vmax','colorbar','zorder'])
    def plotQPlane(self,EMin=None,EMax=None,EBins=None,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=False,ax=None,rlu=True,dataFiles=None,xScale=1.0,yScale=1.0,outputFunction=print,**kwargs):
        """Wrapper for plotting tool to show binned intensities in the Q plane between provided energies.
            
        Kwargs:
            
            - EMin (float): Lower energy limit (Default None).
            
            - EMax (float): Upper energy limit (Default None).

            - EBins (list): List of energy bins (Default None).

            - binning (str): Binning scheme, either 'xy' or 'polar' (default 'xy').
            
            - xBinTolerance (float): bin sizes along x direction (default 0.05). If enlargen is true, this is the minimum bin size.

            - yBinTolerance (float): bin sizes along y direction (default 0.05). If enlargen is true, this is the minimum bin size.
            
            - enlargen (bool): If the bin sizes should be adaptive (default False). If set true, bin tolerances are used as minimum bin sizes.

            - log (bool): Plot intensities as the logarithm (default False).
            
            - ax (matplotlib axes): Axes in which the data is plotted (default None). If None, the function creates a new axes object.

            - rlu (bool): If true and axis is None, a new reciprocal lattice axis is created and used for plotting (default True).

            - dataFiles (DataFile): If set, method uses these converted data files instead of the ones in self (default None)

            - vmin (float): Lower limit for colorbar (default min(Intensity)).
            
            - vmax (float): Upper limit for colorbar (default max(Intensity)).

            - colorbar (bool): If True, a colorbar is created in figure (default False)

            - zorder (int): If provided decides the z ordering of plot (default 10)

            - outputFunction (function): Function called on output string (default print)
            
            - other: Other key word arguments are passed to the pcolormesh plotting algorithm.
            
        Returns:
            
            - dataList (list): List of all data points in format [Intensity, Monitor, Normalization, Normcount]

            - bins (list): List of bin edges as function of plane in format [xBins,yBins].

            - ax (matplotlib axes): Returns provided matplotlib axis
            
        .. note::
            The axes object has a new method denoted 'set_clim' taking two parameters (VMin and VMax) used to change axes colouring.
            
        .. note::
            If a 3D matplotlib axis is provided, the planes are plotted in 3D with the provided energy bins. As the method 
            contourf is used and it needs X,Y, and Z to have same shape, x and y are found as middle of bins. 
            
        """
        

        if EMax is None or EMin is None:
            if EBins is None:
                raise AttributeError('Either minimal/maximal energy or the energy bins is to be given.')
            else:
                if len(EBins)<=1:
                    raise AttributeError('Length of provided energy bins is {}, while at least 2 is needed! Received "{}"'.format(len(EBins),EBins))
                EBins = np.asarray(EBins)
        else:
            if EMin>=EMax:
                raise AttributeError('Provided limits are either wrong or the same. Received EMin={} and EMax={}, expects EMin<EMax.'.format(EMin,EMax))
            EBins = np.array([EMin,EMax])

        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                I = self.I.extractData()#
                qx = self.qx.extractData()#
                qy = self.qy.extractData()#
                energy = self.energy.extractData()#
                Norm = self.Norm.extractData()#
                Monitor = self.Monitor.extractData()#
                samples = self.sample
                maskIndices = self.maskIndices

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor,samples,maskIndices = DS.I,DS.qx,DS.qy,DS.energy,DS.Norm,DS.Monitor,DS.sample,DS.maskIndices
        if ax is None:
            if rlu is True:
                ax = self.createRLUAxes()
            else:
                fig,ax = plt.subplots()

            _3D = False
        else:
            if ax.name =='3d':
                _3D = True
            else:
                _3D = False
        
        #if len(qx.shape)==1 or len(qx.shape)==4:
        #    
        #    qx = np.concatenate(qx,axis=0)
        #    qy = np.concatenate(qy,axis=0)
        #    energy = np.concatenate(energy,axis=0)
        #    I = np.concatenate(I,axis=0)
        #    Norm = np.concatenate(Norm,axis=0)
        #    Monitor = np.concatenate(Monitor,axis=0)
        
        
        if rlu == True: # Rotate positions with taslib.misalignment to line up with RLU
            Q = [[QX,QY] for QX,QY in zip(np.split(qx,maskIndices),np.split(qy,maskIndices))]
            qx,qy = np.concatenate([np.einsum('ij,j...->i...',s.RotMat,q) for s,q in zip(samples,Q)],axis=1)
            
        if 'zorder' in kwargs:
            zorder = kwargs['zorder']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='zorder')
        else:
            zorder = 10

        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='cmap')
        else:
            cmap = None

        intensity = []
        monitorCount = []
        Normalization = []
        NormCount = []
        Int = []
        xBins = []
        yBins = []
        offset = [] # Only used for binning in polar
        pmeshs = []

        binnings = ['xy','polar']
        if not binning in binnings:
            raise AttributeError('The provided binning is not understood, should be {}'.format(', '.join(binnings)))

        for i in range(len(EBins)-1):
            #print('Binning {} to {}.'.format(EBins[i],EBins[i+1]))
            EBinEdges = [EBins[i],EBins[i+1]]
            e_inside = np.logical_and(energy>EBinEdges[0],energy<=EBinEdges[1])
            if np.sum(e_inside)==0:
                continue
            if binning == 'polar':
                
                x = np.arctan2(qy[e_inside],qx[e_inside]) # Gives values between -pi and pi
                bins = 20
                # Following block checks if measured area corresponds to alpha ~pi as arctan2 only gives
                # values back in range -pi to pi.
                if np.max(x.flatten())+xBinTolerance>np.pi and np.min(x.flatten())-xBinTolerance<-np.pi:
                    h = np.histogram(x.flatten(),bins = bins)
                    while np.max(h[0]==0) == False:
                        bins *= 2
                        h = np.histogram(x.flatten(),bins = bins)
                        if bins > 200:
                            break
                    if bins > 200: # If everything has been covered, do nothing.
                        offset.append(0.0)
                    else:
                        offset.append(2*np.pi-h[1][np.argmax(h[0]==0)]) # Move highest value of lump to fit 2pi
                        x = np.mod(x+offset[-1],2*np.pi)-np.pi # moves part above 2pi to lower than 2pi and make data fit in range -pi,pi
                        offset[-1]-=np.pi # As x is moved by pi, so should the offset
                else:
                    offset.append(0.0)

                y = np.linalg.norm([qx[e_inside],qy[e_inside]],axis=0)  
                if not enlargen:
                    xBins.append(np.arange(-np.pi,np.pi+xBinTolerance*0.999,xBinTolerance)) # Add tolerance as to ensure full coverage of parameter
                    yBins.append(np.arange(0,np.max(y)+yBinTolerance*0.999,yBinTolerance)) # Add tolerance as to ensure full coverage of parameter
                else:
                    xBins.append(_tools.binEdges(x,tolrance=xBinTolerance))
                    yBins.append(_tools.binEdges(y,tolerance=yBinTolerance))

            elif binning == 'xy':
                x = qx[e_inside]
                y = qy[e_inside]
                if not enlargen:
                    xBins.append(np.arange(np.min(x),np.max(x)+0.999*xBinTolerance,xBinTolerance)) # Add tolerance as to ensure full coverage of parameter
                    yBins.append(np.arange(np.min(y),np.max(y)+0.999*yBinTolerance,yBinTolerance)) # Add tolerance as to ensure full coverage of parameter
                else:
                    xBins.append(_tools.binEdges(x,tolerance=xBinTolerance))
                    yBins.append(_tools.binEdges(y,tolerance=yBinTolerance))
            
            X = x.flatten()
            Y = y.flatten()
            
            intensity.append(np.histogram2d(X,Y,bins=(xBins[i],yBins[i]),weights=I[e_inside])[0].astype(I.dtype))
            monitorCount.append(np.histogram2d(X,Y,bins=(xBins[i],yBins[i]),weights=Monitor[e_inside])[0].astype(Monitor.dtype))
            Normalization.append(np.histogram2d(X,Y,bins=(xBins[i],yBins[i]),weights=Norm[e_inside])[0].astype(Norm.dtype))
            NormCount.append(np.histogram2d(X,Y,bins=(xBins[i],yBins[i]),weights=np.ones_like(I[e_inside]))[0].astype(I.dtype))
                
                

            warnings.simplefilter('ignore')
            Int.append(np.divide(intensity[-1]*NormCount[-1],monitorCount[-1]*Normalization[-1]))
            warnings.simplefilter('once')

        if binning == 'polar':
            Qx = [np.outer(np.cos(xBins[i]-offset[i]),yBins[i]) for i in range(len(intensity))]
            Qy = [np.outer(np.sin(xBins[i]-offset[i]),yBins[i]) for i in range(len(intensity))]

        elif binning == 'xy':
            Qx =[np.outer(xBins[i],np.ones_like(yBins[i])) for i in range(len(intensity))]
            Qy =[np.outer(np.ones_like(xBins[i]),yBins[i]) for i in range(len(intensity))]
            
        
        if 'vmin' in kwargs:
            vmin = kwargs['vmin']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='vmin')
        else:
            vmin = np.min([np.nanmin(intens) for intens in Int])

        if 'vmax' in kwargs:
            vmax = kwargs['vmax']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='vmax')
        else:
            vmax = np.max([np.nanmax(intens) for intens in Int])

        if 'colorbar' in kwargs:
            colorbar = kwargs['colorbar']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='colorbar')
        else:
            colorbar = False
        pmeshs = []
        if log:
            Int = [np.log10(1e-20+np.array(intens)) for intens in Int]

        for i in range(len(EBins)-1):
            if _3D:
                QX = 0.25*np.array(np.array(Qx[i])[1:,1:]+np.array(Qx[i])[:-1,1:]+np.array(Qx[i])[1:,:-1]+np.array(Qx[i])[:-1,:-1])/xScale
                QY = 0.25*np.array(np.array(Qy[i])[1:,1:]+np.array(Qy[i])[:-1,1:]+np.array(Qy[i])[1:,:-1]+np.array(Qy[i])[:-1,:-1])/yScale
                #QY = np.array(np.array(Qy[i])[1:,1:])
                I = np.array(Int[i])
                levels = np.linspace(vmin,vmax,50)
                pmeshs.append(ax.contourf3D(QX,QY,I,zdir = 'z',offset=np.mean(EBins[i:i+2]),levels=levels,cmap=cmap,**kwargs))
            else:
                pmeshs.append(ax.pcolormesh(Qx[i],Qy[i],Int[i],zorder=zorder,cmap=cmap,**kwargs))
        if not _3D:
            ax.set_aspect('equal')
            ax.grid(True, zorder=0)
        def set_clim(pmeshs,vmin,vmax):
            for pmesh in pmeshs:
                pmesh.set_clim(vmin,vmax)
        

        if 'pmeshs' in ax.__dict__:
            ax.pmeshs = np.concatenate([ax.pmeshs,np.asarray(pmeshs)],axis=0)
        else:
            ax.pmeshs = pmeshs

        ax.set_clim = lambda vMin,vMax: set_clim(ax.pmeshs,vMin,vMax)

        if colorbar:
            ax.colorbar = ax.get_figure().colorbar(ax.pmeshs[0],pad=0.1)
            ax.colorbar.set_label('$I$ [arb.u.]', rotation=270)

        ax.set_clim(vmin,vmax)
        if _3D:
            minEBins = np.min(EBins)
            maxEBins = np.max(EBins)
            if not np.isclose(minEBins,maxEBins):
                ax.set_zlim(minEBins,maxEBins)
            else:
                ax.set_zlim(minEBins-0.1,maxEBins+0.1)
        else:
            def onclick(ax, event,Qx,Qy,data,outputFunction): # pragma: no cover
                if event.xdata is not None and ax.in_axes(event):
                    try:
                        C = ax.get_figure().canvas.cursor().shape() # Only works for pyQt5 backend
                    except:
                        pass
                    else:
                        if C != 0:
                            return
                if not ax.in_axes(event):
                    return
                printString = ''
                printString+=ax.format_coord(event.xdata, event.ydata)+', '
                QX = np.array(Qx[0])
                QY = np.array(Qy[0])
                
                Qx = 0.25*(QX[1:,1:]+QX[:-1,:-1]+QX[1:,:-1]+QX[:-1,1:])
                Qy = 0.25*(QY[1:,1:]+QY[:-1,:-1]+QY[1:,:-1]+QY[:-1,1:])

                arg = np.argmin(np.linalg.norm(np.array([Qx,Qy])-np.array([event.xdata,event.ydata]).reshape(2,1,1),axis=0))
                arg2D = np.unravel_index(arg,Qx.shape)
                
                cts = data[0][0][arg2D[0],arg2D[1]]
                Norm = data[1][0][arg2D[0],arg2D[1]]
                Mon = data[2][0][arg2D[0],arg2D[1]]
                NC = data[3][0][arg2D[0],arg2D[1]]
                warnings.simplefilter('ignore')
                Intensity = np.divide(cts*NC,Norm*Mon)
                warnings.simplefilter('once')
                localSize = np.linalg.norm(np.array([QX[arg2D[0],arg2D[1]]-Qx[arg2D[0],arg2D[1]],QY[arg2D[0],arg2D[1]]-Qy[arg2D[0],arg2D[1]]]))

                if not np.isfinite(Intensity) or np.linalg.norm(np.array([Qx[arg2D[0],arg2D[1]]-event.xdata,Qy[arg2D[0],arg2D[1]]-event.ydata]))>localSize:
                    printString+='I = NaN'
                else:
                    printString+='I = {:.4E}'.format(Intensity)
                    printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)

                outputFunction(printString)
            ax.cid = ax.figure.canvas.mpl_connect('button_press_event', lambda x: onclick(ax,x,Qx,Qy,[intensity,monitorCount,Normalization,NormCount],outputFunction=outputFunction))
            
        if len(Qx)!=0:
            xmin = np.min([np.min(qx) for qx in Qx])
            xmax = np.max([np.max(qx) for qx in Qx])
            ax.set_xlim(xmin,xmax)#np.min(Qx),np.max(Qx))
        
        if len(Qy)!=0:
            ymin = np.min([np.min(qy) for qy in Qy])
            ymax = np.max([np.max(qy) for qy in Qy])
            ax.set_ylim(ymin,ymax)#np.min(Qy),np.max(Qy))
        return [intensity,monitorCount,Normalization,NormCount],[Qx,Qy],ax

    @_tools.KwargChecker()
    def plotA3A4(self,dataFiles=None,ax=None,planes=[],log=False,returnPatches=False,binningDecimals=3,singleFigure=False,plotTessellation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
        """Plot data files together with pixels created around each point in A3-A4 space. Data is binned in the specified planes through their A3 and A4 values. 
        This can result in distorted binning when binning across large energy regions. Data is plotted using the pixels calculated for average plane value, i.e. 
        binning 7,8,9,10, and 11 patches for plane 9 are used for plotting.

        Kwargs:
            - dataFiles (DataFiles): single file or list of files to be binned together (Default self.convertedFiles)

            - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)

            - planes (list (of lists)): Planes to be plotted and binned (default [])

            - log (bool): Whether or not to plot intensities as logarithm (default False)

            - returnPatches (bool): If true the method returns the patches otherwise plotted in the given axis (default False)

            - binningDecimals (int): Number of decimal places Q positions are rounded before binning (default 3)

            - singleFigure (bool): If true, all planes are plotted in same figure (default False)

            - plotTessellation (bool): Plot Tessellation of points (default False)

            - Ei_err (float): Tolerance of E_i for which the values are equal (default = 0.05)

            - temperature_err (float): Tolerance of temperature for which the values are equal (default = 0.2)
            
            - magneticField_err (float): Tolerance of magnetic field for which the values are equal (default = 0.2)
            
            - electricField_err (float): Tolerance of electric field for which the values are equal (default = 0.2)

        Returns:
            
            - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.

        Raises:

            - NotImplimentedError

        Examples:

        The following example will combine the two files and plot all of the available planes in different figures.

        >>> DS = DataSet.DataSet(convertedFiles=[--.nxs,---.nxs])
        >>> plt.figure()
        >>> ax = plt.gca()
        >>>
        >>> DataSet.plotA3A4(DS.convertedFiles,ax=ax)

        If only a subset of planes or different planes are to be combined the following will achieve this:

        >>> DataSet.plotA3A4(DS.convertedFiles,ax=ax,planes=[0,1,2,3,[4,5,6],[8,9]])

        Here planes 0 through 3 are plotted separately while 4,5, and 6 as well as 8 and 9 are binned.

        .. note::
            Binning planes from different analysers might result in nonsensible binnings.

        """
        if dataFiles is None: 
            dataFiles = self.convertedFiles
        
        return plotA3A4(dataFiles,ax=ax,planes=planes,log=log, returnPatches=returnPatches,binningDecimals=binningDecimals,
        singleFigure=singleFigure,plotTessellation=plotTessellation,Ei_err=Ei_err,temperature_err=temperature_err,\
        magneticField_err=magneticField_err,electricField_err=electricField_err)
#
#    def plotQPatches(self,dataFiles=None,ax=None,planes=[],binningDecimals=3,log=False,returnPatches=False,A4Extend=0.2,A3Extend=0.5,singleFigure=False,plotTessellation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
#        """Plot data files together with pixels created around each point in Q space. 
#
#        .. warning::
#           This method plots all measurement points unless they are literaly on top of each other and is thus really slow! Binning 8 planes for two files takes approximately
#           3.5 minutes. Alternatively use binning, i.e. plotQPlane.
#
#        Kwargs:
#
#            - dataFiles (DataFiles): single file or list of files to be binned together (Default self.convertedFiles)
#
#            - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)
#
#            - planes (list (of lists)): Planes to be plotted and binned (default [])
#
#            - binningDecimals (int): Number of decimal places Q positions are rounded before binning (default 3)
#            
#            - log (bool): Whether or not to plot intensities as logarithm (default False)
#
#            - returnPatches (bool): If true the method returns the patches otherwise plotted in the given axis (default False)
#
#            - A4Extend (float): Angle value with which the boundary is extended away from points in A4 direction (default 0.2)
#            
#            - A3Extend (float): Angle value with which the boundary is extended away from points in A3 direction (default 0.5)
#
#            - singleFigure (bool): If true, all planes are plotted in same figure (default False)
#
#            - plotTessellation (bool): Plot Tessellation of points (default False)
#
#            - Ei_err (float): Tolerance of E_i for which the values are equal (default = 0.05)
#
#            - temperature_err (float): Tolerance of temperature for which the values are equal (default = 0.2)
#            
#            - magneticField_err (float): Tolerance of magnetic field for which the values are equal (default = 0.2)
#            
#            - electricField_err (float): Tolerance of electric field for which the values are equal (default = 0.2)
#
#        Returns:
#            
#            - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.
#
#        Raises:
#
#            - AttributeError
#
#        The following example will combine the two files and plot all of the available planes in different figures.
#
#        >>> DS = DataSet.DataSet(convertedFiles=[--.nxs,---.nxs])
#        >>> plt.figure()
#        >>> ax = plt.gca()
#        >>>
#        >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax)
#
#        If only a subset of planes or different planes are to be combined the following will achieve this:
#
#        >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax,planes=[0,1,2,3,[4,5,6],[8,9]])
#
#        Here planes 0 through 3 are plotted separately while 4,5, and 6 as well as 8 and 9 are binned.
#
#        .. note::
#            Binning planes from different analysers might result in nonsensible binnings.
#
#        """
#        if dataFiles is None:
#            dataFiles = self.convertedFiles
#        
#        return plotQPatches(dataFiles,ax=ax,planes=planes,binningDecimals=binningDecimals,log=log,returnPatches=returnPatches,A4Extend=A4Extend,A3Extend=A3Extend,singleFigure=singleFigure,\
#        plotTessellation=plotTessellation,Ei_err=Ei_err,temperature_err=temperature_err,\
#        magneticField_err=magneticField_err,electricField_err=electricField_err)

    
    @_tools.KwargChecker()
    def cutQELine(self,QPoints,EnergyBins,width=0.1,minPixel=0.01,rlu=True,dataFiles=None,constantBins=False):
        """
        Method to perform Q-energy cuts from a variable number of points. The function takes both qx/qy or hkl positions. In the case of using only two Q points,
        the method is equivalent to cutQE.
        
        Args:
            
            - QPoints (list of points): Q positions between which cuts are performed. Can be specified with both qx, qy or hkl positions dependent on the choice of format.
            
            - EnergyBins (list of floats): Energy bins for which the cuts are performed
            
        Kwargs:
        
            - width (float): Width of the cut in 1/AA (default 0.1).
            
            - minPixel (float): Minimal size of binning along the cutting directions. Points will be binned if they arecloser than minPixel (default=0.01)
        
            - rlu (bool): If True, provided QPoints are interpreted as (h,k,l) otherwise as (qx,qy), (default True).
        
            - dataFiles (list): List of dataFiles to cut. If none, the ones in the object will be used (default None).

            - constantBins (bool): If True only bins of size minPixel is used (default False)
        
        .. warning::
            The way the binning works is by extending the end points with 0.5*minPixel, but the method sorts away points not between the two Q points given and thus the start and end
            bins are only half filled. This might result in discrepancies between a single cut and the same cut split into different steps. Further, splitting lines into sub-cuts 
            forces a new binning to be done and the bin positions can then differ from the case where only one cut is performed.

        
        Returns: m = Q points, n = energy bins
                
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for all 2D cuts.
            
            - Bin list (m * n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (m * n * 3D arrays): n instances of center positions for the bins.

            - binDistance (m * n arrays): n instances of arrays holding the distance in q to q1.

        .. note::
            If an HKL point outside of the scattering plane is given, the program will just take the projection onto the scattering plane.
            
        """
        if not isinstance(QPoints,np.ndarray):
            QPoints = np.array(QPoints)

        if(len(QPoints)<2):
            raise AttributeError('Number of Q points given is less than 2.')
        if rlu==True: # Recalculate q points into qx and qy points
        #    sample =self.sample[0]
        #    positions = self.convertToQxQy(QPoints)
            pass
            
        elif rlu==False: # RLU is false
        #    positions = QPoints
            if QPoints.shape[1]!=2:
                raise AttributeError('Provide Q list is not 2 dimensional, should have shape (n,2) in QxQy mode but got shape {}.'.format(QPoints.shape))
        else:
            raise AttributeError('Given Q mode not understood. Got {} but must be either "RLU", "HKL" or "QxQy"')

        if EnergyBins.shape == ():
            EnergyBins = np.array([EnergyBins])

        if len(EnergyBins.shape)==1 and not isinstance(EnergyBins[0],(list,np.ndarray)):
            EnergyBins = np.array([EnergyBins for _ in range(len(QPoints)-1)]).reshape(len(QPoints)-1,-1)

        if not isinstance(width,(list,np.ndarray)):
            width = np.array([width for _ in range(len(QPoints)-1)]).reshape(len(QPoints)-1)

        if not isinstance(minPixel,(list,np.ndarray)):
            minPixel = np.array([minPixel for _ in range(len(QPoints)-1)]).reshape(len(QPoints)-1)

        DataList = []
        BinList = []
        centerPosition = []
        binDistance = []

        for cutIndex,[pStart,pStop,w,mP,EB] in enumerate(zip(QPoints,QPoints[1:],width,minPixel,EnergyBins)):
            _DataList,_BinList,_centerPosition,_binDistance = self.cutQE(q1=pStart,q2=pStop,width=w,minPixel=mP,EnergyBins=EB,rlu=rlu,
                                                                         dataFiles=dataFiles,extend=False,constantBins=constantBins)
            _DataList['qCut']=cutIndex
            DataList.append(_DataList)
            if rlu:
                UB2D = self.sample[0].convertHKLINV # Matrix to calculate HKL from Qx,Qy 
                _BinListUpdated = []
                _centerPositionUpdated = []
                for i,[Position,ortho,E] in enumerate(_BinList):
                    pos = np.array([np.concatenate([np.dot(UB2D,x[:2]),[x[2]]],axis=0) for x in Position])
                    orthogonal = [np.dot(UB2D,x) for x in ortho]
                    _BinListUpdated.append([pos,orthogonal,E])

                    cPos = np.array([np.concatenate([np.dot(UB2D,x[:2]),[x[2]]],axis=0) for x in _centerPosition[i]])
                    _centerPositionUpdated.append(cPos)
                _BinList = _BinListUpdated
                _centerPosition = _centerPositionUpdated
            BinList.append(_BinList)
            centerPosition.append(_centerPosition)
            binDistance.append(_binDistance)
            
        DataList = pd.concat(DataList)
        return DataList,np.array(BinList,dtype=object),np.array(centerPosition,dtype=object),np.array(binDistance,dtype=object)

    
    @_tools.KwargChecker(include=np.concatenate([_tools.MPLKwargs,['vmin','vmax','log','ticks','seperatorWidth','tickRound','plotSeperator','cmap','colorbar','edgecolors']]))
    def plotCutQELine(self,QPoints,EnergyBins,width=0.1,minPixel=0.01,rlu=True,ax=None,dataFiles=None,constantBins=False,outputFunction=print,**kwargs):
        """Plotting wrapper for the cutQELine method. Plots the scattering intensity as a function of Q and E for cuts between specified Q-points.
        
        Args:
            
            - QPoints (list): List of Q points in either RLU (3D) or QxQy (2D).
            
            - EnergyBins (list): List of bin edges in the energy direction.
        
        Kwargs:
            
            - width (float): Width perpendicular to Q-direction for cuts (default 0.1)

            - minPixel (float): Minimum size of pixel for cut (default 0.01)
            
            - rlu (bool): If True, provided points are interpreted as (h,k,l) otherwise (qx,qy), (default RLU)
            
            - ax (matplotlib axis): Axis into whiht the data is plotted. If None a new will be created (default None).
            
            - dataFiles (DataFile(s)): DataFile or list of, from which data is to be taken. If None all datafiles in self is taken (default None).
            
            - vmin (float): Lower limit for colorbar (default min(Intensity)).
            
            - vmax (float): Upper limit for colorbar (default max(Intensity)).
            
            - tickRound (int): Number of decimals ticks are rounded to (default 3).
            
            - ticks (int): Number of ticks in plot, minimum equal to number of Q points (default 8).
            
            - plotSeperator (bool): If true, vertical lines are plotted at Q points (default True).
            
            - seperatorWidth (float): Width of seperator line (default 2).
            
            - log (bool): If true the plotted intensity is the logarithm of the intensity (default False)

            - constantBins (bool): If True only bins of size minPixel is used (default False)

            - outputFunction (function): Function called on output string (default print)

        Return:  m = Q points, n = energy bins
            
            - ax: matplotlib axis in which the data is plotted
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for all 2D cuts.
                
            - Bin list (m * n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (m * n * 3D arrays): n instances of center positions for the bins.

            - binDistance (m * n arrays): n instances of arrays holding the distance in q to q1.

        .. note::
            
            The ax.set_clim function is created to change the colour scale. It takes inputs vmin,vmax. This function does however not work in 3D....

        """
        if not isinstance(EnergyBins,np.ndarray):
            EnergyBins = np.array(EnergyBins)

        DataList,BinListTotal,centerPositionTotal,binDistanceTotal = self.cutQELine(QPoints=QPoints,EnergyBins=EnergyBins,width=width,minPixel=minPixel,rlu=rlu,dataFiles=dataFiles,constantBins=constantBins)
        if rlu==True: # Recalculate q points into qx and qy points
            positions = self.convertToQxQy(QPoints)
        else: # Do nothing
            positions = QPoints
        
        if ax is None:
            plt.figure()
            ax = plt.gca()
            _3D = False
        else:
            ax = ax
            if not hasattr(ax,'name'):
                _3D = False
            else:
                if ax.name == '3d':
                    _3D = True
                else:
                    _3D = False

        if not 'log' in kwargs:
            log = False
        else:
            log = kwargs['log']
            del kwargs['log']

        if 'vmin' in kwargs:
            vmin = kwargs['vmin']
        elif hasattr(ax,'vmin'):
            vmin = ax.vmin
        else:
            if log==True:
                vmin = np.log10(DataList['Int'].min()+1e-20)
            else:
                vmin = DataList['Int'].min()
            
        
        if 'vmax' in kwargs:
            vmax = kwargs['vmax']
        elif hasattr(ax,'vmax'):
            vmax = ax.vmax
        else:
            if log==True:
                vmax = np.log10(DataList['Int'].max()+1e-20)
            else:
                vmax = DataList['Int'].max()


        if _3D == False:
            if len(EnergyBins.shape)==1 and not isinstance(EnergyBins[0],(list,np.ndarray)): # If a common energy binning is requested
                if not len(QPoints)==2:
                    EnergyBins = np.array([EnergyBins for x in range(len(QPoints)-1)])
                else:
                    EnergyBins = np.array(EnergyBins).reshape(1,-1)

            emptyCuts = [len(cut)==0 for cut in binDistanceTotal]
            emptyIndex = np.arange(len(binDistanceTotal),dtype=int)[emptyCuts]
            if len(emptyIndex)!=0:
                string = ['No data points found between {} and {} with energies between {} and {}.'.format(QPoints[idx],QPoints[idx+1],EnergyBins[0],EnergyBins[-1]) for idx in emptyIndex]
                raise AttributeError('\n'.join([x for x in string]))
            
            BinNums = [np.max([np.max(binDistanceTotal[i][j]) for j in range(len(binDistanceTotal[i]))]) for i in range(len(binDistanceTotal))]
            
            if not 'ticks' in kwargs:
                ticks = 8
            else:
                ticks = kwargs['ticks']
                kwargs = _tools.without_keys(dictionary=kwargs,keys='ticks')
                
            
            NumQPointTicks = len(QPoints)
            if NumQPointTicks > ticks: # If there are more QPoints than ticks, make them match
                ticks=NumQPointTicks
            
            freeTicks = ticks-NumQPointTicks
            
            # Figure out how many ticks are to be placed in each segment for quasi-equal spacing    
            ticksInSegments = np.ones((len(BinNums)),dtype=int)
            for i in range(freeTicks):
                TickDistance = BinNums/ticksInSegments
                arg = np.argmax(TickDistance)
                ticksInSegments[arg]+=1

            if not 'tickRound' in kwargs:
                tickRound = 3
            else:
                tickRound = kwargs['tickRound']
                del kwargs['tickRound']
                
            if not 'plotSeperator' in kwargs:
                plotSeperator = True
            else:
                plotSeperator = kwargs['plotSeperator']
                del kwargs['plotSeperator']
                
            if not 'seperatorWidth' in kwargs:
                seperatorWidth = 2
            else:
                seperatorWidth = kwargs['seperatorWidth']
                del kwargs['seperatorWidth']
                
            if not 'colorbar' in kwargs:
                colorbar = True
            else:
                colorbar = kwargs['colorbar']
                del kwargs['colorbar']
                
            

            pmeshs = []
            xticks = []
            xticklabels = []
            IntTotal = []
            
            offset = [0.0]
            idmax = len(BinListTotal) # Number of segments, when assuming no reordering of points in each segment
            #DataList['Int'] = DataList['Intensity']*DataList['BinCount']/(DataList['Monitor']*DataList['Normalization'])
            edgeQDistance = []
            actualEnergy = []
            qstart = []
            qstop = []
            direction = []
            distanceChange = []

            for BLT in BinListTotal:
                localE = []
                for BLTSub in BLT:
                    localE.append(BLTSub[2][0])
                
                localE.append(BLT[-1][2][1])
                
                actualEnergy.append(localE)
                
            for segID in range(idmax): # extract relevant data for current segment
                #[intensityArray,monitorArray,normalizationArray,normcountArray] = DataList[segID]#[i] for i in range(4)]
                BinList = BinListTotal[segID]
                centerPosition = centerPositionTotal[segID]
                binDistance = binDistanceTotal[segID]
                q1 = positions[segID]
                q2 = positions[segID+1]
                if rlu:
                    q1 = np.dot(self.sample[0].convertHKLINV,q1) 
                    q2 = np.dot(self.sample[0].convertHKLINV,q2)
                dirvec = np.array(q2)-np.array(q1)
                
                leftEdgeBins = np.array([np.dot(BL[0][0,:len(q1)],dirvec) for BL in BinList])
            
                leftEdgeIndex = np.argmin(leftEdgeBins)

                binedges = [np.linalg.norm(BL[0][:,:len(q1)]-BinList[leftEdgeIndex][0][0,:len(q1)],axis=1) for BL in BinList]
                
                binenergies = [BL[2] for BL in BinList]

                edgeQDistanceLocal = []
                for BL in BinList:
                    p = BL[0][:,:len(q1)] - q1
                    #print(p)
                    #print(q1)
                    q = np.dot(p, dirvec)
                    #print(q)
                    if not (np.sort(q) == q).all():
                        raise RuntimeError("edgeQDistance[{}] is not sorted".format(segID))
                    edgeQDistanceLocal.append(q)
                edgeQDistance.append(edgeQDistanceLocal)

                localDataList = DataList[DataList['qCut']==segID]
                # plot intensity with provied kwargs
                for [_ ,intensity],binEdge,binEnergy in zip(localDataList[['Int','energyCut']].groupby('energyCut'),binedges,binenergies):
                    if log==True:
                        pmeshs.append(ax.pcolormesh(binEdge+offset[-1],binEnergy,np.log10(intensity.T+1e-20),**kwargs))
                    else:
                        pmeshs.append(ax.pcolormesh(binEdge+offset[-1],binEnergy,intensity.T,**kwargs))   
                if plotSeperator == True:
                    plt.plot([offset[-1],offset[-1]],[np.min(EnergyBins[segID]),np.max(EnergyBins[segID])],'k',linewidth=seperatorWidth)

                # Find bounding Q-points
                minimalDistanceIDEnergy = np.argmin([np.min(x) for x in binDistance])
                minimalDistanceID = np.argmin(binDistance[minimalDistanceIDEnergy])
                maximalDistanceIDEnergy = np.argmax([np.max(x) for x in binDistance])
                maximalDistanceID = np.argmax(binDistance[maximalDistanceIDEnergy])
                
                
                qstart.append(centerPosition[minimalDistanceIDEnergy][minimalDistanceID][:len(q1)])
                qstop.append(centerPosition[maximalDistanceIDEnergy][maximalDistanceID][:len(q1)])
                # Prepare the calculation of tick markers
                direction.append(qstop[segID]-qstart[segID])
                distanceChange.append(np.max(binDistance[maximalDistanceIDEnergy])-np.min(binDistance[minimalDistanceIDEnergy]))
                if rlu:
                    qstartQ,qstopQ = [np.dot(self.sample[0].convertHKL,x) for x in [qstart[segID],qstop[segID]]]
                    dirLen = np.linalg.norm(direction[segID])
                    dirLenQ = np.linalg.norm(qstopQ-qstartQ)
                    distanceChange[-1]*=dirLen/dirLenQ
                
                binminmaxList = np.linspace(0,1,100)
                ticksCurrentSeg = ticksInSegments[segID]
                

                num=len(binminmaxList)
                
                if segID==idmax-1: ## Find best bins for tick marking
                    xvalues = np.round(np.linspace(0,num-1,ticksCurrentSeg+1)).astype(int) 
                else:
                    xvalues = np.round(np.linspace(0,num-1-int(num/ticksCurrentSeg),ticksCurrentSeg)).astype(int)

                my_xticks=[]
                for i in xvalues:
                    positionValues = binminmaxList[i]*direction[segID]+qstart[segID]
                    my_xticks.append('\n'.join([('{:.'+str(tickRound)+'f}').format(x+0.0) for x in positionValues]))

                
                
                xticks.append(binminmaxList[xvalues]*distanceChange[segID] +offset[segID]) # binDistanceAll[xvalues]
                xticklabels.append(my_xticks)
                offset.append(offset[segID]+np.max([np.max(binedge) for binedge in binedges]))
            def calculateXPosition(x,ID):
                return np.dot(x-qstart[ID],direction[ID])/(np.linalg.norm(direction[ID])**2)*distanceChange[ID] +offset[ID]
            ax.converterFunction = calculateXPosition
            if plotSeperator == True: # plot last seperator
                plt.plot([offset,offset],[np.min(EnergyBins[-1]),np.max(EnergyBins[-1])],'k',linewidth=seperatorWidth)
            
            ax.set_xticks(np.concatenate(xticks))
            ax.set_xticklabels(np.concatenate(xticklabels), multialignment="center",ha="center")
            ax.EnergyBins = EnergyBins
            if rlu==True:
                ax.set_xlabel('$Q_h$ [RLU]\n$Q_k$ [RLU]\n$Q_l$ [RLU]')
            else:
                ax.set_xlabel(r'$Q_x [\AA^{-1}]$'+'\n'+r'$Q_y [\AA^{-1}]$')
            
            ax.xaxis.set_label_coords(1.15, -0.025) 
            ax.set_ylabel('E [meV]')
            if colorbar:
                ax.colorbar = ax.get_figure().colorbar(pmeshs[0],pad=0.1,format='%.2E')
                ax.colorbar.set_label('$I$ [arb.u.]', rotation=270)
            #plt.tight_layout()
            if 'pmeshs' in ax.__dict__:
                ax.pmeshs = np.concatenate([ax.pmeshs,pmeshs],axis=0)
            else:
                ax.pmeshs = pmeshs
        
            # Create mouse over function
            def format_coord(x,y,edgeQDistance,centerPos,EnergyBins,DataList,rlu,offset,self):# pragma: no cover
                val = calculateIndex(x,y,offset,EnergyBins,edgeQDistance,DataList,textReturn=True)
                if type(val)==str:
                    return val
                segID,Eindex,index = val
                dfOffset = DataList['energyCut'][DataList['qCut']==segID].min()

                Intensity = DataList['Int'][np.logical_and(DataList['qCut']==segID,DataList['energyCut']==Eindex+dfOffset)][index]  #Int[segID][Eindex][index][0]
                if rlu==False: 
                    qx, qy, E = centerPos[segID][Eindex][index]
                    return "qx = {0:.3f}, qy = {1:.3f}, E = {2:.3f}, I = {3:.3e}".format(qx+0.0,qy+0.0,E,Intensity)
                else:
                    H,K,L, E = centerPos[segID][Eindex][index]
                    return "h = {0:.3f}, h = {1:.3f}, l = {2:.3f}, E = {3:.3f}, I = {4:.3e}".format(H+0.0,K+0.0,L+0.0,E,Intensity)

            def calculateIndex(x,y,offset,EnergyBins,edgeQDistance,DataList,textReturn): # pragma: no cover
                if x<offset[0] or x>=offset[-1]:
                    if textReturn == True:
                        return "x out of range: {:.3}".format(x)
                    else:
                        return -1,-1,-1
                try:
                    segID = (np.arange(len(offset)-1)[[x>offStart and x<offStop for offStart,offStop in zip(offset,offset[1:])]])[0]
                except IndexError: # index error returned if no segment is found
                    if textReturn == True:
                        return "x out of range: {:.3}".format(x)
                    else:
                        return -1,-1,-1
                Eindex = np.array(EnergyBins[segID]).searchsorted(y) - 1
                minspan = np.min(np.concatenate(edgeQDistance[segID]))
                maxspan = np.max(np.concatenate(edgeQDistance[segID]))
                xInSegment = (x-offset[segID])/(offset[segID+1]-offset[segID])*(maxspan-minspan)+minspan
                x = xInSegment
                if len(EnergyBins[segID]) < 2:
                    if textReturn == True:
                        return "len(EnergyBins[{}]) < 2".format(segID)
                    else:
                        return -1,-1,-1
                    
                if y < EnergyBins[segID][0] or y >= EnergyBins[segID][-1]:
                    if textReturn == True:
                        return "E out of range {:.3}  {}..{}".format(y, EnergyBins[0], EnergyBins[-1])
                    else:
                        return -1,-1,-1
                if len(edgeQDistance[segID][Eindex]) < 2:
                    raise RuntimeError("len(edgeQDistance[{}][{}]) < 2".format(segID,Eindex))
                if x < edgeQDistance[segID][Eindex][0] or x >= edgeQDistance[segID][Eindex][-1]:
                    if textReturn == True:
                        return "x out of range: {:.3}".format(x)
                    else:
                        return -1,-1,-1

    
                index = edgeQDistance[segID][Eindex].searchsorted(xInSegment) - 1
                
                return segID,Eindex,index
            ax.calculateIndex = lambda x,y: calculateIndex(x,y,offset,actualEnergy,edgeQDistance,DataList,textReturn=False)#EnergyBins
            def onclick(event,ax,DataList,outputFunction): # pragma: no cover
                if ax.in_axes(event):
                    try:
                        C = ax.get_figure().canvas.cursor().shape() # Only works for pyQt5 backend
                    except:
                        pass
                    else:
                        if C != 0:
                            return
                    x = event.xdata
                    y = event.ydata
                    printString = ax.format_coord(x,y)
                    segID,Eindex,index = ax.calculateIndex(x,y)
                    dfOffset = DataList['energyCut'][DataList['qCut']==segID].min()
                    Eindex+=dfOffset
                    if index < len(DataList[np.logical_and(DataList['qCut']==segID,DataList['energyCut']==Eindex)]) and index>=0:
                        dataPoint = DataList[np.logical_and(DataList['qCut']==segID,DataList['energyCut']==Eindex)][index]
                        if not np.any([x==-1 for x in [segID,Eindex,index]]):
                            cts = int(dataPoint['Intensity'])
                            Mon = int(dataPoint['Monitor'])
                            Norm = float(dataPoint['Normalization'])
                            NC = int(dataPoint['BinCount'])
                            printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
                    
                    outputFunction(printString)
            ax.format_coord = lambda x,y: format_coord(x,y,edgeQDistance,centerPositionTotal,actualEnergy,DataList,rlu,offset,self)#EnergyBins
            ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,DataList,outputFunction=outputFunction))
            ax.edgeQDistance = edgeQDistance   
            ax.offset = offset  
            
        else: 
            def set_clim_local(self,vmin,vmax):
                color=list(ax.cmap(ax.norm(self.value,vmin,vmax)))
                self.set_facecolor(color)
                self.set_edgecolor(color)
            
            def norm(x,vmin,vmax):
                return np.divide(x-vmin,vmax)
            
            ax.norm = norm
            if not 'cmap' in kwargs:
                #from matplotlib.colors import ListedColormap
                ax.cmap = plt.cm.coolwarm
            else:
                ax.cmap = kwargs['cmap']
                kwargs = _tools.without_keys(dictionary=kwargs,keys='cmap')
            sfp = []
            for bins,[_,datlist] in zip(BinListTotal,DataList.groupby('qCut')):
                energies = len(bins)
                
                energyEdges = np.array([bins[idx][2] for idx in range(energies)])
                ELength = np.array([len(x[0][:,0]) for x in bins])
                ELengthCummu = np.concatenate([[0],np.cumsum(ELength)],axis=0)
                H = np.concatenate([bins[idx][0][:,0] for idx in range(energies)],axis=0)
                K = np.concatenate([bins[idx][0][:,1] for idx in range(energies)],axis=0)
                L = np.concatenate([bins[idx][0][:,2] for idx in range(energies)],axis=0)
                P0,P1 = self.sample[0].calculateHKLToQxQy(H,K,L)  
                P0,P1 = np.einsum('mj,j...->m...',self.sample[0].RotMat,[P0,P1])
            
                IntCommu = np.concatenate([[0],np.cumsum(ELength-1)],axis=0)
                
                E = np.concatenate([bins[idx][0][:,3] for idx in range(energies)],axis=0)
                EBins = np.array([bins[idx][2] for idx in range(energies)])
                
                for E in range(energies):
                    X = P0[ELengthCummu[E]:ELengthCummu[E+1]].reshape(-1,1).repeat(2,axis=1)
                    Y = P1[ELengthCummu[E]:ELengthCummu[E+1]].reshape(-1,1).repeat(2,axis=1)
                    Z = np.ones_like(X)*EBins[E].reshape(1,2)
                    
                    if IntCommu[E+1]-IntCommu[E]==0: # If segment is empty
                        continue
                    
                    #normColor = np.divide(INT[IntCommu[i]:IntCommu[i+1]]-vmin,vmax-vmin).reshape(-1,1).repeat(2,axis=1)
                    normColor = datlist['Int'][IntCommu[E]:IntCommu[E+1]].values#.reshape(-1,1).repeat(2,axis=1)
                    
                    #if len(normColor)!=X.shape[0]:
                    #    continue
                    color=list(ax.cmap(ax.norm(normColor,vmin,vmax)))
                    
                    
                    sf = ax.plot_surface(X,Y,Z,rstride=1,cstride=1,shade=False)
                    sf.value = normColor
                    sf.set_facecolor(color)
                    sf.set_edgecolor(color)
                    sf.ax = ax
                    # TODO: make colorlimits changeable after creation of axis
                    sf.name = str(E)
                    
                    #sf.set_clim = lambda vmin,vmax: set_clim_local(self=sf,vmin=vmin,vmax=vmax,ax=ax)
                    sfp.append(sf)
            
            #if rlu:
            #    ax.set_xlabel(_tools.generateLabel(self.sample[0].projectionVector1)+ ' [RLU]')
            #    ax.set_ylabel(_tools.generateLabel(self.sample[0].projectionVector2)+ ' [RLU]')
            #else:
            ax.set_xlabel(r'Qx [$\AA^{-1}$]')
            ax.set_ylabel(r'Qy [$\AA^{-1}$]')
            ax.set_zlabel('E [meV]')
    
            setattr(sf.__class__,'set_clim',set_clim_local)
    
            pmeshs = np.array(sfp).flatten()
            if 'pmeshs' in ax.__dict__:
                ax.pmeshs = np.concatenate([ax.pmeshs,pmeshs],axis=0)
            else:
                ax.pmeshs = pmeshs

        ax.vmin = vmin
        ax.vmax = vmax
        def set_clim(vmin,vmax,ax):
            ax.vmin = vmin
            ax.vmax = vmax
            for pm in ax.pmeshs:
                pm.set_clim(vmin,vmax)

        ax.get_clim = lambda: (ax.vmin,ax.vmax)

        ax.set_clim = lambda vmin,vmax:set_clim(vmin,vmax,ax)
        ax.set_clim(*ax.get_clim())

        
        
        return ax,DataList,BinListTotal,centerPositionTotal,binDistanceTotal

    @_tools.KwargChecker()
    def extractData(self, A4 = None, A4Id = None, Ef = None, EfId = None, raw = False, A4Tolerance = 0.1, EfTolerance = 0.1):
        """Extract data given A4 value and Ef (or the corresponding indices).
        
        Kwargs:
            
            - A4 (float): Wanted A4 value in degrees (default None)
            
            - A4Id (int): Id of wedge which is a number between 0 and 103 (default None)
            
            - Ef (float): Wanted Ef value in meV (default None)
            
            - EfId (int): Wanted Id of analyser energy, number between 0-7 (default None)
            
            - raw (bool): If true method returns Intensity,Normalization,Monitor, else returns Intensity/(Norm*Monitor) (default False)
            
            - A4Tolerance (float): Tolerance between found and wanted A4 value in degrees (default 0.1)

            - EfTolerance (float): Tolerance between found and wanted Ef value in meV (default 0.1)
            
        .. note::
            If A4 or Ef is provided, then these will be used instead of A4Id or EfId.
            
        """
        if raw: # Shape is (1 or 3, no Files, steps, 104, binning)
            Data = np.array([self.I,self.Norm,self.Monitor])
        else:
            Data = np.divide(self.I,self.Norm*self.Monitor)
        
        if A4 is None and A4Id is None and Ef is None and EfId is None:
            return Data
        returnData = []
        for file in self.convertedFiles:
            
            binning = file.binning
            if not binning == 1:
                raise AttributeError('Provided file has a binning different from 1. This is not supported for data extraction as one is not allowed to use the prismatic concept for alignment...')
                
            # Find A4 id if not given
            if not A4 is None:
                NominalA4 = file.instrumentCalibrationA4.reshape(104,8*binning)
                A4Id,A4Analyser = np.unravel_index(np.argmin(np.abs(NominalA4-file.A4Off-A4)),NominalA4.shape) # Extract only wedge number
                A4Found = (NominalA4[A4Id,A4Analyser]-file.A4Off)[0]
                if np.abs(A4Found-A4)>A4Tolerance:
                    raise AttributeError('Difference between wanted ({}) and found A4 value ({}) is {} for file {}.. Maybe the sign of A4 should be changed.'.format(A4,A4Found,A4-A4Found,file.name))
            

            # Find A4 id if not given
            if not Ef is None:
                if not A4Id is None:
                    NominalEf = file.instrumentCalibrationEf[:,1].reshape(104,8*binning)[A4Id].reshape(1,-1)
                else:
                    NominalEf = file.instrumentCalibrationEf[:,1].reshape(104,8*binning)
                EfDetector,EfId = np.unravel_index(np.argmin(np.abs(NominalEf-Ef)),NominalEf.shape) # Extract only wedge number
                EfFound = NominalEf[EfDetector,EfId]
                if np.abs(EfFound-Ef)>EfTolerance:
                    raise AttributeError('Difference between wanted ({}) and found Ef value ({}) is {} for file {}.. Maybe the sign of A4 should be changed.'.format(Ef,EfFound,Ef-EfFound,file.name))

            
            if raw: # Shape is (1 or 3, no Files, steps, 104, binning)
                Data = np.array([file.I,file.Norm,file.Monitor])
            else:
                Data = np.divide(file.I,file.Norm*file.Monitor)
            
            if not A4Id is None:
                rData = Data[:,A4Id].reshape(-1,1,8*binning)
            else:
                rData = Data
                
            if not EfId is None:
                rData = rData[:,:,EfId]
                
            returnData.append(np.squeeze(rData))

        return returnData

    
    @_tools.KwargChecker()
    def cut1DE(self,E1,E2,q,rlu=True,width=0.02, minPixel = 0.1, dataFiles = None,constantBins=False,ufit=False):
        """Perform 1D cut through constant Q point returning binned intensity, monitor, normalization and normcount. The width of the cut is given by 
        the width attribute.
        
        
        Args:
            
            - E1 (float): Start energy.
            
            - E2 (float): End energy.

            - q (3D or 2D vector): Q point 
        
        Kwargs:
            
            - rlu (bool): If True, provided Q point is interpreted as (h,k,l) otherwise as (qx,qy), (Default true)

            - width (float): Full width of cut in q-plane (default 0.02).
            
            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel (default 0.1).
            
            - dataFiles (list): Data files to be used. If none provided use the ones in self (default None)

            - constantBins (bool): If True only bins of size minPixel is used (default False)

            - ufit (bool): If True a uFit Dataset object is returned in stead of pandas data frame
            
        Returns:
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for 1D cut.
            
            - Bin list (1 array): Bin edge positions in energy

        """
        if dataFiles is None:
                if len(self.convertedFiles)==0:
                    raise AttributeError('No data file to be binned provided in either input or DataSet object.')
                else:
                    I = self.I
                    qx = self.qx
                    qy = self.qy
                    energy = self.energy
                    Norm = self.Norm
                    Monitor = self.Monitor
                    sample = self.convertedFiles[0].sample

        else: 
            #dataFiles = isListOfDataFiles(dataFiles)
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor = DS.I,DS.qx,DS.qy,DS.energy,DS.Norm,DS.Monitor
            sample = DS.convertedFiles[0].sample
            
        I = np.concatenate([x.flatten() for x in I])
        qx = np.concatenate([x.flatten() for x in qx])#self.qx
        qy = np.concatenate([x.flatten() for x in qy])#self.qy
        energy = np.concatenate([x.flatten() for x in energy])#self.energy
        Norm = np.concatenate([x.flatten() for x in Norm])#self.Norm
        Monitor = np.concatenate([x.flatten() for x in Monitor])#self.Monitor
        positions = [qx,qy,energy]

        if rlu==True: # Recalculate q points into qx and qy points
            Q = self.convertToQxQy(q).flatten()
            variables = ['H','K','L']
        else: # Do nothing
            Q = np.array(q).flatten()
            variables = ['Qx','Qy']
        variables.append('Energy')
        [intensity,MonitorCount,Normalization,normcounts],bins  = cut1DE(positions = positions, I=I, Norm=Norm,Monitor=Monitor,E1=E1,E2=E2,q=Q,width=width,minPixel=minPixel,constantBins=constantBins)
        data = pd.DataFrame()
        

        HKL = self.convertToHKL(Q.flatten())
        data['Qx'] = Q[0]*np.ones_like(intensity)
        data['Qy'] = Q[1]*np.ones_like(intensity)
        data['H'] = HKL[0]*np.ones_like(intensity)
        data['K'] = HKL[1]*np.ones_like(intensity)
        data['L'] = HKL[2]*np.ones_like(intensity)
        data['Energy'] = 0.5*(bins[0][1:]+bins[0][:-1])
        data['Intensity'] = intensity.astype(int)
        data['Monitor'] = MonitorCount.astype(int)
        data['Normalization'] = Normalization.astype(float)
        data['BinCount'] = normcounts.astype(int)
        data['binDistance'] = np.linalg.norm(data[variables]-np.array(data[variables].iloc[1]),axis=1)
        
        data['Int'] = data['Intensity']*data['BinCount']/(data['Normalization']*data['Monitor'])
        if not ufit:
            return data,bins
        
        ufitData = self.generateUFitDataset(data,q1=q,q2=None,rlu=rlu,width=width,minPixel=minPixel,Emin=E1,Emax=E2,QDirection=False)
        
        return ufitData

    @_tools.KwargChecker(function=plt.errorbar,include=np.concatenate([_tools.MPLKwargs,['ticks','tickRound','mfc','markeredgewidth','markersize']])) #Advanced KWargs checker for figures
    def plotCut1DE(self,E1,E2,q,rlu=True,width=0.02, minPixel = 0.1, dataFiles = None,constantBins=False,ax=None,ufit=False,**kwargs):
        """Perform 1D cut through constant Q point returning binned intensity, monitor, normalization and normcount. The width of the cut is given by 
        the width attribute.
        
        
        Args:
            
            - E1 (float): Start energy.
            
            - E2 (float): End energy.

            - q (3D or 2D vector): Q point 
        
        Kwargs:
            
            - rlu (bool): If True, provided Q point is interpreted as (h,k,l) otherwise as (qx,qy), (Default true)

            - width (float): Full width of cut in q-plane (default 0.02).
            
            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel (default 0.1).
            
            - dataFiles (list): Data files to be used. If none provided use the ones in self (default None)

            - constantBins (bool): If True only bins of size minPixel is used (default False)

            - ufit (bool): If True a uFit Dataset object is returned in stead of pandas data frame

            - ax (matplotlib.axes): If None, new axis is created in which to plot (default None)

            - **kwargs: Pass on to the ax.errorbar method used to plot
            
        Returns:
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for 1D cut.
            
            - Bin list (1 array): Bin edge positions in energy

        """
        Data,bins = self.cut1DE(E1=E1,E2=E2,q=q,rlu=rlu,width=width, minPixel =  minPixel,dataFiles = dataFiles,constantBins=constantBins)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            INT = np.divide(Data['Intensity']*Data['BinCount'],Data['Monitor']*Data['Normalization'])
            INT_err = np.divide(np.sqrt(Data['Intensity'])*Data['BinCount'],Data['Monitor']*Data['Normalization'])
        


        if rlu:
            variables = ['H','K','L']
        else:
            variables = ['Qx','Qy']
        
        variables = variables+['Energy']

        num = len(Data)
        if not 'ticks' in kwargs:
            ticks = 5
        else:
            ticks = kwargs['ticks']
            kwargs = _tools.without_keys(dictionary=kwargs, keys='ticks')

        if not 'tickRound' in kwargs:
            tickRound = 3
        else:
            tickRound = kwargs['tickRound']
            kwargs = _tools.without_keys(dictionary=kwargs,keys='tickRound')

        if not 'fmt' in kwargs:
            kwargs['fmt'] = '.'

        xvalues = np.round(np.linspace(0,num-1,ticks)).astype(int)
        my_xticks=[]
        for i in xvalues:
            my_xticks.append('\n'.join(map(lambda x:('{:.'+str(tickRound)+'f}').format(x),[np.round(Data[var][i],tickRound) for var in variables])))
        

        Data['binDistance'] = np.linalg.norm(Data[variables]-np.array(Data[variables].iloc[1]),axis=1)

        if ax is None:
            plt.figure()
            ax = plt.gca()

        if not 'label' in kwargs:
            kwargs['label'] = 'Data'
        ax.errorbar(Data['binDistance'],INT,yerr=INT_err,**kwargs)

        ax.set_xticks(Data['binDistance'].iloc[xvalues])
        ax.set_xticklabels(my_xticks, multialignment="center",ha="center")

        def calculateIndex(binDistance,x):
            idx = np.argmin(np.abs(binDistance-x))
            return idx

        ax.calculateIndex = lambda x: calculateIndex(Data['binDistance'],x)

        if rlu==False:
            ax.set_xlabel(r'$Q_x [\AA^{-1}]$'+'\n'+r'$Q_y [\AA^{-1}]$'+'\n'+'E/meV')
            def format_coord(x,y,ax,binCenter):# pragma: no cover
                index = ax.calculateIndex(x)
                qx,qy,E = binCenter[index]
                return  "qx = {0:.3e}, qy = {1:.3e}, E = {2:.3f}, I = {3:0.4e}".format(qx,qy,E,y)
        else:
            def format_coord(x,y,ax,binCenter):# pragma: no cover
                index = ax.calculateIndex(x)
                h,k,l,E = binCenter[index]
                return  "H = {0:.3e}, K = {1:.3e}, L = {2:.3e}, E = {3:.3f}, I = {4:0.4e}".format(h,k,l,E,y)
            ax.set_xlabel('$Q_h$ [RLU]\n$Q_k$ [RLU]\n$Q_l$ [RLU]\nE [meV]')

        
        ax.xaxis.set_label_coords(1.15, -0.025)
        ax.set_ylabel('$I$ [arb.u.]')
        plt.tight_layout()


        ax.format_coord = lambda x,y: format_coord(x,y,ax,np.array(Data[variables]))

        if not ufit:
            return ax,Data,bins
        
        # Create meta data for uFit dataset
        meta = dict()
        
        meta['instrument'] = self[0].instrument
        meta['experiment'] = ', '.join(d.experimentIdentifier for d in self)
        meta['title'] = self[0].title # TODO: Should be a collection of titles for all files?
        meta['datafilename'] = ', '.join(d.name for d in self)
        
        dist,Int = np.array(Data[['Energy','Int']]).T
        err = np.sqrt(Data['Intensity'])*Data['BinCount']/(Data['Monitor']*Data['Normalization'])
        data = np.array([dist,Int,err]).T
        xcol = 'E [meV]'
        ycol = 'Intensity'
        name = 'Intensity'
        ufitData = Dataset(meta=meta,data=data,xcol=xcol,ycol=ycol,name=name)
        
        return ax,ufitData


    def cutELine(self, Q1, Q2, Emin=None, Emax=None, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=True, dataFiles=None, constantBins=False):
        """Perform cut along energy in steps between two Q Point 
        
        Args:

            - Q1 (3D or 2D vector): Starting point for energy cut.

            - Q2 (3D or 2D vector): End of energy cut
            

        Kwargs: 

            - Emin (float): Start energy (default is self.Energy.min() for data in cut).
            
            - Emax (float): End energy (default is self.Energy.max() for data in cut).

            - energyWidth (float): Height of energy bins (default 0.05 meV)

            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel (default 0.1).

            - width (float): Full width of cut in q-plane (default 0.02).

            - rlu (bool): If True, provided Q point is interpreted as (h,k,l) otherwise as (qx,qy), (Default true)
            
            - dataFiles (list): Data files to be used. If none provided use the ones in self (default None)

            - constantBins (bool): If True only bins of size minPixel is used (default False)
            
        Returns:
            
            - Data (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for 1D cut.
            
            - Bin (1 array): Bin edge positions in energy

        """


        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                I = self.I.extractData()
                qx = self.qx.extractData()
                qy = self.qy.extractData()
                energy = self.energy.extractData()
                Norm = self.Norm.extractData()
                Monitor = self.Monitor.extractData()
                samples = self.sample
                maskIndices = self.maskIndices
                DS = self

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor,samples,maskIndices = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData(),DS.sample,DS.maskIndices
            



        Q1 = np.asarray(Q1,dtype=float)
        Q2 = np.asarray(Q2,dtype=float)

        dirvec = Q2-Q1
        
        # Copy the original mask to be reapplied later
        originalMask = DS.mask

        # Cut out relevant part of data set
        Q1re = Q1.copy().reshape(-1,1,1)
        
        if rlu:
            normal = DS.sample[0].planeNormal
            normal*=1.0/np.linalg.norm(normal)
            perp = np.cross(dirvec,normal)/np.linalg.norm(dirvec)

            mask = [np.array([np.logical_or(np.abs(np.einsum('i...,i->...',(A-Q1re),dirvec)-0.5)>0.5,np.abs(np.einsum('i...,i->...',A-Q1re,perp))>width) for A in zip(d.h,d.k,d.l)]) for d in DS]
        else:
            perp = dirvec[[1,0]]

            mask = [np.array([np.logical_or(np.abs(np.einsum('i...,i->...',(A-Q1re),dirvec)-0.5)>0.5,np.abs(np.einsum('i...,i->...',A-Q1re,perp))>width) for A in zip(d.qx,d.qy)]) for d in DS]

        # Combine the old mask and new with points close to cut
        DS.mask = [np.logical_or(mNew,mOld) for mNew,mOld in zip(mask,self.mask)]

        # Number of points along the Q cutting direction to reach minPixel
        points = np.linalg.norm(dirvec)/minPixel

        # Find minimum and maximum for newly masked data
        if Emin is None:
            Emin = DS.energy.min()
        if Emax is None:
            Emax = DS.energy.max()

        # Points for which constant Q cut in energy is to be performed
        QPoints = np.array([Q1+dirvec*x for x in np.linspace(0.0,1.0,int(np.floor(points)))])


        if rlu==True: # Recalculate H,K,L to qx
            rotationMatrices = [np.dot(samples[0].RotMat.T,s.RotMat) for s in samples]#[_tools.Rot(theta,deg=False) for theta in thetaDifference]
            Q = [[QX,QY] for QX,QY in zip(np.split(qx,maskIndices),np.split(qy,maskIndices))]
            qx,qy = np.concatenate([np.einsum('ij,j...->i...',rot,q) for rot,q in zip(rotationMatrices,Q)],axis=1)

            positions = np.array([qx,qy,energy])
            
            Qs = np.array([DS.convertToQxQy(q) for q in QPoints])
        else:
            positions = np.array([qx,qy,energy])
            Qs = QPoints

        
        Data = []
        Bins = []

        # Perform actual binning
        for i,Q in enumerate(Qs):
            Q = Q.flatten()
            [intensity,MonitorCount,Normalization,normcounts],bins  = cut1DE(positions = positions, I=I, Norm=Norm,Monitor=Monitor,E1=Emin,E2=Emax,q=Q,width=width,minPixel=energyWidth,constantBins=constantBins)
            data = pd.DataFrame()
            
            HKL = self.convertToHKL(Q.flatten())
            data['Qx'] = Q[0]*np.ones_like(intensity)
            data['Qy'] = Q[1]*np.ones_like(intensity)
            data['H'] = HKL[0]*np.ones_like(intensity)
            data['K'] = HKL[1]*np.ones_like(intensity)
            data['L'] = HKL[2]*np.ones_like(intensity)
            data['Energy'] = 0.5*(bins[0][1:]+bins[0][:-1])
            data['Intensity'] = intensity.astype(int)
            data['Monitor'] = MonitorCount.astype(int)
            data['Normalization'] = Normalization.astype(int)
            data['BinCount'] = normcounts.astype(int)
            data['QCut'] = i*np.ones_like(intensity).astype(int)
            
            data['Int'] = data['Intensity']*data['BinCount']/(data['Normalization']*data['Monitor'])
            Data.append(data)
            Bins.append(bins)

        Data = pd.concat(Data)
        DS.mask = originalMask

        return Data,Bins


    def plotCutELine(self, Q1, Q2, ax=None, Emin=None, Emax=None, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=True, dataFiles=None, constantBins=False, Vmin=None, Vmax = None, **kwargs):
        """Perform cut along energy in steps between two Q Point 
        
        Args:

            - Q1 (3D or 2D vector): Starting point for energy cut.

            - Q2 (3D or 2D vector): End of energy cut
            

        Kwargs: 

            - ax (matplotlib axis): Axis into which the plot is to go (default None, new created)

            - Emin (float): Start energy (default is self.Energy.min() for data in cut).
            
            - Emax (float): End energy (default is self.Energy.max() for data in cut).

            - energyWidth (float): Height of energy bins (default 0.05 meV)

            - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel (default 0.1).

            - width (float): Full width of cut in q-plane (default 0.02).

            - rlu (bool): If True, provided Q point is interpreted as (h,k,l) otherwise as (qx,qy), (Default true)
            
            - dataFiles (list): Data files to be used. If none provided use the ones in self (default None)

            - constantBins (bool): If True only bins of size minPixel is used (default False)

            - Vmin (float): Lower limit for colorbar (default min(Intensity)).
            
            - Vmax (float): Upper limit for colorbar (default max(Intensity)).
            
        Returns:
            
            - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for 1D cut.
            
            - Bin list (1 array): Bin edge positions in energy

        """
        
        Data,Bins = self.cutELine(Q1=Q1, Q2=Q2,Emin=Emin, Emax=Emax, energyWidth=energyWidth, minPixel =minPixel, width = width, rlu=rlu, dataFiles=dataFiles, constantBins=constantBins)

        
        if Vmin is None:
            Vmin = Data['Int'].min()
        if Vmax is None:
            Vmax = Data['Int'].max()

        dirvec = np.asarray(Q2)-np.asarray(Q1)


        if np.sign(dirvec[np.argmax(np.abs(dirvec))])==-1:
            dirvec = -dirvec
            
        if rlu:
            QPointColumns = ['H','K','L']
            visualizationBinPosition = np.array([-1,1])*0.5*np.linalg.norm(self.convertToQxQy(dirvec))/(len(Bins)-1)
        else:
            QPointColumns = ['Qx','Qy']
            visualizationBinPosition = np.array([-1,1])*0.5*np.linalg.norm(dirvec)/(len(Bins)-1)

        meshs = []

        if ax is None:
            if rlu:
                ortho = np.cross(self.sample[0].planeNormal,dirvec)
                ax = RLUAxes.createQEAxes(self,projectionVector1 = dirvec, projectionVector2 = ortho)
            else:
                fig, ax = plt.subplots()
        else:
            fig, ax = plt.subplots()
                
        columns = QPointColumns+['Int','QCut']
        for (_,_data),_bins in zip(Data[columns].groupby('QCut'),Bins):
            Q = np.array(_data[QPointColumns].iloc[0])
            if rlu:
                position = visualizationBinPosition+np.dot(Q,dirvec)/np.linalg.norm(dirvec)
            else:
                position = visualizationBinPosition+np.linalg.norm(Q)
                
            x,y = np.meshgrid(position,_bins[0])
            
            pmesh = ax.pcolormesh(x,y,_data['Int'].values.reshape(-1,1))
            meshs.append(pmesh)
            
        ax.meshs = meshs

        def set_clim(vmin,vmax,meshs):
            for p in meshs:
                p.set_clim(vmin,vmax)
                
        ax.set_clim = lambda Vmin,Vmax: set_clim(Vmin,Vmax,ax.meshs)

        ax.set_clim(Vmin,Vmax)

        return ax, Data, Bins


    @_tools.KwargChecker(function=createRLUAxes)
    def View3D(self,dQx,dQy,dE,rlu=True, log=False,grid=False,axis=2,counts=False,adjustable=True,customSlicer=False,outputFunction=print,**kwargs):
        """View data in the Viewer3D object. 

        Args:

            - dQx (float): step size in qx

            - dQy (float): step size in qy

            - dE (float): step size in E

        Kwargs:

            - rlu (Bool): If true a rlu axis is used for plotting otherwise qx,qy (Default True).

            - log (Bool): If true logarithm of intensity is plotted

            - grid (Bool): If true, grid is plotted. If float or integer, value is used as zorder of grid (Default False)

            - axis (int): Axis shown initially (default 2)

            - counts (bool): If set true, data shown is number of neutrons/pixel

            - adjustable (bool): If set true, 2 sliders will be present allowing to fine tune the c-axis (Default True)

            - customSlicer (bool): If true, utilize the interactive viewer based on PyQtGraph

            - kwargs: The remaining kwargs are given to the createRLUAxes method, intended for tick mark positioning (see createRLUAxes)

        If one plots not using RLU, everything is plotted in real units (1/AA), and the Qx and QY is not rotated. That is, the
        x axis in energy is not along the projection vector. The cuts of constant Qx and Qy does not represent any symmetry directions in
        the sample.
        However, if one utilizes the RLU flag, first Qx and Qy are rotated with first HKL vector along the x-axis. This thus means that 
        cuts of constant Qx (or more correctly along the principal HKL vector) represents a symmetry direction. However, as the data is 
        binned in equal sized voxels, constant Qy does not necessarily correspond to HKL vector 2 (it will in systems with 90 degrees 
        between the two vectors). 
        """

        if rlu:
            rluax = self.createRLUAxes(**kwargs)
            figure = rluax.get_figure()
            figure.delaxes(rluax)
            qxEax = self.createQEAxes(axis=1,figure=figure)
            figure.delaxes(qxEax)
            qyEax = self.createQEAxes(axis=0,figure=figure)
            figure.delaxes(qyEax)
            
            axes = [qxEax,qyEax,rluax]

        else:
            axes = None

        Data,bins = self.binData3D(dx=dQx,dy=dQy,dz=dE,rlu=rlu)
        
        if counts:
            Intensity = Data[0]/Data[3]
            Data = Intensity
            
        if customSlicer == True:
            
            if QtWidgets.QApplication.instance() is None:
                _cache.append(QtWidgets.QApplication(sys.argv))
            win = QtGui.QMainWindow()
            win.resize(800,800)
            win.setWindowTitle('Interactive View3D')
            win.setAttribute(QtCore.Qt.WA_DeleteOnClose)
            win.destroyed.connect(lambda: _cache.remove(win))
            _cache.append(win)
            
            Viewer = Viewer3DPyQtGraph.Interactive3DViewer(Data,bins,self.sample[0],log=log)
            win.setCentralWidget(Viewer)
            win.show()
            
            

        else:
            Viewer = Viewer3D.Viewer3D(Data=Data,bins=bins,axis=axis,ax=axes,grid=grid,log=log,adjustable=adjustable,outputFunction=outputFunction)
        return Viewer

    def cutRaw1D(self,detectorSelection=None,analyzerSelection=None):
        """cut 1D data to be used for 1D raw plot using specified DASEL.

        kwargs:
            detectorSelection (int): Detector to be used (default from file)
            
            analyzerSelection (int): Analyzer segment to be used (default from file)

        returns:

        """
        if len(self)>1:
            if not np.all([self[-1].scanParameters == d.scanParameters for d in self[:-1]]):
                raise AttributeError('Provided data files do not have the same scan variables!')
        
        # If no detector or analyzer is provided use dasel from file
        if analyzerSelection is None:
            analyzerSelection = [d.analyzerSelection for d in self]
        elif len(analyzerSelection) != len(self):
            raise AttributeError('Provided analyzerSelection list does not match length of dataset. Expected {}, but got {}'.format(len(self),len(analyzerSelection)))
            
        
        if detectorSelection is None:
            detectorSelection = [d.detectorSelection for d in self]
        elif len(detectorSelection) != len(self):
            raise AttributeError('Provided detectorSelection list does not match length of dataset. Expected {}, but got {}'.format(len(self),len(detectorSelection)))
        
        
        dataFiles = [d.original_file if d.type=='nxs' else d for d in self]
        
        for d in dataFiles:
            d.loadBinning(1)

        X = [d.scanValues for d in dataFiles]
       
        analyzerPixels = [d.instrumentCalibrationEdges[analyzer+detector*8] for d,analyzer,detector in zip(dataFiles,analyzerSelection,detectorSelection)] # Find pixel limits for specified selection (8 Energies/detector)
        
        I = [np.sum(d.I[:,detector,pixels[0]:pixels[1]],axis=1) for d,detector,pixels in zip(dataFiles,detectorSelection,analyzerPixels)]
        
        # Names to be changed for readability
        nameChange = {'polar_angle':'A4',
                'rotation_angle':'A3'}
        parameter = dataFiles[0].scanParameters
        # If parameter is in nameChange, change it. Otherwise keep it
        parameter = [nameChange.get(p, p) for p in parameter]
        unit = dataFiles[0].scanUnits

        return X,I,parameter,unit
         

    def plotRaw1D(self,detectorSelection=None,analyzerSelection=None, legend=True, grid=-10):
        """plot 1D figures of data in the specified DASEL.
        
        kwargs:
            detectorSelection (int): Detector to be used (default from file)
            
            analyzerSelection (int): Analyzer segment to be used (default from file)

            legend (bool list): Insert legend with provided names or file names (default True)

            grid (bool, Int): If true plot at grid on figure. If integer set zorder of grid (default -10)
            
        returns:
            Ax (list): List of axes in which data are plotted.
        """

        
        def intrextrapolate(oldPosition,oldValues,newValues,outputFunction=print):# pragma: no cover
            """interpolates between old and new through linear regression and returns best estimate at old positions
            
            arg:
                - oldPosition (list): List of position to estimate the new value at (in coordinates of oldValues)
            
                - oldValues (list): List of old values 
            
                - newValues (list): List of new values

            kwargs:

                - outputFunction (function): Function called on output string (default print)

            return:
                newPosition (list): estimate of newValue at oldPosition
            """
                
            oldOffset = oldValues[0]
            oldSize = np.diff([oldValues[0],oldValues[-1]])
            
            newOffset = newValues[0]
            newSize = np.diff([newValues[0],newValues[-1]])
            
            X = (oldPosition-oldOffset)/oldSize  # Linear interpolation
        
            newPosition = X*newSize+newOffset # Convert to new values
            return newPosition
        
        
        
        def format_coord(x,y,X,labels):# pragma: no cover
            """format coordinates according to x and y. When X is multidimensional an approximation for its value at x is estimated"""
            fmtOrder = np.ceil(-np.log(np.mean(np.abs(np.diff(X,axis=1)),axis=1)))
            fmtOrder = [o if o>0 else 0 for o in (fmtOrder).astype(int)]
            
            xFormatString = ['{:.'+str(fO)+'f}' for fO in fmtOrder] # Format to use enough decimals to show all changes
            
            # Calculate estimates for values of X[i] at x
            newX = [intrextrapolate(x,X[0],XX)[0] for XX in X[1:]]
             
            xs = [fstring.format(XX) for fstring,XX in zip(xFormatString,np.concatenate([[x],newX],axis=0))]
            return ', '.join([label+' = '+str(X) for X,label in zip(np.concatenate([xs,[y]],axis=0),labels)])
        
        def onclick(event,ax,outputFunction):# pragma: no cover
            if ax.in_axes(event):
                try:
                    C = ax.get_figure().canvas.cursor().shape() # Only works for pyQt5 backend
                except:
                    pass
                else:
                    if C != 0: # Cursor corresponds to arrow
                        return
        
                x = event.xdata
                y = event.ydata
                if hasattr(ax,'__format_coord__'):
                    outputFunction(ax.__format_coord__(x,y))
                else:
                    outputFunction(ax.format_coord(x,y))
        
    
        
        
        _,ax = plt.subplots()
        
        
        if legend:
            if legend is True:
                legend = [d.name for d in self]
            else:
                if len(legend) != len(self):
                    raise AttributeError('Provided list of legends does not match length of dataset. Expected {}, but got {}'.format(len(self),len(legend)))

        ax.X,ax.I,ax.parameter,ax.unit = self.cutRaw1D(detectorSelection=detectorSelection,analyzerSelection=analyzerSelection)

        ax.xlabels = ['{} [{}s]'.format(p,u) for p,u in zip(ax.parameter,ax.unit)]
        ax.__labels__ = np.concatenate([ax.xlabels,['Int [count]']],axis=0)
        
        plots = []
        for X,I in zip(ax.X,ax.I):
            plots.append(ax.scatter(X[0],I))
            
        
        ax.set_ylabel('Int [counts]')
        
        
        ax.__format_coord__ = lambda x,y: format_coord(x,y,X=np.concatenate(ax.X,axis=0),labels=ax.__labels__)
        ax.format_coord = lambda x,y: ax.__format_coord__(x,y)
        ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,outputFunction=outputFunction))
        #        ax.format_xdata = lambda x: format_xdata(x,ax)#ax.X,ax.xlabels)
        ax.set_xlabel(ax.xlabels[0])
        if legend:
            for plot,lab in zip(plots,legend):
                plot.set_label(lab)
            ax.legend()
        
        ax.grid(grid)
        return ax

    def convertToQxQy(self,HKL):
        """Convert array or vector of HKL point(s) to corresponding Qx and QY

        Args:

            - HKL (array): array or vector of HKL point(s)

        Returns

            - Q (array): Converted HKL points in Qx QY of un-rotated coordinate system.
        """
        return convertToQxQy(self.sample[0],HKL)

    def convertToHKL(self,QxQy):
        """Convert array or vector of QxQy point(s) to corresponding HKL

        Args:

            - QxQy (array): array or vector of QxQy point(s)

        Returns

            - HKL (array): Converted QxQy points in HKL
        """
        return convertToHKL(self.sample[0],QxQy)

    def updateCalibration(self,calibFiles,overwrite=False):
        """Update calibrations for all data files in data set. Does not save the changes.
        
        Args:
            
            - calibrationFile (string or list): calibration file, as generated from MJOLNIR.Geometry.Instrument or list of these.
            
        Kwargs:
            
            - overwrite (bool): If true, previous binnings will be overwritten if new files contain same binning (default False)
        
    .. note::
        Changes performed by this method is not saved to disk.     
            
        """
        for d in self:
            d.updateCalibration(calibFiles,overwrite=overwrite)


    def generateUFitDataset(self, pdData,q1,q2,rlu,width,minPixel,Emin,Emax,QDirection=True):
        """Generate uFitDataset from cut.

        Args:

            - pdData (pandas dataframe): Data generated from 1D cut

            - q1 (array): Start point for cut

            - q2 (array): End point for cut

            - rlu (bool): If in reciprocal lattice unites or not

            - width (float): Width size (used for rounding of labels)

            - minPixel (float): Minimum pixel size (used for rounding of labels)

            - Emin (float): Minimum energy
             
            - Emax (float): Maximum energy

        Kwargs:

            - QDirection (bool): If true ufitdata is created along Q, otherwise energy (default True)

        """

        if rlu:
            variables = ['H','K','L']
        else:
            variables = ['Qx','Qy']
        variables = variables

        if QDirection:
            QRounding = int(-np.round(np.log10(minPixel)))
            ERounding = int(np.round(6/(np.linalg.norm(Emin-Emax))))
            ERounding = np.max([ERounding,1])
            QRounding = np.max([QRounding,1])

            dirVec = np.array(q2)-np.array(q1)
            if rlu:
                dirVec = _tools.LengthOrder(dirVec)
            else:
                dirVec = _tools.Norm2D(dirVec)

            pdData['binDistance'] = np.linalg.norm(pdData[variables]-np.array(pdData[variables].iloc[1]),axis=1)/np.linalg.norm(dirVec)

            startPos = pdData[variables].iloc[0]   

            if rlu:
                xdirection = _tools.generateLabel(np.round(dirVec,QRounding))[1:-1].split(', ')
            else:
                xdirection = _tools.generateLabel(np.round(dirVec,QRounding),labels=['Qx','Qy'])[1:-1].split(', ')
            xConstantOffset = np.round(startPos,QRounding)
            xlabel = []
            for l,val in zip(xdirection,xConstantOffset):
                if np.isclose(val,0):
                    xlabel.append(l)
                elif l=='0':
                    xlabel.append('{}'.format(val))
                else:
                    xlabel.append(l+'{:+}'.format(val))

        else:
            ERounding = int(-np.round(np.log10(minPixel)))
            QRounding = int(-np.round(np.log10(width)))

            ERounding = np.max([ERounding,1])
            QRounding = np.max([QRounding,1])
            pdData['binDistance'] = pdData['Energy']


            xlabel = [str(x) for x in np.array(q1,dtype=float).round(QRounding)]

        
        x = np.array(pdData['binDistance'])

        # Calcualte mean energy from bins (last return value)
        Energy = (Emin+Emax)*0.5
        # Create meta data for uFit dataset
        meta = dict()

        meta['instrument'] = self[0].instrument
        meta['experiment'] = ', '.join(d.experimentIdentifier for d in self)
        meta['title'] = self[0].title # TODO: Should be a collection of titles for all files?
        meta['datafilename'] = ', '.join(d.name for d in self)

        Int = np.array(pdData['Int'])
        err = np.sqrt(pdData['Intensity'])*pdData['BinCount']/(pdData['Monitor']*pdData['Normalization'])
        data = np.array([x,Int,err]).T

        

        xcol = '\n'.join(xlabel)+'\n'+'{:.3}'.format(np.round(Energy,ERounding))
        if rlu:
            xcol+='\n[RLU,meV]'
        else:
            xcol+='\n'+r'[$\AA^{-1}$,meV]'
        ycol = 'Intensity'
        name = 'Intensity'
        ufitData = Dataset(meta=meta,data=data,xcol=xcol,ycol=ycol,name=name)
        return ufitData

    def updateSampleParameters(self,unitCell):
        """Update unit cell parameters and corresponding UB matrix

        Args:

            - unitCell (list): List of cell parameters (a,b,c,alpha,beta,gamma)

        """
        for d in self:
            d.sample.updateSampleParameters(unitCell=unitCell)


    def __sub__(self,other):

        filesSelf = [d for d in self]
        filesOther = [d for d in other]
        if not np.all([x.type == 'nxs' for x in np.concatenate([filesSelf,filesOther])]):
            raise AttributeError('Data files have to be converted!')
        monoSelf = [x.MonitorPreset for x in self]
        monoOther = [x.MonitorPreset for x in other]

        data = []
        for s,o in zip(filesSelf,filesOther):
            sMono = s.MonitorPreset
            oMono = o.MonitorPreset
            if sMono>oMono:
                s.I = s.I-o.I*(oMono/sMono)
            elif sMono<oMono:
                s.I = s.I*(sMono/oMono)-o.I
            else:
                s.I = s.I-o.I
            data.append(s)

        newFile = DataSet(data)
        return newFile

    def undoAbsolutNormalize(self):
        """Undo normalization previously performed"""
        if self.absolutNormalized!=0:
            normFactor = self.absolutNormalized
            for d in self:
                d.Norm *= normFactor
                if hasattr(d,'absolutNormalizationFactor'):
                    d.absolutNormalizationFactor*= 1.0/normFactor
                else:
                    d.absolutNormalizationFactor= 1
                d.absolutNormalized = False

            self.absolutNormalized = False




    def absolutNormalize(self,sampleMass=None,sampleMolarMass=None,sampleChemicalFormula=None,formulaUnitsPerUnitCell=1.0,
                         sampleGFactor=2.0, correctVanadium=False,vanadiumMass=15.25,
                         vanadiumMonitor=100000,vanadiumSigmaIncoherent=5.08,vanadiumChemicalFormula='V',vanadiumGFactor=2.0,
                         vanadiumUnitsPerUnitCell=1.0,vanadiumMolarMass=None):
        """Normaliza dataset to absolut units () by 

        Kwargs:

            - sampleMass (float): Mass of sample in gram

            - sampleChemicalFormula (string): Chemical formula of sample

            - formulaUnitsPerUnitCell (float): Number of formula units per unit cell (default 1.0)

            - sampleDebyeWaller (float): Debye Waller factor (default 1.0)

            - sampleGFactor (float): Magnetic G factor for sample (defalt 2.0)
        
            - sampleFormFactor (float): Formfactor of sample (default 1.0)

            - correctVanadium (bool): If normalization files have not been updated set this to True (default False)

            - vanadiumMass (float): Mass of vanadium used in normalization in gram (default 15.25)

            - vanadiumMonitor (int): Monitor count used in normalization scan (default 100000)

            - vanadiumSigmaIncoherent (float): Incoherent scattering strength of Vanadium (default 5.08)

        """

        if len(self.convertedFiles) == 0:
            raise AttributeError("Data set needs to be converted before absolut normalization can be applied.")
        
        

        normFactor = \
        _tools.calculateAbsolutNormalization(sampleChemicalFormula=sampleChemicalFormula,sampleMolarMass=sampleMolarMass,sampleMass=sampleMass,
                                             formulaUnitsPerUnitCell=formulaUnitsPerUnitCell,sampleGFactor=sampleGFactor,
                                             correctVanadium=correctVanadium,vanadiumMass=vanadiumMass,vanadiumChemicalFormula=vanadiumChemicalFormula,
                                             vanadiumMonitor=vanadiumMonitor,vanadiumSigmaIncoherent=vanadiumSigmaIncoherent,vanadiumMolarMass=vanadiumMolarMass,
                                             vanadiumGFactor=vanadiumGFactor,vanadiumUnitsPerUnitCell=vanadiumUnitsPerUnitCell)
            
        if self.absolutNormalized != 0:
            warnings.warn("\nAlready Normalized\nDataSet seems to already have beeen normalized absolutly. Reverting previous normalization...")
            normFactor /= self.absolutNormalized

        for d in self:
            d.Norm *= normFactor
            if hasattr(d,'absolutNormalizationFactor'):
                d.absolutNormalizationFactor*= normFactor
            else:
                d.absolutNormalizationFactor= normFactor
            d.absolutNormalized = True

        if self.absolutNormalized != 0:
            self.absolutNormalized *= normFactor
        else:
            self.absolutNormalized = normFactor




def load(filename):
    """Function to load an object from a pickled file.

    .. note::
        It is not possible to un-pickle an object created in python 3 in python 2 or vice versa.
        
    """
    try:                                # Opening the given file with an error catch
        fileObject = open(filename, 'rb')
    except IOError as e: # pragma: no cover
        print("Error in opening file:\n{}".format(e))
    else:
        tmp_dict = pickle.load(fileObject)
        
        fileObject.close()
        return tmp_dict

@_tools.KwargChecker()
def cut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False,extend=True,constantBins=False):
    """Perform 1D cut through constant energy plane from q1 to q2 returning binned intensity, monitor, normalization and normcount. The full width of the line is width while height is given by Emin and Emax. 
    the minimum step sizes is given by minPixel.
    
    .. note::
        Can only perform cuts for a constant energy plane of definable width.
    
    Args:
        
        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.
        
        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - q1 (2D array): Start position of cut in format (qx,qy).
        
        - q2 (2D array): End position of cut in format (qx,qy).
        
        - width (float): Full width of cut in q-plane.
        
        - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.
        
        - Emin (float): Minimal energy to include in cut.
        
        - Emax (float): Maximal energy to include in cut
        
    Kwargs:
        
        - plotCoverage (bool): If True, generates plot of all points in the cutting plane and adds bounding box of cut in new figure, if axis, plots on top (default False).

        - extend (bool): Whether or not the cut from q1 to q2 is to be extended throughout the data (default true)

        - constantBins (bool): If True only bins of size minPixel is used (default False)
    
    Returns:
        
        - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
        
        - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
    """
    dirvec = np.array(q2,dtype=float)-np.array(q1,dtype=float)
    dirLength = np.linalg.norm(dirvec)
    dirvec/=dirLength
    orthovec=np.array([dirvec[1],-dirvec[0]])
    
    ProjectMatrix = np.array([dirvec,orthovec])
    insideEnergy = np.logical_and(positions[2]<=Emax,positions[2]>=Emin)
    if(np.sum(insideEnergy)==0):
        #raise AttributeError('No points are within the provided energy limits.')
        return [np.array(np.array([])),np.array([]),np.array([]),np.array([])],[np.array([]),np.array([]),[Emin,Emax]]
        

    positions2D = np.array([positions[0][insideEnergy], positions[1][insideEnergy]])
    propos = np.dot(ProjectMatrix,positions2D-q1.reshape(2,1))
    
    if extend==False: # Only take points between the given q points
        if constantBins==False:
            insideQ = np.logical_and(propos[0]>-0.05,propos[0]<dirLength*1.05)
        else:
            insideQ = np.logical_and(propos[0]>0.0,propos[0]<dirLength)
        propos = propos[:,insideQ]


    orthobins = [-width/2.0,width/2.0]
    insideWidth = np.logical_and(propos[1]<orthobins[1],propos[1]>orthobins[0])
    
    if constantBins==False:
        lenbins = np.array(_tools.binEdges(propos[0][insideWidth],minPixel,startPoint=0.0,endPoint=dirLength))
    else:
        Min,Max = _tools.minMax(propos[0][insideWidth])
        lenbins = np.arange(Min,Max+0.5*minPixel,minPixel)
    orthopos = np.outer(orthobins,orthovec)
    binpositions = np.outer(lenbins,dirvec)+q1
    
    if len(lenbins)==0:
        return [np.array(np.array([])),np.array([]),np.array([]),np.array([])],[np.array([]),orthopos,[Emin,Emax]]
    
    normcounts = np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=np.ones((propos.shape[1])).flatten())[0]

    if extend==False: # Test both inside energy range AND inside q-limits
        intensity = np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=I[insideEnergy][insideQ].flatten())[0]
        MonitorCount=  np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=Monitor[insideEnergy][insideQ].flatten())[0]
        Normalization= np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=Norm[insideEnergy][insideQ].flatten())[0]
    else:
        intensity = np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=I[insideEnergy].flatten())[0]
        MonitorCount=  np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=Monitor[insideEnergy].flatten())[0]
        Normalization= np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=Norm[insideEnergy].flatten())[0]
    
    EmeanVec = np.ones((len(binpositions),1))*(Emin+Emax)*0.5
    binpositionsTotal = np.concatenate((binpositions,EmeanVec),axis=1)
   
    if not plotCoverage is False: # pragma: no cover
        if not isinstance(plotCoverage,bool): # Assuming matplotlib axis
            ax = plotCoverage
            plotPoints = False
        else:
            fig,ax = plt.subplots() # Generate new figure
            plotPoints = True
            ax.scatter(positions2D[0],positions2D[1],s=0.5,zorder=100)
        ax.plot([binpositions[0][0]+orthopos[0][0],binpositions[-1][0]+orthopos[0][0]],[binpositions[0][1]+orthopos[0][1],binpositions[-1][1]+orthopos[0][1]],c='w',zorder=100)
        ax.plot([binpositions[0][0]+orthopos[1][0],binpositions[-1][0]+orthopos[1][0]],[binpositions[0][1]+orthopos[1][1],binpositions[-1][1]+orthopos[1][1]],c='w',zorder=100)
        for i in [0,-1]:
            ax.plot([binpositions[i][0]+orthopos[0][0],binpositions[i][0]+orthopos[1][0]],[binpositions[i][1]+orthopos[0][1],binpositions[i][1]+orthopos[1][1]],c='w',zorder=100)
        for binPos in binpositions:#i in range(len(binpositions)):
            ax.plot([binPos[0]+orthopos[0][0],binPos[0]+orthopos[1][0]],[binPos[1]+orthopos[0][1],binPos[1]+orthopos[1][1]],c='w',linewidth=0.5,zorder=100)
        if extend==False and plotPoints:
            ax.scatter(positions2D[0][insideQ][insideWidth],positions2D[1][insideQ][insideWidth],s=0.5,zorder=100)
        elif plotPoints:
            ax.scatter(positions2D[0][insideWidth],positions2D[1][insideWidth],s=0.5,zorder=100)
        if plotPoints:
            ax.set_aspect('equal', 'datalim')
            ax.set_xlabel(r'Qx [$\AA^{-1}$]')
            ax.set_ylabel(r'Qy [$\AA^{-1}$]')

    return [intensity,MonitorCount,Normalization,normcounts],[binpositionsTotal,orthopos,np.array([Emin,Emax])]



def cut1DE(positions,I,Norm,Monitor,E1,E2,q,width,minPixel,constantBins=False):#,plotCoverage=False):
    """Perform 1D cut through constant Q point returning binned intensity, monitor, normalization and normcount. The width of the cut is given by 
    the width attribute. 
    
    .. note::
        Can only perform cuts for a constant energy plane of definable width.
    
    Args:
        
        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.
        
        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - E1 (float): Start energy.
        
        - E2 (float): End energy.

        - q (2d vector): Q point in (qx,qy)
        
        - width (float): Full width of cut in q-plane.
        
        - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.
        
        - Emin (float): Minimal energy to include in cut.
        
        - Emax (float): Maximal energy to include in cut

    Kwargs:

        - constantBins (bool): If True only bins of size minPixel is used (default False)
        
    Returns:
        
        - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
        
        - Bin list (1 array): Bin edge positions in energy
        
    """
    if len(q.shape)==1:
        q.shape = (2,1)
    distToQ = np.linalg.norm(positions[:2]-q,axis=0)

    inside = distToQ<width
    
   
    
    insideEnergy = np.logical_and(positions[2]<=E2,positions[2]>=E1)
    if(np.sum(insideEnergy)==0):
        raise AttributeError('No points are within the provided energy limits.')
    elif(np.sum(inside)==0):
        raise AttributeError('No points are inside selected q range.')

    allInside = np.logical_and(inside,insideEnergy)
    Energies = positions[2][allInside]
    
    if constantBins==False:
        bins = np.array(_tools.binEdges(Energies,tolerance=minPixel))
    else:
        Min,Max = _tools.minMax(Energies)
        bins = np.arange(Min,Max+0.5*minPixel,minPixel)
    
    if len(bins)==0:
        return [np.array(np.array([])),np.array([]),np.array([]),np.array([])],[[E1,E2]]
    
    normcounts = np.histogram(Energies,bins=bins,weights=np.ones_like(Energies).flatten())[0]
    intensity = np.histogram(Energies,bins=bins,weights=I[allInside].flatten())[0]
    MonitorCount=  np.histogram(Energies,bins=bins,weights=np.array(Monitor[allInside].flatten(),dtype=np.float))[0] # Need to change to int64 to avoid overflow
    Normalization= np.histogram(Energies,bins=bins,weights=np.array(Norm[allInside].flatten(),dtype=np.float))[0]
    

    return [intensity,MonitorCount,Normalization,normcounts],[bins]



@_tools.KwargChecker()
def cutPowder(positions,I,Norm,Monitor,EBinEdges,qMinBin=0.01,constantBins=False):
    """Cut data powder map with intensity as function of the length of q and energy. 

    Args:

        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.

        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - EBinEdges (list): Bin edges between which the cuts are performed.

    Kwargs:

        - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.

        - constantBins (bool): If True only bins of size minPixel is used (default False)

    Returns:
        
        - Data list (pandas DataFrame): DataFrame containing qx,qy,H,K,L,Intensity,Normalization,Monitor,BinCount,Int,binDistance for powder cut.
        
        - qbins (n arrays): n arrays holding the bin edges along the length of q

    """
    qx,qy,energy = positions

    q = np.linalg.norm([qx,qy],axis=0)

    qbins = []
    
    data = []
    #for i in range(len(EBinEdges)-1):
    for energyBin,[binStart,binEnd] in enumerate(zip(EBinEdges,EBinEdges[1:])):
        e_inside = np.logical_and(energy>binStart,energy<=binEnd)
        q_inside = q[e_inside]
        if constantBins==False:
            qbins.append(np.array(_tools.binEdges(q_inside,tolerance=qMinBin)))
        else:
            Min,Max = _tools.minMax(q_inside)
            qbins.append(np.arange(Min,Max+0.5*qMinBin,qMinBin))
            
        intensity,bins = np.histogram(q_inside,bins=qbins[-1],weights=I[e_inside].flatten())
        monitorCount = np.histogram(q_inside,bins=qbins[-1],weights=Monitor[e_inside].flatten())[0].astype(Monitor.dtype)
        Normalization = np.histogram(q_inside,bins=qbins[-1],weights=Norm[e_inside].flatten())[0].astype(Norm.dtype)
        NormCount = np.histogram(q_inside,bins=qbins[-1],weights=np.ones_like(I[e_inside]).flatten())[0].astype(I.dtype)

        _data = pd.DataFrame(np.array([intensity,monitorCount,Normalization,NormCount],dtype=np.int).T,
                        columns=['Intensity','Monitor','Normalization','BinCount'],dtype=np.int)
        _data['q'] = 0.5*(bins[:-1]+bins[1:])
        _data['Energy'] = np.ones_like(_data['q'])*0.5*(binStart+binEnd)
        _data['EnergyCut'] = np.ones_like(_data['q'],dtype=np.int)*energyBin
        data.append(_data)
    data = pd.concat(data)
    data['Int'] = data['Intensity']*data['BinCount']/(data['Normalization']*data['Monitor'])
    return data,qbins




#@_tools.KwargChecker(function=plt.pcolormesh)
#def plotQPlane(I,Monitor,Norm,pos,EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=False,ax=None,**kwargs):
#    """Plotting tool to show binned intensities in the Q plane between provided energies.
#    
#    Args:
#        
#        - I (array): Intensity of data.
#        
#        - Monitor (array): Monitor of data.
#        
#        - Norm (array): Normalization of data.
#        
#        - pos (3 array): Position of data in qx, qy, and energy.
#        
#        - EMin (float): Lower energy limit.
#        
#        - EMax (float): Upper energy limit.
#        
#    Kwargs:
#        
#        - binning (str): Binning scheme, either 'xy' or 'polar' (default 'xy').
#        
#        - xBinTolerance (float): bin sizes along x direction (default 0.05). If enlargen is true, this is the minimum bin size.
#
#        - yBinTolerance (float): bin sizes along y direction (default 0.05). If enlargen is true, this is the minimum bin size.
#        
#        - enlargen (bool): If the bin sizes should be adaptive (default False). If set true, bin tolerances are used as minimum bin sizes.
#
#        - log (bool): Plot intensities as the logarithm (default False).
#        
#        - ax (matplotlib axes): Axes in which the data is plotted (default None). If None, the function creates a new axes object.
#        
#        - other: Other key word arguments are passed to the pcolormesh plotting algorithm.
#        
#    Returns:
#        
#        - ax (matplotlib axes)
#        
#    .. note::
#        The axes object gets a new method denoted 'set_clim' taking two parameters (VMin and VMax) used to change axes coloring.
#        
#        
#    """
#    qx,qy,energy=pos
#
#
#    
#    if ax is None:
#    #        if RLUPlot:
#    #            ax = self.createRLUAxes()
#    #        else:
#        plt.figure()
#        ax = plt.gca()
#            
#    
#    
#    binnings = ['xy','polar']#,'rlu']
#    if not binning in binnings:
#        raise AttributeError('The provided binning is not understood, should be {}'.format(', '.join(binnings)))
#    if binning == 'polar':
#        x = np.arctan2(qy,qx)
#        y = np.linalg.norm([qx,qy],axis=0)  
#    
#    elif binning == 'xy':
#        x = qx
#        y = qy
#        
#    #elif binning == 'rlu':
#    #    raise NotImplementedError('Currently the RLU binning is not implimented')
#    
#    
#    EBinEdges = [EMin,EMax]
#
#    intensity = []
#    monitorCount = []
#    Normalization = []
#    NormCount = []
#    bins = []
#    
#    
#    for i in range(len(EBinEdges)-1):
#        e_inside = np.logical_and(energy>EBinEdges[i],energy<=EBinEdges[i+1])
#        if enlargen:
#            yBins = _tools.binEdges(y[e_inside],yBinTolerance)
#        else:
#            yBins = np.arange(np.min(y[e_inside]),np.max(y[e_inside]),yBinTolerance)
#        for j in range(len(yBins)-1):
#            ey_inside = np.logical_and(np.logical_and(e_inside,np.logical_and(y>yBins[j],y<yBins[j+1])),(1-np.isnan(Norm)).astype(bool))
#            
#            x_inside = x[ey_inside]
#            #y_inside = y[ey_inside]
#            
#            if enlargen:
#                xbins = _tools.binEdges(x_inside,tolerance=xBinTolerance)
#            else:
#                xbins = np.arange(np.min(x),np.max(x),xBinTolerance)
#                
#            if len(xbins)==0:
#                continue
#            bins.append(np.array([xbins,np.array([yBins[j],yBins[j+1]])]))
#            
#            intensity.append(np.histogram(x_inside,bins=bins[-1][0],weights=I[ey_inside].flatten())[0].astype(I.dtype))
#            monitorCount.append(np.histogram(x_inside,bins=bins[-1][0],weights=Monitor[ey_inside].flatten())[0].astype(Monitor.dtype))
#            Normalization.append(np.histogram(x_inside,bins=bins[-1][0],weights=Norm[ey_inside].flatten())[0].astype(Norm.dtype))
#            NormCount.append(np.histogram(x_inside,bins=bins[-1][0],weights=np.ones_like(I[ey_inside]).flatten())[0].astype(I.dtype))
#
#    warnings.simplefilter('ignore')
#    Int = [np.divide(intensity[i]*NormCount[i],monitorCount[i]*Normalization[i]) for i in range(len(intensity))]
#    warnings.simplefilter('once')
#    
#    if binning == 'polar':
#        Qx = [np.outer(bins[i][1],np.cos(bins[i][0])).T for i in range(len(intensity))]
#        Qy = [np.outer(bins[i][1],np.sin(bins[i][0])).T for i in range(len(intensity))]
#    
#    elif binning == 'xy':
#        Qx = [np.outer(bins[i][0],np.ones_like(bins[i][1])) for i in range(len(intensity))]
#        Qy = [np.outer(np.ones_like(bins[i][0]),bins[i][1]) for i in range(len(intensity))]
#        
#   
#    pmeshs = []
#    if log:
#        Int = [np.log10(1e-20+np.array(Int[i])) for i in range(len(Int))]
#    for i in range(len(intensity)):
#        pmeshs.append(ax.pcolormesh(Qx[i],Qy[i],Int[i].reshape((len(Int[i]),1)),zorder=10,**kwargs))
#    ax.set_aspect('equal', 'datalim')
#    ax.grid(True, zorder=0)
#    ax.set_clim = lambda VMin,VMax: [pm.set_clim(VMin,VMax) for pm in pmeshs]
#    ax.pmeshs = pmeshs
#    return ax

@_tools.KwargChecker()
def plotA3A4(files,ax=None,planes=[],binningDecimals=3,log=False,returnPatches=False,singleFigure=False,plotTessellation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2): # pragma: no cover
    """Plot data files together with pixels created around each point in A3-A4 space. Data is binned in the specified planes through their A3 and A4 values. 
    This can result in distorted binning when binning across large energy regions. Data is plotted using the pixels calculated for average plane value, i.e. 
    binning 7,8,9,10, and 11 patches for plane 9 are used for plotting.

    Args:
        
        - files (DataFiles): single file or list of files to be binned together

    Kwargs:

        - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)

        - planes (list (of lists)): Planes to be plotted and binned (default [])

        - binningDecimals (int): Number of decimal places A3-A4 positions are rounded before binning (default 3)
        
        - log (bool): Whether or not to plot intensities as logarithm (default False)

        - returnPatches (bool): If true the method returns the patches otherwise plotted in the given axis (default False)

        - singleFigure (bool): If true, all planes are plotted in same figure (default False)

        - plotTessellation (bool): Plot Tessellation of points (default False)

        - Ei_err (float): Tolerance of E_i for which the values are equal (default = 0.05)

        - temperature_err (float): Tolerance of temperature for which the values are equal (default = 0.2)
        
        - magneticField_err (float): Tolerance of magnetic field for which the values are equal (default = 0.2)
        
        - electricField_err (float): Tolerance of electric field for which the values are equal (default = 0.2)

    Returns:
        
        - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.

    Raises:

        - AttributeError

    Examples:

    The following example will combine the two files and plot all of the available planes in different figures.

    >>> DS = DataSet.DataSet(convertedFiles=[--.nxs,---.nxs])
    >>> plt.figure()
    >>> ax = plt.gca()
    >>>
    >>> DataSet.plotA3A4(DS.convertedFiles,ax=ax)

    If only a subset of planes or different planes are to be combined the following will achieve this:

    >>> DataSet.plotA3A4(DS.convertedFiles,ax=ax,planes=[0,1,2,3,[4,5,6],[8,9]])

    Here planes 0 through 3 are plotted separately while 4,5, and 6 as well as 8 and 9 are binned.

    .. note::
        Binning planes from different analysers might result in nonsensible binnings.

    """
   
    if not isinstance(ax, (list,)) and ax is not None:
        ax = np.array([ax])
    
    if not isinstance(planes, (list,)):
        planes = np.array([planes])
        
    if not ax is None:
        if singleFigure and np.array([ax]).size != 1:
            raise AttributeError('Single figure chosen but multiple axes given ({}).'.format(np.array([ax]).size))
        
        elif not singleFigure and len(ax) != len(planes) and not len(planes)==0:
            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,len(planes)))
            
    
    if not isinstance(files,(list,np.ndarray)):
        files = [files]
    numFiles = len(files)

    @_tools.my_timer_N()
    def testFiles(files,numFiles):
        
        if numFiles>1:
            comparison = np.array([np.all([np.isclose(files[0].Ei,files[i+1].Ei,atol=Ei_err) for i in range(numFiles-1)]),\
                    np.all([compareNones(files[0].temperature,files[i+1].temperature,temperature_err) for i in range(numFiles-1)]),\
                    np.all([compareNones(files[0].magneticField,files[i+1].magneticField,magneticField_err) for i in range(numFiles-1)]),\
                    np.all([compareNones(files[0].electricField,files[i+1].electricField,electricField_err) for i in range(numFiles-1)]),\
                    np.all([files[0].binning==files[i+1].binning for i in range(numFiles-1)])])
            
            tests = np.array(['Ei','Temperature','Magnetic Field','Electric Field','Binning'])
            
            if not np.all(comparison):
                errors = np.array(1-comparison,dtype=bool)
                raise AttributeError('Attributes for the datafiles are not the same! Difference is in :\n'+','.join([x for x in tests[errors]])+'\nIf the files are to be binned anyway change the tolerence limits.')
    testFiles(files,numFiles)

    #@_tools.my_timer_N()
    #def getA3A4(files,numFiles):
    #print(numFiles)
    #type(files)
    #for i in range(numFiles):
       # print(type(files[i]))

    #print('_____************_____')
    A4All = np.array([files[i].A4 for i in range(numFiles)])
    A3All = np.array([files[i].A3 for i in range(numFiles)])
    #    return A3All,A4All

    
    #A3All,A4All = getA3A4(files,numFiles)

    #@_tools.my_timer_N()
    #def getData(files,numFiles):
    Ishape = files[0].I.shape
    IAll = np.array([files[i].I for i in range(numFiles)]) # into shape sum(A3),104,64 for CAMEA ## np.array([files[i].I[:,0,0,:,:].reshape((A3All[i].size,Ishape[3],Ishape[4])) for i in range(numFiles)])
    NormAll = np.array([files[i].Norm for i in range(numFiles)]) ## np.array([files[i].Norm[:,0,0,:,:].reshape((A3All[i].size,Ishape[3],Ishape[4])) for i in range(numFiles)])
    MonitorAll = np.array([files[i].Monitor for i in range(numFiles)]) ## np.array([files[i].Monitor[:,0,0,:,:].reshape((A3All[i].size,Ishape[3],Ishape[4])) for i in range(numFiles)])
    
    if not ax is None:
        if not singleFigure and len(ax) != Ishape[2] and len(planes) == 0: # Plot all planes in provided axes
            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,Ishape[2]))

    #@_tools.my_timer_N()
    #def concatINormMon(IAll, NormAll,MonitorAll):
    I = np.concatenate(IAll,axis=0)
    Norm = np.concatenate(NormAll,axis=0)
    Mon = np.concatenate(MonitorAll,axis=0)
    #    return I,Norm,Mon
    #I,Norm,Mon = concatINormMon(IAll, NormAll,MonitorAll)


    #@_tools.my_timer_N()
    #def A4Instr(files,numFiles):
    A4InstrAll = -( np.array([files[i].instrumentCalibrationA4+A4All[i] for i in range(numFiles)]))
    
    # Find binning (All are equal through testing)
    binning = files[0].binning
    
    if binning==1:
        if A4InstrAll.shape[1]==155: #MULTIFLEXX
            A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,5,binning))
        elif A4InstrAll.shape[1]==32: # FLATCONE
            A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,1,binning))
        else:
            A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,8,binning))
    else:
        A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,8,binning))
    
    ####################################################################### Assume that all energies have same A4
    A4InstrAll = A4InstrAll.reshape(numFiles,A4InstrAll[0].shape[0],-1)[:,:,0]

    # Generate measured points in A3-A4 space

    #@_tools.my_timer_N()
    #def genPointsAndBoundary(A3All,A4InstrAll,numFiles,I,Norm,Mon):
    points = []

    for i in range(numFiles):
        X,Y = [x.flatten() for x in np.meshgrid(A3All[i],A4InstrAll[i],indexing='ij')]
        points.append([X,Y])
    
    PosAll = np.concatenate(points,axis=1)
    unique,uindex,count = np.unique(PosAll,axis=1,return_index=True,return_counts=True)
    
    if np.sum(count>1)>0: # If there is any duplicate points
        BoundPoly= [convexHullPoints(points[i][0].flatten(),points[i][1].flatten()) for i in range(numFiles)]

        mask = np.ones(PosAll.shape[1],dtype=bool)
        mask[uindex] = False
        doublePoints = PosAll[:,mask]
        kdtree = KDTree(unique.T)
        
        doubleIndex = kdtree.query(np.round(doublePoints,binningDecimals).T,distance_upper_bound=np.power(10,-binningDecimals*1.0)*1.1)[1]
        #doubleIndex = np.concatenate([np.where(np.all(x==unique.T,axis=1)) for x in doublePoints.T]).reshape(-1)
        
        points = unique
        shape = I.shape[2]

        IReshape = I.reshape(-1,shape)
        NormReshape = Norm.reshape(-1,shape)
        MonReshape = Mon.reshape(-1,shape)

        doubleI = IReshape[mask,:]
        doubleNorm = NormReshape[mask,:]
        doubleMon = MonReshape[mask,:]

        Isorted = IReshape[uindex,:]
        Normsorted = NormReshape[uindex,:]
        Monsorted = MonReshape[uindex,:]

        Isorted[doubleIndex,:]+=doubleI
        Normsorted[doubleIndex,:]=np.nanmean([Normsorted[doubleIndex,:],doubleNorm],axis=0)
        Monsorted[doubleIndex,:]+=doubleMon
        
    else:
        BoundPoly = False
        # Sort measured points first in y and then x direction
        index = np.lexsort((unique[1], unique[0]))
        shape = I.shape[2] #(64 or 8 depending on instrument and binning)
        Isorted =  np.concatenate(I,axis=0)[index,:] #.reshape(-1,shape)
        Normsorted = np.concatenate(Norm,axis=0)[index,:]#Norm.reshape(-1,shape)[index,:]
        Monsorted = np.concatenate(Mon,axis=0)[index,:]#Mon.reshape(-1,shape)[index,:]

    #   return points,BoundPoly,Isorted,Normsorted,Monsorted

    #points,BoundPoly,Isorted,Normsorted,Monsorted = genPointsAndBoundary(A3All,A4InstrAll,numFiles,I,Norm,Mon)

    if numFiles==1:
        points = [np.array(points).reshape(2,-1)]
        numGroups = 1
    else:
        numGroups = False
    polygons,GoodPolyPoints = voronoiTessellation(points=points ,plot = plotTessellation,Boundary = BoundPoly, numGroups=numGroups)


    # Sort centroids (i.e. polygons) like measurement points
    #@_tools.my_timer_N()
    #def calcCentroids(GoodPolyPoints):
    centroids = np.array([centeroidnp(x) for x in GoodPolyPoints]).T
    #    return centroids

    #centroids = calcCentroids(GoodPolyPoints)  


    #@_tools.my_timer_N()
    #def sortPoints(points,centroids):
    if isinstance(points,list):
        X = np.concatenate(points,axis=1).T
    else:
        X = points.T

    Y = centroids.T

    kdtree = KDTree(X)
    _,A = kdtree.query(Y)

    _,SortUindex,SortCount = np.unique(A,return_index=True,return_counts=True)

    if np.sum(SortCount>1)!=0:
        raise AttributeError('The number of points connecting the centroids from Tessellation and points are not equal...')
    centInd = SortUindex

    #@_tools.my_timer_N()
    #def calculateQ(GoodPolyPoints,centInd,files):
    sortedPolyPoints = GoodPolyPoints[centInd]
    factorsqrtEK = 0.694692
    
    # Calcualte k vectors
    Ei = files[0].Ei
    ki = np.sqrt(Ei)*factorsqrtEK
    kf = np.sqrt(Ei-files[0].energy[0,:,:].mean(axis=0))*factorsqrtEK
    
    
    # Convert to Q-space
    ## Qx = ki-kf*cos(A4), Qy = -kf*sin(A4)
    QX = np.array([ki-np.outer(np.cos(np.deg2rad(p[:,1])),kf) for p in sortedPolyPoints])
    QY = np.array([-np.outer(np.sin(np.deg2rad(p[:,1])),kf) for p in sortedPolyPoints])
        
    Theta = np.array([p[:,0].reshape((-1,1))*np.pi/180.0 for p in sortedPolyPoints])
    
    QRX = np.array([QX[i]*np.cos(Theta[i])-QY[i]*np.sin(Theta[i]) for i in range(QX.shape[0])])
    QRY = np.array([QY[i]*np.cos(Theta[i])+QX[i]*np.sin(Theta[i]) for i in range(QX.shape[0])])
    
    # Find common axis limits
    qxmin = np.min([np.min(val) for val in QRX])
    qymin = np.min([np.min(val) for val in QRY])
    qxmax = np.max([np.max(val) for val in QRX])
    qymax = np.max([np.max(val) for val in QRY])
    
    QXlim = np.max(np.abs([qxmin,qxmax]))
    QYlim = np.max(np.abs([qymin,qymax]))
    E = np.mean(files[0].energy,axis=(0,1))
    #    return QRX,QRY,E,QXlim,QYlim

    #QRX,QRY,E,QXlim,QYlim = calculateQ(GoodPolyPoints,centInd,files)
    
    if len(planes)==0:
        planes = range(len(E))
        
    plots = len(planes)
    if not returnPatches:
        if ax is None: # Create needed axes
            if singleFigure: # pragma: no cover
                # create only one
                rows,cols = figureRowColumns(plots)
                fig,ax = plt.subplots(nrows=rows, ncols=cols)
                ax = np.array(ax).flatten()
        if singleFigure: # pragma: no cover
            if ax is None:
                ax = plt.figure().gca()
        else:
            if ax is None: # pragma: no cover
                ax = [plt.figure().gca() for _ in range(plots)]
            
    counter = 0
    if returnPatches:
        ReturnPatches = []
        Energies = []
    for plane in planes:
        
        #@_tools.my_timer_N()
        #def binPlanes(plane,Isorted,Normsorted,Monsorted):
        subplanes = len(np.array([plane]).flatten())
        # Check if plane inpu is single plane
        if subplanes==1:
            plotPlane = plane
            IntensityBin = np.divide(Isorted[:,plane],Normsorted[:,plane]*Monsorted[:,plane])+1e-20
            IntensityBin = np.ma.masked_invalid(IntensityBin)

        else:
            plotPlane = int(np.mean(plane))
            IntensityBin = np.divide(np.nansum(Isorted[:,plane],axis=1),np.nanmean(Normsorted[:,plane],axis=1)*np.nansum(Monsorted[:,plane],axis=1))+1e-20
            IntensityBin = np.ma.masked_invalid(IntensityBin)
        #    return plotPlane,IntensityBin,subplanes
        
        #plotPlane,IntensityBin,subplanes = binPlanes(plane,Isorted,Normsorted,Monsorted)
         # Generate polygons in Qspace
        

        #@_tools.my_timer_N()
        #def genPatchesAndCollection(QRX,QRY,plotPlane):
        #patches = [Polygon(np.array([QRX[i][:,plotPlane],QRY[i][:,plotPlane]]).T) for i in range(len(QRX))]
        patches = [Polygon(np.array([qrx[:,plotPlane],qry[:,plotPlane]]).T) for qrx,qry in zip(QRX,QRY)]
        pcollection = PatchCollection(patches)
        #    return pcollection

        #pcollection = genPatchesAndCollection(QRX,QRY,plotPlane)
        currentInt = IntensityBin
        

        #@_tools.my_timer_N()
        #def plotter(pcollection,currentInt,counter,ax,QXlim,QYlim,E,plotPlane,plane,subplanes):
        if log==True:
            pcollection.set_array(np.log10(currentInt+1e-20))
        else:
            pcollection.set_array(currentInt)
        if returnPatches:
            pcollection.set_edgecolor('None')
            ReturnPatches.append(pcollection)
            Energies.append(np.mean(E[plane]))
            #continue
        else:
            pcollection.set_edgecolor('face')
            currIntMin = np.max([np.nanmin(currentInt),0.0])
            pcollection.set_clim(currIntMin,np.nanmax(currentInt))
            
            ax[counter].add_collection(pcollection)
            ax[counter].set_xlim(-QXlim,QXlim)
            ax[counter].set_ylim(-QYlim,QYlim)
            ax[counter].colorbar = ax[counter].get_figure().colorbar(ax[counter].collections[0], ax=ax[counter],format=ticker.FuncFormatter(fmt))
            
            ax[counter].collections[0].set_clim(currIntMin,np.max(currentInt))
            if subplanes==1:
                ax[counter].set_title('Energy {0:.3f} meV - plane {1}'.format(E[plotPlane],plane))
            else:
                ax[counter].set_title('Energy {0:.3f} meV - planes '.format(np.mean(E[plane]))+\
                    ','.join([str(x) for x in plane]))
            counter +=1

            #return ax,counter
        
        #ax,counter = plotter(pcollection,currentInt,counter,ax,QXlim,QYlim,E,plotPlane,plane,subplanes)

    if returnPatches:
        return ReturnPatches,Energies
    else:
        return ax


#def plotQPatches(dataFiles,ax=None,planes=[],binningDecimals=3,log=False,returnPatches=False,A4Extend=0.2,A3Extend=0.5,singleFigure=False,plotTessellation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
#    """Plot data files together with pixels created around each point in Q space. See :doc:`Voronoi Tessellation<../../InDepthDocumentation/VoronoiTessellation>` for further information.
#
#    .. warning::
#        This method plots all measurement points unless they are literally on top of each other and is thus really slow! Binning 8 planes for two files takes approximately
#        3.5 minutes. Alternatively use binning, i.e. plotQPlane.
#
#
#    Args:
#        
#        - dataFiles (DataFiles): single file or list of files to be binned together
#
#    Kwargs:
#
#        - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)
#
#        - planes (list (of lists)): Planes to be plotted and binned (default [])
#
#        - binningDecimals (int): Number of decimal places Q positions are rounded before binning (default 3)
#        
#        - log (bool): Whether or not to plot intensities as logarithm (default False)
#
#        - returnPatches (bool): If true the method returns the patches otherwise plotted in the given axis (default False)
#        
#        - A4Extend (float): Angle value with which the boundary is extended away from points in A4 direction (default 0.2)
#        
#        - A3Extend (float): Angle value with which the boundary is extended away from points in A3 direction (default 0.5)
#
#        - singleFigure (bool): If true, all planes are plotted in same figure (default False)
#
#        - plotTessellation (bool): Plot Tessellation of points (default False)
#
#        - Ei_err (float): Tolerance of E_i for which the values are equal (default = 0.05)
#
#        - temperature_err (float): Tolerance of temperature for which the values are equal (default = 0.2)
#        
#        - magneticField_err (float): Tolerance of magnetic field for which the values are equal (default = 0.2)
#        
#        - electricField_err (float): Tolerance of electric field for which the values are equal (default = 0.2)
#
#    Returns:
#        
#        - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.
#
#    Raises:
#
#        - AttributeError
#
#    Examples: 
#
#    The following example will combine the two files and plot all of the available planes in different figures.
#
#    >>> DS = DataSet.DataSet(convertedFiles=[--.nxs,---.nxs])
#    >>> plt.figure()
#    >>> ax = plt.gca()
#    >>>
#    >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax)
#
#    If only a subset of planes or different planes are to be combined the following will achieve this:
#
#    >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax,planes=[0,1,2,3,[4,5,6],[8,9]])
#
#    Here planes 0 through 3 are plotted separately while 4,5, and 6 as well as 8 and 9 are binned.
#
#    .. note::
#        Binning planes from different analysers might result in nonsensible binnings.
#
#    """
#    #if dimension!='2D':
#    #    raise NotImplementedError('Only 2D plotting is currently supported')
#    
#    if not isinstance(ax, (list,)) and ax is not None:
#        ax = np.array([ax])
#    
#    if not isinstance(planes, (list,)):
#        planes = np.array([planes])
#        
#    if not ax is None:
#        if singleFigure and np.array([ax]).size != 1:
#            raise AttributeError('Single figure chosen but multiple axes given ({}).'.format(np.array([ax]).size))
#        
#        elif not singleFigure and len(ax) != len(planes) and not len(planes)==0:
#            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,len(planes)))
#            
#    
#    dataFiles = np.asarray(dataFiles)
#    numFiles = len(dataFiles)
#
#    
#    if numFiles>1:
#        comparison = np.array([np.all([np.isclose(dataFiles[0].Ei,dataFiles[i+1].Ei,atol=Ei_err) for i in range(numFiles-1)]),\
#                  np.all([compareNones(dataFiles[0].temperature,dataFiles[i+1].temperature,temperature_err) for i in range(numFiles-1)]),\
#                  np.all([compareNones(dataFiles[0].magneticField,dataFiles[i+1].magneticField,magneticField_err) for i in range(numFiles-1)]),\
#                  np.all([compareNones(dataFiles[0].electricField,dataFiles[i+1].electricField,electricField_err) for i in range(numFiles-1)]),\
#                  np.all([dataFiles[0].binning==dataFiles[i+1].binning for i in range(numFiles-1)])])
#        
#        tests = np.array(['Ei','Temperature','Magnetic Field','Electric Field','Binning'])
#        
#        if not np.all(comparison):
#            errors = np.array(1-comparison,dtype=bool)
#            raise AttributeError('Attributes for the datafiles are not the same! Difference is in :\n'+','.join([x for x in tests[errors]])+'\nIf the files are to be binned anyway change the tolerence limits.')
#    
#    Ishape = dataFiles[0].I.shape
#    if not ax is None:
#        if not singleFigure and len(ax) != Ishape[4] and len(planes) == 0: # Plot all planes in provided axes
#            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,Ishape[4]))
#
#    
#    IAll = np.array([dataFiles[i].I[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)]) # into shape sum(A3),104,64 for CAMEA
#    NormAll = np.array([dataFiles[i].Norm[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
#    MonitorAll = np.array([dataFiles[i].Monitor[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
#  
#    I = np.concatenate(IAll,axis=0)
#    Norm = np.concatenate(NormAll,axis=0)
#    Mon = np.concatenate(MonitorAll,axis=0)
#    
#    QxAll = np.array([dataFiles[i].qx[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
#    QyAll = np.array([dataFiles[i].qy[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
#    Qx = np.concatenate(QxAll,axis=0)
#    Qy = np.concatenate(QyAll,axis=0)
#
#    
#    if len(planes)==0:
#        planes = range(len(I.shape[-1]))
#    
#    plots = len(planes)
#    if not returnPatches: # only check axes if the user wants to plot in these
#        if ax is None: # Create needed axes
#            if singleFigure: # create only one
#                rows,cols = figureRowColumns(plots)
#                fig,ax = plt.subplots(nrows=rows, ncols=cols)
#                ax = np.array(ax).flatten()
#        if singleFigure:
#            if ax is None:
#                ax = plt.figure().gca()
#        else:
#            if ax is None:
#                ax = [plt.figure().gca() for _ in range(plots)]
#    counter = 0
#
#    if returnPatches:
#        ReturnPatches = []
#        Energies = []
#    for plane in planes:
#        mp = []
#        for i in range(len(dataFiles)):
#            xx = boundaryQ(dataFiles[i],plane,A4Extend=A4Extend,A3Extend=A3Extend)
#            polygons = [PolygonS(x.T) for x in xx.transpose(1,0,2)]
#            if isinstance(plane,list):
#                if len(plane)>1:
#                    mplocal = polygons[0]
#                    for j in range(len(polygons)-1):
#                        mplocal = mplocal.union(polygons[j+1])
#                    mp.append(mplocal)
#                else:
#                    mp.append(polygons[0])
#            else:
#                mp.append(polygons[0])
#        
#        
#        if len(mp)>1:
#            boundary = mp[0]
#            for i in range(len(mp)-1):
#                boundary = boundary.union(mp[i+1])
#            boundary = [boundary]
#        else:
#            boundary = mp
#
#        if isinstance(plane,list) or isinstance(plane,np.ndarray):         
#            IAlive = []
#            NormAlive = []
#            MonAlive = []
#            QxAlive = []
#            QyAlive = []
#            for i in range(len(plane)):
#                alive = np.logical_not(np.isnan(Norm[:,:,plane[i]]))
#                IAlive.append(I[alive,plane[i]])
#                NormAlive.append(Norm[alive,plane[i]])
#                MonAlive.append(Mon[alive,plane[i]])
#                
#                QxAlive.append(Qx[alive,plane[i]])
#                QyAlive.append(Qy[alive,plane[i]])
#            IAlive = np.concatenate(IAlive)
#            NormAlive = np.concatenate(NormAlive)
#            MonAlive = np.concatenate(MonAlive)
#            QxAlive = np.concatenate(QxAlive)
#            QyAlive = np.concatenate(QyAlive)
#        else:
#            alive = np.logical_not(np.isnan(Norm[:,:,plane]))
#            IAlive = I[alive,plane]
#            NormAlive = Norm[alive,plane]
#            MonAlive = Mon[alive,plane]
#            QxAlive = Qx[alive,plane]
#            QyAlive = Qy[alive,plane]
#            
#
#        points = np.array([QxAlive,QyAlive])
#        unique,uindex = np.unique(np.round(points,binningDecimals),axis=1,return_index=True)
#        if unique.shape[1]!=points.shape[1]:
#            #print('BINNING!')
#            mask = np.ones(points.shape[1],dtype=bool)
#            mask[uindex] = False
#            doublePoints = points[:,mask]
#            kdtree = KDTree(unique.T)
#            doubleIndex = kdtree.query(np.round(doublePoints,binningDecimals).T,distance_upper_bound=np.power(10,-binningDecimals*1.0)*1.1)[1]
#
#            points = unique
#           
#            Isorted = IAlive[uindex]
#            Normsorted = NormAlive[uindex]
#            Monsorted = MonAlive[uindex]
#
#            IAliveDouble = IAlive[mask]
#            NormAliveDouble = NormAlive[mask]
#            MonAliveDouble = MonAlive[mask]
#
#            Isorted[doubleIndex]+=IAliveDouble
#            Normsorted[doubleIndex]=np.mean([Normsorted[doubleIndex],NormAliveDouble],axis=0)
#            Monsorted[doubleIndex]+=MonAliveDouble
#            currentInt = np.divide(Isorted,Normsorted*Monsorted)
#
#        else:
#            #print('BINNING! 2')
#            #pointIndex = np.lexsort((unique[1], unique[0]))
#            currentInt = np.divide(IAlive,NormAlive*MonAlive)[uindex]
#            
#         
#        polygons,GoodPolyPoints = voronoiTessellation([unique],plot = plotTessellation,Boundary = boundary)
#        
#        centroids = np.array([np.array(x.centroid.coords).reshape(2) for x in polygons]).T
#        
#
#        X = unique.T
#        Y = Y = centroids.T
#        #print(X.shape)
#        #print(Y.shape)
#
#        #plt.figure()
#        #plt.scatter(X[:,0],X[:,1],c='r')
#        #plt.scatter(Y[:,0],Y[:,1],c='b')
#
#
#        kdtree = KDTree(X)
#
#        A = kdtree.query(Y)[1]#,distance_upper_bound=0.02)[1]
#        #print(A.shape)
#        #plt.scatter(X[A,0],X[A,1],c=np.linspace(0,1,len(A)),s=5)
#    
#
#        _,SortUindex,SortCount = np.unique(A,return_index=True,return_counts=True)
#        #print(_.shape)
#        if np.sum(SortCount>1)!=0:
#            #plt.scatter(X[_,0][SortCount>1],X[_,1][SortCount>1],c='k')
#            raise AttributeError('The number of points tieing the centroids and Q poinst together are not equal, difference is {}. Try extending A3 and A4.'.format(np.sum(SortCount>1)))
#        patchIndex = SortUindex
#        E = dataFiles[0].energy
#        #patches = [Polygon(np.array([np.array(x.boundary.coords)[:,0],np.array(x.boundary.coords)[:,1]]).T) for x in polygons[patchIndex]]
#        pcollection = PolyCollection([np.array([np.array(x.boundary.coords)[:,0],np.array(x.boundary.coords)[:,1]]).T for x in polygons[patchIndex]])
#        #pcollection = PatchCollection(patches)
#        
#        try:
#            bpoints = np.array(boundary[0].boundary.coords)
#        except:
#            bpoints = np.concatenate([np.array(x.boundary.coords) for x in boundary[0]])
#            
#        qxmin = np.min(bpoints[:,0])
#        qymin = np.min(bpoints[:,1])
#        qxmax = np.max(bpoints[:,0])
#        qymax = np.max(bpoints[:,1])
#        
#        QXlim = np.max(np.abs([qxmin,qxmax]))
#        QYlim = np.max(np.abs([qymin,qymax]))
#        
#        if log==True:
#            pcollection.set_array(np.log10(currentInt+1e-20))
#        else:
#            pcollection.set_array(currentInt)
#        
#        currIntMin = np.max([np.nanmin(currentInt),0.0])
#        pcollection.set_clim(currIntMin,np.nanmax(currentInt))
#        #return pcollection
#        if returnPatches:
#            pcollection.set_edgecolor('None')
#            ReturnPatches.append(pcollection)
#            counter +=1
#            Energies.append(np.mean(E[:,:,:,:,plane]))
#            continue
#        else:
#            pcollection.set_edgecolor('face')
#        ax[counter].add_collection(pcollection)
#        ax[counter].set_xlim(-QXlim,QXlim)
#        ax[counter].set_ylim(-QYlim,QYlim)
#        ax[counter].axes.grid(True)
#        ax[counter].get_figure().colorbar(ax[counter].collections[0], ax=ax[counter],format=ticker.FuncFormatter(fmt))
#        
#        ax[counter].collections[0].set_clim(currIntMin,np.max(currentInt))
#        if not isinstance(plane,list):
#            ax[counter].set_title('Energy {0:.3f} meV - plane {1}'.format(np.mean(E[:,:,:,:,plane]),plane))
#        else:
#            if len(plane) == 1:
#                ax[counter].set_title('Energy {0:.3f} meV - plane {1}'.format(np.mean(E[:,:,:,:,plane]),plane))
#            else:
#                ax[counter].set_title('Energy {0:.3f} meV - planes '.format(np.mean(E[:,:,:,:,plane]))+\
#                  ','.join([str(x) for x in plane]))
#        counter +=1
#    
#    if returnPatches:
#        return ReturnPatches,Energies
#    else:
#        return ax
#
#
#@_tools.my_timer_N()

# @_tools.KwargChecker()
# def voronoiTessellation(points,plot=False,Boundary=False,numGroups=False):
#     """Generate individual pixels around the given datapoints.

#     Args:

#         - points (list of list of points): Data points to generate pixels in shape [files,XY,N] i.e. [1,2,N] for one file with N points

#     Kwargs:

#         - plot (bool): If True, method plots pixels created with green as edge bins and red as internal (default False)

#         - Boundary (list of Polygons): List of Shapely polygons constituting the boundaries (Default False)


#     """

#     if numGroups == False:
#         numGroups = len(points)

#     if Boundary==False:
#         BoundPoly= [convexHullPoints(np.array(points[i][0]).flatten(),np.array(points[i][1]).flatten()) for i in range(numGroups)]
#     else:
#         BoundPoly = Boundary#[PolygonS(x.T) for x in Boundary]

#     if numGroups == 1:

#         combiPoly = BoundPoly[0]
#         pointsX = np.array([points[0][0].flatten()])[0]
#         pointsY = np.array([points[0][1].flatten()])[0]
#     else: # Combine all files
#         combiPoly = BoundPoly[0].union(BoundPoly[1])
#         for i in range(len(BoundPoly)-2):
#             combiPoly = combiPoly.union(BoundPoly[i+2])
#         if Boundary==False:
#             pointsX = np.concatenate([points[i][0].flatten() for i in range(numGroups)])
#             pointsY = np.concatenate([points[i][1].flatten() for i in range(numGroups)])
#         else:
#             pointsX = points[0]
#             pointsY = points[1]
        
#     containsAllPoints=np.all([combiPoly.contains(PointS(pointsX[i],pointsY[i])) for i in range(len(pointsX))])
#     if not containsAllPoints:
#         plt.figure()
#         plt.scatter(pointsX,pointsY,c='b')
#         boundaryXY = np.array(combiPoly.boundary.coords)
#         plt.plot(boundaryXY[:,0],boundaryXY[:,1],c='r')
#         raise AttributeError('The provided boundary does not contain all points')
#     # Add extra points to ensure that area is finite
#     extraPoints = np.array([[np.mean(pointsX),np.max(pointsY)+50],[np.mean(pointsX),np.min(pointsY)-50],\
#                              [np.min(pointsX)-50,np.mean(pointsY)],[np.max(pointsX)+50,np.mean(pointsY)],\
#                              [np.min(pointsX)-50,np.max(pointsY)+50],[np.min(pointsX)-50,np.min(pointsY)-50],\
#                              [np.max(pointsX)+50,np.max(pointsY)+50],[np.max(pointsX)+50,np.min(pointsY)-50]])

#     AllPoints = np.array([np.concatenate([pointsX,extraPoints[:,0]]),np.concatenate([pointsY,extraPoints[:,1]])])
    

#     vor = Voronoi(AllPoints.T)
#     regions = np.array([reg for reg in vor.regions])
#     boolval = np.array([len(x)>2 and not -1 in x for x in regions]) # Check if region has at least 3 points and is not connected to infinity (-1))
        
#     PolyPoints = np.array([vor.vertices[reg,:] for reg in regions[boolval]])
#     polygons = np.array([PolygonS(X) for X in PolyPoints])

#     insidePolygonsBool = np.array([combiPoly.contains(P) for P in polygons])

#     edgePolygonsBool = np.logical_not(insidePolygonsBool)
    
#     intersectionPolygon = []
#     for poly in polygons[edgePolygonsBool]:
#         inter = poly.intersection(combiPoly)
#         if not isinstance(inter,PolygonS): # Not a simple polygon
#             inter = inter[np.argmax([x.area for x in inter])] # Return the polygon with biggest area inside boundary
#         intersectionPolygon.append(inter)
    
#     Polygons = np.concatenate([polygons[np.logical_not(edgePolygonsBool)],intersectionPolygon])
    
    
#     if plot or len(pointsX)!=len(Polygons): # pragma: no cover
#         plt.figure()
#         insiders = np.logical_not(edgePolygonsBool)
        
#         [plt.plot(np.array(inter.boundary.coords)[:,0],np.array(inter.boundary.coords)[:,1],c='r') for inter in polygons[insiders]]
#         [plt.plot(np.array(inter.boundary.coords)[:,0],np.array(inter.boundary.coords)[:,1],c='g') for inter in intersectionPolygon]
#         [plt.plot(np.array(bound.boundary.coords)[:,0],np.array(bound.boundary.coords)[:,1],'-.',c='r') for bound in BoundPoly]
#         plt.scatter(extraPoints[:,0],extraPoints[:,1])

#         from scipy.spatial import voronoi_plot_2d
#         voronoi_plot_2d(vor)
#     if not len(pointsX)==len(Polygons):
#         raise AttributeError('The number of points given({}) is not the same as the number of polygons created({}). This can be due to many reasons, mainly:\n - Points overlap exactly\n - Points coinsides with the calculated edge\n - ??'.format(len(pointsX),len(Polygons)))


#     return Polygons,np.array([np.array(P.boundary.coords[:-1]) for P in Polygons])
#
#
#
#
#@_tools.my_timer_N(N=0)
#def voronoiTessellationOPTIMIZED2(points,plot=False,Boundary=False): # pragma: no cover
#    """Generate individual pixels around the given datapoints.
#
#    Args:
#
#        - points (list of list of points): Data points to generate pixels in shape [files,XY,N] i.e. [1,2,N] for one file with N points
#
#    Kwargs:
#
#        - plot (bool): If True, method plots pixels created with green as edge bins and red as internal (default False)
#
#        - Boundary (lost of Polygons): List of Shapely polygons constituting the boundaries (Default False)
#
#
#    """
#    numGroups = len(points)
#    
#    if Boundary==False:
#        BoundPoly= [convexHullPoints(points[i][0].flatten(),points[i][1].flatten()) for i in range(numGroups)]
#    else:
#        BoundPoly = Boundary 
#        
#    if numGroups == 1:
#        combiPoly = BoundPoly[0]
#        pointsX = np.array([points[0][0].flatten()])[0]
#        pointsY = np.array([points[0][1].flatten()])[0]
#    else: # Combine all files
#        combiPoly = BoundPoly[0].union(BoundPoly[1])
#        for i in range(len(BoundPoly)-2):
#            combiPoly = combiPoly.union(BoundPoly[i+2])
#        if Boundary==False:
#            pointsX = np.concatenate([points[i][0].flatten() for i in range(numGroups)])
#            pointsY = np.concatenate([points[i][1].flatten() for i in range(numGroups)])
#        else:
#            pointsX = points[0]
#            pointsY = points[1]
#        
#    containsAllPoints=np.all(contains(combiPoly,pointsX,pointsY))
#    
#    
#    if not containsAllPoints:
#        raise AttributeError('The provided boundary does not contain all points')
#    # Add extra points to ensure that area is finite
#    extraPoints = np.array([[np.mean(pointsX),np.max(pointsY)+50],[np.mean(pointsX),np.min(pointsY)-50],\
#                             [np.min(pointsX)-50,np.mean(pointsY)],[np.max(pointsX)+50,np.mean(pointsY)],\
#                             [np.min(pointsX)-50,np.max(pointsY)+50],[np.min(pointsX)-50,np.min(pointsY)-50],\
#                             [np.max(pointsX)+50,np.max(pointsY)+50],[np.max(pointsX)+50,np.min(pointsY)-50]])
#    
#    
#    AllPoints = np.array([np.concatenate([pointsX,extraPoints[:,0]]),np.concatenate([pointsY,extraPoints[:,1]])])
#    
#    
#    vor = Voronoi(AllPoints.T)
#    
#    regions = np.array(vor.regions)
#       
#    boolval = np.array([len(x)>2 and not -1 in x for x in regions]) # Check if region has at least 3 points and is not connected to infinity (-1))
#    
#    PolyPoints = np.array([vor.vertices[reg,:] for reg in regions[boolval]])
#    
#    
#    def genPolygon(PolyPoints):
#        return PolygonS(PolyPoints)
#
#    genPolygon_vectorized = np.vectorize(genPolygon,otypes=[PolygonS])
#
#    polygons = genPolygon_vectorized(PolyPoints)#
#
#    insidePolygonsBool = np.array([combiPoly.contains(P) for P in polygons])
#    edgePolygonsBool = np.logical_not(insidePolygonsBool)
#    
#    def intersectionLoop(polygons,edgePolygonsBool,combiPoly):
#        intersectionPolygon = []
#        for poly in polygons[edgePolygonsBool]:
#            inter = poly.intersection(combiPoly)
#            if not isinstance(inter,PolygonS): # Not a simple polygon
#                inter = inter[np.argmax([x.area for x in inter])] # Return the polygon with biggest area inside boundary
#            intersectionPolygon.append(inter)
#        return intersectionPolygon
#    intersectionPolygon = intersectionLoop(polygons,edgePolygonsBool,combiPoly)
#
#    Polygons = np.concatenate([polygons[np.logical_not(edgePolygonsBool)],intersectionPolygon])
#
#    if not len(pointsX)==len(Polygons):
#        raise AttributeError('The number of points given({}) is not the same as the number of polygons created({}). This can be due to many reasons, mainly:\n - Points overlap exactly\n - Points coinsides with the calulated edge\n - ??'.format(len(pointsX),len(Polygons)))
#
#    return Polygons,np.concatenate([PolyPoints[insidePolygonsBool],np.array([np.array(P.boundary.coords[:-1]) for P in intersectionPolygon])],axis=0)














@_tools.KwargChecker() # Following function is not used
def boundaryQ(file,plane,A4Extend=0.0,A3Extend=0.0): # pragma: no cover
    """Calculate the boundary of a given scan in Q space
    A4Extend: in degrees
    A3Extend: in degrees
    """
    energy = file.energy[:,0,0,:,:]
    A3 = file.A3+file.A3Off
    
    A4 = file.A4-file.A4Off
    Ei = file.Ei
        
    InstrumentA4 = file.instrumentCalibrationA4.reshape(energy.shape[1],-1)[:,plane]
    
    factorsqrtEK = 0.694692
    InstA4 = (InstrumentA4-A4)*np.pi/180.0 
    A4Min = np.min(InstA4,axis=0)
    A4Max = np.max(InstA4,axis=0)
    
    InstrumentEnergy = IC[:,4].reshape(energy.shape[1],-1)[:,plane] # TODO: IC is not defined before usage! Should be loaded from instrument?
    
    kf = np.sqrt(InstrumentEnergy)*factorsqrtEK
       
    kfmin = np.min(kf,axis=0)
    kfmax = np.max(kf,axis=0)
    
    if not isinstance(kfmin,list): # If only one plane, reshape kfmin/max
        kfmin=np.array([kfmin])
        kfmax=np.array([kfmax])
        A4Min = np.array([A4Min])
        A4Max = np.array([A4Max])
    kfmin.shape= (-1)
    kfmax.shape= (-1)
    A4Min.shape= (-1)
    A4Max.shape= (-1)
    
    A4Min-=A4Extend*np.pi/180.0
    A4Max+=A4Extend*np.pi/180.0
    A3 = np.linspace(np.min(A3)-A3Extend,np.max(A3)+A3Extend,len(A3))
    
    
    ki = np.sqrt(Ei)*factorsqrtEK
    
    #### Qx = ki-kf*cos(A4), Qy = -kf*sin(A4)
    
    # inner line
    qxInner = ki-kfmin*np.cos(A4Min)
    qyInner = -kfmin*np.sin(A4Min)
    qxOuter = ki-kfmax*np.cos(A4Max)
    qyOuter = -kfmax*np.sin(A4Max)
    QInner = np.zeros((qxInner.shape[0],2,len(A3)))
    QOuter = np.zeros_like(QInner)
    for i in range(len(qxInner)):
        QInner[i,0] = qxInner[i]*np.cos(A3*np.pi/180.0)-qyInner[i]*np.sin(A3*np.pi/180.0)
        QInner[i,1] = qyInner[i]*np.cos(A3*np.pi/180.0)+qxInner[i]*np.sin(A3*np.pi/180.0)
        QOuter[i,0] = qxOuter[i]*np.cos(A3*np.pi/180.0)-qyOuter[i]*np.sin(A3*np.pi/180.0)
        QOuter[i,1] = qyOuter[i]*np.cos(A3*np.pi/180.0)+qxOuter[i]*np.sin(A3*np.pi/180.0)
        
    
    A4Values = np.array([np.linspace(A4Min[i],A4Max[i],50) for i in range(len(A4Min))])
    kfValues = np.array([np.linspace(kfmin[i],kfmax[i],50) for i in range(len(kfmin))])
    
    
    QStart = np.array([ki-kfValues*np.cos(A4Values),-kfValues*np.sin(A4Values)])
    
    QxStartA3 = QStart[0]*np.cos(A3[0]*np.pi/180.0)-QStart[1]*np.sin(A3[0]*np.pi/180.0)
    QxStopA3 = QStart[0]*np.cos(A3[-1]*np.pi/180.0)-QStart[1]*np.sin(A3[-1]*np.pi/180.0)
    QyStartA3 = QStart[1]*np.cos(A3[0]*np.pi/180.0)+QStart[0]*np.sin(A3[0]*np.pi/180.0)
    QyStopA3 = QStart[1]*np.cos(A3[-1]*np.pi/180.0)+QStart[0]*np.sin(A3[-1]*np.pi/180.0)    
    return np.array([np.concatenate([np.flip(QInner[:,0,:],axis=1)[:,:-1],QxStartA3[:,:-1],QOuter[:,0,:-1],np.flip(QxStopA3,axis=1)[:,:-1]],axis=-1),np.concatenate([np.flip(QInner[:,1,:],axis=1)[:,:-1],QyStartA3[:,:-1],QOuter[:,1,:-1],np.flip(QyStopA3,axis=1)[:,:-1]],-1)])



def convexHullPoints(A3,A4):
    """Calculate the convex hull of rectangularly spaced A3 and A4 values"""
    A3Unique = np.unique(A3)
    A4Unique = np.unique(A4)
    
    A3Step = np.diff(A3Unique)[[0,-1]]*0.5
    A4Step = np.diff(A4Unique)[[0,-1]]*0.5

    addLeft = np.array(np.meshgrid(A3Unique[0]-A3Step[0],A4Unique)).reshape((2,-1))
    addRight= np.array(np.meshgrid(A3Unique[-1]+A3Step[1],A4Unique)).reshape((2,-1))
    addBottom=np.array(np.meshgrid(A3Unique,A4Unique[0]-A4Step[0])).reshape((2,-1))
    addTop  = np.array(np.meshgrid(A3Unique,A4Unique[-1]+A4Step[1])).reshape((2,-1))
    corners = np.array([[addLeft[0,0],addBottom[1,0]],[addLeft[0,0],addTop[1,-1]],[addRight[0,-1],addBottom[1,0]],[addRight[0,-1],addTop[1,-1]]]).T
    boundary = np.concatenate([addLeft,addRight,addBottom,addTop,corners],axis=1)
    hullPoints = ConvexHull(np.array([boundary[0],boundary[1]]).T)
    bound = hullPoints.points[hullPoints.vertices].T
    return PolygonS(bound.T)

def isListOfStrings(object):
    if isinstance(object, list):
        isListOfStr = True
        for item in object:
            if not isinstance(item, str):
                isListOfStr=False
                break
        if isListOfStr:
            return object
        else:
            raise AttributeError('Data files provided are not a list of strings or string!')
    elif isinstance(object,str):
        return [object]
    else:
        raise AttributeError('Data files provided are not a list of strings or string!')
    
def isListOfDataFiles(inputFiles):
    returnList = []
    if isinstance(inputFiles,(list,np.ndarray)):
        for file in inputFiles:
            if isinstance(file,MJOLNIR.Data.DataFile.DataFile):
                returnList.append(file)
            elif isinstance(file,str):
                # Check if file exists
                if not os.path.isfile(file):
                    raise AttributeError('Following file does not exist:\n{}'.format(file))
                returnList.append(MJOLNIR.Data.DataFile.DataFile(file))
    elif isinstance(inputFiles,MJOLNIR.Data.DataFile.DataFile):
        returnList.append(inputFiles)
    elif isinstance(inputFiles,str):
        returnList.append(MJOLNIR.Data.DataFile.DataFile(inputFiles))
    else:
        raise AttributeError('File provided is not of type string, list, or DataFile')
    if len(returnList)>1:
        sameSample = [returnList[0].sample==file.sample for file in returnList]
        if not np.all(sameSample):
            raise AttributeError('Files does not have the same sample! Compared to first entry: {}'.format(sameSample))
    return returnList




def calculateGrid3D(X,Y,Z):
    """Generate 3D grid with centers given by X,Y, and Z.
     Args:
        
        X (3D array): 3D array of x values generated by np.meshgrid.
                
        Y (3D array): 3D array of y values generated by np.meshgrid.
                
        Z (3D array): 3D array of z values generated by np.meshgrid.
        
    Example:

    >>> x = np.linspace(-1.5,1.5,20)
    >>> y = np.linspace(0,1.5,10)
    >>> z = np.linspace(-1.0,5.5,66)
    >>> X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
    >>> XX,YY,ZZ = calculateGrid3D(X,Y,Z)

    Now XX is a 21x11x67 array containing all x coordinates of the edges exactly midway between the points. Same goes for YY and ZZ with y and z coordinates respectively.
    """

    xshape = np.array(X.shape)
    if np.any(xshape <= 1):
        raise AttributeError('Provided array has dimension(s) {} of size <= 1'.format(xshape))
    XT = np.zeros((xshape[0]+1,xshape[1]+1,xshape[2]+1))
    YT = np.zeros_like(XT)
    ZT = np.zeros_like(XT)
    
    
    
    dx0 = np.diff(X,axis=0)
    dx1 = np.diff(X,axis=1)
    dx2 = np.diff(X,axis=2)
    dy0 = np.diff(Y,axis=0)
    dy1 = np.diff(Y,axis=1)
    dy2 = np.diff(Y,axis=2)
    dz0 = np.diff(Z,axis=0)
    dz1 = np.diff(Z,axis=1)
    dz2 = np.diff(Z,axis=2)
    
    
    XX = X.copy()
    XX[:-1]-=0.5*dx0
    XX[-1]-=0.5*dx0[-1]
    XX[:,:-1]-=0.5*dx1
    XX[:,-1]-=0.5*dx1[:,-1]
    XX[:,:,:-1]-=0.5*dx2
    XX[:,:,-1]-=0.5*dx2[:,:,-1]
    
    YY = Y.copy()
    YY[:-1]-=0.5*dy0
    YY[-1]-=0.5*dy0[-1]
    YY[:,:-1]-=0.5*dy1
    YY[:,-1]-=0.5*dy1[:,-1]
    YY[:,:,:-1]-=0.5*dy2
    YY[:,:,-1]-=0.5*dy2[:,:,-1]
    
    ZZ = Z.copy()
    ZZ[:-1]-=0.5*dz0
    ZZ[-1]-=0.5*dz0[-1]
    ZZ[:,:-1]-=0.5*dz1
    ZZ[:,-1]-=0.5*dz1[:,-1]
    ZZ[:,:,:-1]-=0.5*dz2
    ZZ[:,:,-1]-=0.5*dz2[:,:,-1]
    
    XT[:-1,:-1,:-1]=XX.copy()
    YT[:-1,:-1,:-1]=YY.copy()
    ZT[:-1,:-1,:-1]=ZZ.copy()
    
    
    XT[-1,:-1,:-1]=XT[-2,:-1,:-1]+dx0[-1]
    XT[:-1,-1,:-1]=XT[:-1,-2,:-1]+dx1[:,-1,:]
    XT[:-1,:-1,-1]=XT[:-1,:-1,-2]+dx2[:,:,-1]
    XT[:-1,-1,-1]=0.5*(XT[:-1,-1,-2]+dx2[:,-1,-1]+XT[:-1,-2,-1]+dx1[:,-1,-1])
    XT[-1,:-1,-1]=0.5*(XT[-1,:-1,-2]+dx2[-1,:,-1]+XT[-2,:-1,-1]+dx0[-1,:,-1])
    XT[-1,-1,:-1]=0.5*(XT[-1,-2,:-1]+dx1[-1,-1,:]+XT[-2,-1,:-1]+dx0[-1,-1,:])
    XT[-1,-1,-1]=(XT[-1,-2,-1]+dx1[-1,-1,-1]+XT[-2,-1,-1]+dx0[-1,-1,-1]+XT[-1,-1,-2]+dx2[-1,-1,-1])/3
    
    YT[-1,:-1,:-1]=YT[-2,:-1,:-1]+dy0[-1]
    YT[:-1,-1,:-1]=YT[:-1,-2,:-1]+dy1[:,-1,:]
    YT[:-1,:-1,-1]=YT[:-1,:-1,-2]+dy2[:,:,-1]
    YT[:-1,-1,-1]=0.5*(YT[:-1,-1,-2]+dy2[:,-1,-1]+YT[:-1,-2,-1]+dy1[:,-1,-1])
    YT[-1,:-1,-1]=0.5*(YT[-1,:-1,-2]+dy2[-1,:,-1]+YT[-2,:-1,-1]+dy0[-1,:,-1])
    YT[-1,-1,:-1]=0.5*(YT[-1,-2,:-1]+dy1[-1,-1,:]+YT[-2,-1,:-1]+dy0[-1,-1,:])
    YT[-1,-1,-1]=(YT[-1,-2,-1]+dy1[-1,-1,-1]+YT[-2,-1,-1]+dy0[-1,-1,-1]+YT[-1,-1,-2]+dy2[-1,-1,-1])/3
    
    ZT[-1,:-1,:-1]=ZT[-2,:-1,:-1]+dz0[-1]
    ZT[:-1,-1,:-1]=ZT[:-1,-2,:-1]+dz1[:,-1,:]
    ZT[:-1,:-1,-1]=ZT[:-1,:-1,-2]+dz2[:,:,-1]
    ZT[:-1,-1,-1]=0.5*(ZT[:-1,-1,-2]+dz2[:,-1,-1]+ZT[:-1,-2,-1]+dz1[:,-1,-1])
    ZT[-1,:-1,-1]=0.5*(ZT[-1,:-1,-2]+dz2[-1,:,-1]+ZT[-2,:-1,-1]+dz0[-1,:,-1])
    ZT[-1,-1,:-1]=0.5*(ZT[-1,-2,:-1]+dz1[-1,-1,:]+ZT[-2,-1,:-1]+dz0[-1,-1,:])
    ZT[-1,-1,-1]=(ZT[-1,-2,-1]+dz1[-1,-1,-1]+ZT[-2,-1,-1]+dz0[-1,-1,-1]+ZT[-1,-1,-2]+dz2[-1,-1,-1])/3
    
    
    return XT,YT,ZT



@_tools.KwargChecker()
def binData3D(dx,dy,dz,pos,data,norm=None,mon=None,bins=None):
    """ 3D binning of data.

    Args:

        - dx (float): Step size in x (required).

        - dy (float): Step size in x (required).

        - dz (float): Step size in x (required).

        - pos (2D array): Position of data points as flattened lists (X,Y,Z) (required).

        - data (array): Flattened data array (required).

    Kwargs:

        - norm (array): Flattened normalization array.

        - mon (array): Flattened monitor array.

        - bins (list of arrays): Bins locating edges in the x, y, and z directions.

    returns:

        Re-binned intensity (and if provided Normalization, Monitor, and Normalization Count) and X, Y, and Z bins in 3 3D arrays.


    Example:

    >>> pos = [Qx,Qy,E]
    >>> Data,bins = DataSet.binData3D(0.05,0.05,0.2,pos,I,norm=Norm,mon=Monitor)

    """

    if bins is None:
        bins = calculateBins(dx=dx,dy=dy,dz=dz,pos=pos)
    if len(pos[0].shape)>1: # Flatten positions
        pos = np.array([x.flatten() for x in pos])
    #NonNaNs = 1-np.isnan(data.flatten())

    #pos = [np.array(x[NonNaNs]) for x in pos]
    HistBins = [bins[0][:,0,0],bins[1][0,:,0],bins[2][0,0,:]]
    intensity =    np.histogramdd(np.array(pos).T,bins=HistBins,weights=data.flatten())[0].astype(data.dtype)

    returndata = [intensity]
    if mon is not None:
        MonitorCount=  np.histogramdd(np.array(pos).T,bins=HistBins,weights=mon.flatten())[0].astype(mon.dtype)
        returndata.append(MonitorCount)
    if norm is not None:
        Normalization= np.histogramdd(np.array(pos).T,bins=HistBins,weights=norm.flatten())[0].astype(norm.dtype)
        
        returndata.append(Normalization)
        
    NormCount =    np.histogramdd(np.array(pos).T,bins=HistBins,weights=np.ones_like(data).flatten())[0].astype(int)
    returndata.append(NormCount)
    return returndata,bins

def calculateBins(dx,dy,dz,pos):
    diffx = np.abs(np.max(pos[0])-np.min(pos[0]))
    diffy = np.abs(np.max(pos[1])-np.min(pos[1]))
    diffz = np.abs(np.max(pos[2])-np.min(pos[2]))
    
    xbins = np.round(diffx/dx).astype(int)+1
    ybins = np.round(diffy/dy).astype(int)+1
    zbins = np.round(diffz/dz).astype(int)+1
    
    _X = np.linspace(np.min(pos[0]),np.max(pos[0]),xbins)
    _Y = np.linspace(np.min(pos[1]),np.max(pos[1]),ybins)
    _Z = np.linspace(np.min(pos[2]),np.max(pos[2]),zbins)
    
    X,Y,Z = np.meshgrid(_X,_Y,_Z,indexing='ij')
    
    XX,YY,ZZ = calculateGrid3D(X=X,Y=Y,Z=Z)
    
    bins=[XX,YY,ZZ]
    return bins

def getNX_class(x,y,attribute):
    try:
        variableType = y.attrs['NX_class']
    except:
        variableType = ''
    if variableType==attribute:
        return x

def getInstrument(file):
    location = file.visititems(lambda x,y: getNX_class(x,y,attribute=b'NXinstrument'))
    return file.get(location)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def figureRowColumns(subplots):
    if subplots<1:
        raise AttributeError('Negative or zero number of subplots requested.')
    if subplots==1:
        rows = 1
        cols = 1
    else:
        subplots = float(subplots)
        startGuess = int(np.ceil(np.sqrt(subplots)))
        for i in np.arange(startGuess,subplots+1,dtype=int):
            if int(i*np.ceil(subplots/i))>=subplots:# and int(i*np.ceil(plots/i))-plots<startGuess:#np.mod(plots,i)==0:
                rows = int(np.ceil(subplots/i))
                cols = i
                break
    return rows,cols


def centeroidnp(arr): # Calculated centroid
    length = arr.shape[0]
    Totsum = np.sum(arr,axis=0)
    return Totsum/length

def compareNones(first,second,margin): # Function to compare
    if hasattr(first,'dtype'):
        t1 = first.dtype
    else:
        t1 = type(first)
    if hasattr(second,'dtype'):
        t2 = second.dtype
    else:
        t2 = type(second)
    
    if t1 == type(None) and t2 == type(None):
        return True
    elif t1 == t2:
        return np.isclose(first,second,atol=margin)
    else:
        return False

def OxfordList(list):
    """Create a comma separated string from the strings provided with last comma trailed by 'and'."""
    if len(list)==0:
        return None
    elif len(list)==1:
        return str(list[0])
    elif len(list)==2:
        return ' and '.join([str(x) for x in list])
    else:
        return ', '.join([str(x) for x in list[:-1]])+', and ' +str(list[-1])



def convertToQxQy(sample,QPoints):
    """Convert a given list og QPoints to QxQy from UB matrix of sample

    Args:

        - sample (MJOLNIR.Sample.Sample): Sample from which the UB matrix is to be used

        - QPoints (list): List of HKL points to be converted

    Returns:

        - Q (list): List of QxQy points in same shape as provided


    """

    QPoints = np.asarray(QPoints)
    shape = QPoints.shape

    if len(shape)==1: # One point given as [h,k,l]
        if shape[0]!=3:
            raise AttributeError('Provided HKL point is not 3D. Received: {}'.format(QPoints))
        qx,qy,qz = np.einsum('ij,j->i',sample.orientationMatrix,QPoints)
    else:
        if shape[-1]!=3:
            raise AttributeError('Provided HKL point is not 3D. Received: {}'.format(QPoints))
        qx,qy,qz = np.einsum('ij,...j->i...',sample.orientationMatrix,QPoints)


    return np.array([qx,qy]).T

def convertToHKL(sample,QxQy):
    """Convert a given list og QPoints to QxQy from UB matrix of sample

    Args:

        - sample (MJOLNIR.Sample.Sample): Sample from which the UB matrix is to be used

        - QxQy (list): List of HKL points to be converted

    Returns:

        - HKL (list): List of QxQy points in same shape as provided


    """

    QxQy = np.asarray(QxQy)
    shape = QxQy.shape

    if len(shape)==1: # One point given as [h,k,l]
        if shape[0]!=2:
            raise AttributeError('Provided QxQy point is not 3D. Received: {}'.format(QxQy))
        QxQy = np.pad(QxQy, (0, 1), 'constant')
        H,K,L = np.einsum('ij,j->i',sample.orientationMatrixINV,QxQy)
    else:
        if shape[-1]!=2:
            raise AttributeError('Provided QxQy point is not 2D. Received: {}'.format(QxQy))
        Shape = np.asarray(shape)
        Shape[-1]=1
        z  = np.zeros(Shape)
        QxQy = np.concatenate([QxQy,z],axis=-1)
        H,K,L = np.einsum('ij,...j->i...',sample.orientationMatrixINV,QxQy)


    return np.array([H,K,L]).T

