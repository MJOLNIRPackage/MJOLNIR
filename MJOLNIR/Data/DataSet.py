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
import MJOLNIR.Data.DataFile
import MJOLNIR.Data.Sample
from MJOLNIR import _tools
from mpl_toolkits.axisartist.grid_helper_curvelinear import \
    GridHelperCurveLinear
from mpl_toolkits.axisartist import SubplotHost
import pytest
from scipy.ndimage import filters
import scipy.optimize
from scipy.spatial import Voronoi,ConvexHull,KDTree
from shapely.geometry import Polygon as PolygonS
from shapely.geometry import Point as PointS
from shapely.vectorized import contains
import time
import warnings

pythonVersion = sys.version_info[0]


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
            self.sample = self.convertedFiles[0].sample
        elif len(self.dataFiles)!=0:
            self.sample = self.dataFiles[0].sample

        # Add all other kwargs to settings
        for key in kwargs:
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
            [self._dataFiles.append(file) for file in correctDataFiles if file.type=='hdf']
            [self._convertedFiles.append(file) for file in correctDataFiles if file.type=='nxs']
        except Exception as e:
            raise(e)


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
        if np.sum(mask)==0:
            warnings.warn('Provided mask has no masked elements!')
        elif np.sum(mask)==self.I.size:
            warnings.warn('Provided mask masks all elements!')
        self._mask = mask
        for att in self.__dict__.values():
            if hasattr(att,'extractData'):
                att.mask = mask

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
        try:
            return self.dataFiles[index]
        except IndexError:
            raise IndexError('Provided index {} is out of bounds for DataSet with length {}.'.format(index,len(self)))

    def __len__(self):
        return len(self.dataFiles)
    
    def __iter__(self):
        self._index=0
        return self
    
    def __next__(self):
        if self._index >= len(self):
            raise StopIteration
        result = self.dataFiles[self._index]
        self._index += 1
        return result

    def next(self):
        return self.__next__()
    
    def append(self,item):
        try:
            correctDataFiles = isListOfDataFiles(item)
            [self._dataFiles.append(file) for file in correctDataFiles if file.type=='hdf']
            [self._convertedFiles.append(file) for file in correctDataFiles if file.type=='nxs']
        except Exception as e:
            raise(e)
        self._getData
    
    def __reversed__(self):
        return self[::-1]
    
    def __delitem__(self,index):
        try:
            del self.dataFiles[index]
        except IndexError:
            raise IndexError('Provided index {} is out of bounds for DataSet with length {}.'.format(index,len(self)))
        self._getData




    @_tools.KwargChecker()
    def convertDataFile(self,dataFiles=None,binning=8,saveLocation=None,saveFile=True):
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

            if saveFile:
                if not saveLocation is None:
                    if not os.path.isabs(saveLocation): # if full path is given
                        saveloc = saveLocation
                        if not saveLocation.split('.')[-1] == 'nxs':
                            if saveLocation[-1]!='/':
                                saveLocation+='/'
                            saveloc = saveLocation+rawfile.fileLocation.replace('.hdf','.nxs').split('/')[-1]
                        else:
                            saveloc = saveLocation
                    else: # pragma: no cover
                        if not saveLocation.split('.')[-1] == 'nxs': # is not covered as testing platform is to be used with relative paths
                            if saveLocation[-1]!='/':
                                saveLocation+='/'
                            saveloc = saveLocation+rawfile.fileLocation.replace('.hdf','.nxs').split('/')[-1]
                        else:
                            saveloc = saveLocation
                else:
                    saveloc = rawfile.fileLocation.replace('.hdf','.nxs')
                
                convFile.saveNXsqom(saveloc)
            
            #file.close()
            
            convertedFiles.append(convFile)
        self.convertedFiles = convertedFiles    
        self._getData()
            
    def _getData(self): # Internal method to populate I,qx,qy,energy,Norm and Monitor
        
        if len(self.convertedFiles)!=0:
            self.I,self.qx,self.qy,self.energy,self.Norm,self.Monitor,self.a3,self.a3Off,self.a4,self.a4Off,self.instrumentCalibrationEf, \
            self.instrumentCalibrationA4,self.instrumentCalibrationEdges,self.Ei,self.scanParameters,\
            self.scanParameterValues,self.scanParameterUnits,self.h,self.k,self.l = MJOLNIR.Data.DataFile.extractData(self.convertedFiles)
        else:
            self.I,self.Monitor,self.a3,self.a3Off,self.a4,self.a4Off,self.instrumentCalibrationEf, \
            self.instrumentCalibrationA4,self.instrumentCalibrationEdges,self.Ei,self.scanParameters,\
            self.scanParameterValues,self.scanParameterUnits = MJOLNIR.Data.DataFile.extractData(self.dataFiles)

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

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor, = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData()
        if rlu:
            qx,qy = np.einsum('ij,j...->i...',self.sample.RotMat,np.array([qx,qy]))
        pos=[qx,qy,energy]
        returnData,bins = binData3D(dx=dx,dy=dy,dz=dz,pos=pos,data=I,norm=Norm,mon=Monitor)

        return returnData,bins

    @_tools.KwargChecker()
    def cut1D(self,q1,q2,width,minPixel,Emin,Emax,rlu=True,plotCoverage=False,extend=True,dataFiles=None,constantBins=False):
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
        
        
        Returns:
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
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
            I,qx,qy,energy,Norm,Monitor, = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData()
        positions = np.array([qx,qy,energy])

        if rlu==True: # Recalculate H,K,L to qx
            q1,q2 = self.convertToQxQy([q1,q2])
            Data,[binpositionsTotal,orthopos,Earray] = cut1D(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,
                                                            minPixel=minPixel,Emin=Emin,Emax=Emax,plotCoverage=plotCoverage,
                                                            extend=extend,constantBins=constantBins)
            binpositionsTotal = np.concatenate([self.convertToHKL(binpositionsTotal[:,:2]),binpositionsTotal[:,-1].reshape(-1,1)],axis=1)
            orthopos = self.convertToHKL(orthopos)
            return Data,[binpositionsTotal,orthopos,Earray]
            

       
        return cut1D(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,minPixel=minPixel,
                     Emin=Emin,Emax=Emax,plotCoverage=plotCoverage,extend=extend,constantBins=constantBins)

    @_tools.KwargChecker(function=plt.errorbar,include=[_tools.MPLKwargs,'ticks','tickRound']) #Advanced KWargs checker for figures
    def plotCut1D(self,q1,q2,width,minPixel,Emin,Emax,rlu=True,ax=None,plotCoverage=False,extend=True,dataFiles=None,constantBins=False,**kwargs):  
        """Plotting wrapper for the cut1D method. Generates a 1D plot with bins at positions corresponding to the distance from the start point. 
        Adds the 3D position on the x axis with ticks.
        
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
            
            - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
            
            - kwargs: All other keywords will be passed on to the ax.errorbar method.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.

            - ticks (int): Number of tick marks to be used

            - tickRound (int): Decimals to be used when creating ticks

            - constantBins (bool): If True only bins of size minPixel is used (default False)
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - binCenter (3D array): Array containing the position of the bin centers of size (n,3)
            
            - binDistance (array): Distance from centre of bins to start position.
        """
        
        
        
        D,P = self.cut1D(q1=q1,q2=q2,width=width,minPixel=minPixel,Emin=Emin,Emax=Emax,\
        plotCoverage=plotCoverage,extend=extend,rlu=rlu,dataFiles=dataFiles,constantBins=constantBins)

        INT = np.divide(D[0]*D[3],D[1]*D[2])
        INT_err = np.divide(np.sqrt(D[0])*D[3],D[1]*D[2])
        
        
        binCenter = 0.5*(P[0][:-1]+P[0][1:])
        num=len(binCenter)
        
       
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

        xvalues = np.round(np.linspace(0,num-1,ticks)).astype(int)
        my_xticks=[]
        for i in xvalues:
            my_xticks.append('\n'.join(map(lambda x:('{:.'+str(tickRound)+'f}').format(x),[np.round(binCenter[i,j],tickRound) for j in range(len(binCenter[i]))])))
        
        binDistance = np.linalg.norm(binCenter[:,:-1]-P[0][0,:-1],axis=1)
        
        if ax is None:
            plt.figure()
            ax = plt.gca()
        
        ax.errorbar(binDistance,INT,yerr=INT_err,**kwargs)

        ax.set_xticks(binDistance[xvalues])
        ax.set_xticklabels(my_xticks,fontsize=10, multialignment="center",ha="center")
        
        def calculateIndex(binDistance,x):
            return np.argmin(np.abs(binDistance-x))

        ax.calculateIndex = lambda x: calculateIndex(binDistance,x)
        
        if rlu==False:
            ax.set_xlabel('$Q_x/A$\n$Q_y/A$\nE/meV', fontsize=8)
            def format_coord(x,y,ax,binCenter):# pragma: no cover
                index = ax.calculateIndex(x)
                qx,qy,E = binCenter[index]
                return  "qx = {0:.3e}, qy = {1:.3e}, E = {2:.3f}, I = {3:0.4e}".format(qx,qy,E,y)
        else:
            def format_coord(x,y,ax,binCenter):# pragma: no cover
                index = ax.calculateIndex(x)
                h,k,l,E = binCenter[index]
                return  "H = {0:.3e}, K = {1:.3e}, L = {2:.3e}, E = {3:.3f}, I = {4:0.4e}".format(h,k,l,E,y)
            ax.set_xlabel('$Q_h/A$\n$Q_k/A$\n$Q_l/A$\nE/meV', fontsize=8)

        
        def onclick(event,ax,DataList):# pragma: no cover
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
                index = ax.calculateIndex(x)
            
                cts = int(DataList[0][index])
                Mon = int(DataList[1][index][0])
                Norm = float(DataList[2][index][0])
                NC = int(DataList[3][index][0])
                printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
            print(printString)


        ax.xaxis.set_label_coords(1.15, -0.025)
        ax.set_ylabel('Int [arb]')
        plt.tight_layout()
        
        ax.format_coord = lambda x,y: format_coord(x,y,ax,binCenter)
        ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,D))
        return ax,D,P,binCenter,binDistance

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
            
            - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
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

        else: 
            #dataFiles = isListOfDataFiles(dataFiles)
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData()
        
        positions = np.array([qx,qy,energy])
        if rlu==True: # Recalculate H,K,L to qx
            q1,q2 = self.convertToQxQy([q1,q2])

        return cutQE(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,
                    minPixel=minPixel,EnergyBins=EnergyBins,extend=extend,constantBins=constantBins)

 
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
            
            - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
            - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (n * 3D arrays): n instances of center positions for the bins.

            - binDistance (n arrays): n instances of arrays holding the distance in q to q1.
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

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor = DS.I.extractData(),DS.qx.extractData(),DS.qy.extractData(),DS.energy.extractData(),DS.Norm.extractData(),DS.Monitor.extractData()
            
        if rlu==True: # Recalculate H,K,L to qx
            q1,q2 = self.convertToQxQy([q1,q2])

        positions = np.array([qx,qy,energy])
        return plotCutQE(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,
                        minPixel=minPixel,EnergyBins=EnergyBins,rlu=rlu,ax = ax,constantBins=constantBins,**kwargs)

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
            
            - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
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

    @_tools.KwargChecker(function=plt.pcolormesh,include=['vmin','vmax'])
    def plotCutPowder(self,EBinEdges,qMinBin=0.01,ax=None,dataFiles=None,constantBins=False,**kwargs):
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

            - kwargs: All other keywords will be passed on to the ax.pcolormesh method.
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).

        """

        DataList,qbins = self.cutPowder(EBinEdges=EBinEdges,qMinBin=qMinBin,dataFiles=dataFiles,constantBins=constantBins)
        intensity,monitorCount,Normalization,NormCount = DataList
        warnings.simplefilter('ignore')
        Int = [np.divide(Int*NC,mC*N) for Int,NC,mC,N in zip(intensity,NormCount,monitorCount,Normalization)]
        warnings.simplefilter('once')
        eMean = 0.5*(EBinEdges[:-1]+EBinEdges[1:])
        
        if ax is None:
            plt.figure()
            ax = plt.gca()
        pmeshs = []
        
        for i in range(len(EBinEdges)-1):
            if len(Int[i])==0:
                continue
            pmeshs.append(ax.pcolormesh(qbins[i],[EBinEdges[i],EBinEdges[i+1]],Int[i].reshape((len(qbins[i])-1,1)).T,**kwargs))
        def calculateIndex(x,y,eMean,qBin):# pragma: no cover
            EIndex = np.argmin(np.abs(y-eMean))
            qIndex = np.argmin(np.abs(x-0.5*(qBin[EIndex][:-1]+qBin[EIndex][1:])))
            return EIndex,qIndex
        ax.calculateIndex = lambda x,y: calculateIndex(x,y,eMean,qbins)

        def format_coord(x,y,ax,Int,qbins):# pragma: no cover
            EIndex,qIndex = ax.calculateIndex(x,y)
            Intensity = Int[EIndex][qIndex] 
            return  "|q| = {0:.3f}, E = {1:.3f}, I = {2:0.4e}".format(qbins[EIndex][qIndex],eMean[EIndex],Intensity)
            
        ax.format_coord = lambda x,y: format_coord(x,y,ax,Int,qbins)
        ax.set_xlabel('|q| [1/A]')
        ax.set_ylabel('E [meV]')
        
        def set_clim(VMin,VMax,pmesh):
            for pm in pmeshs:
                pm.set_clim(VMin,VMax)

        ax.set_clim = lambda VMin,VMax: set_clim(VMin,VMax,pmeshs)
        
        def onclick(event,ax,DataList):# pragma: no cover
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
                Eindex,index = ax.calculateIndex(x,y)
                cts = int(DataList[0][Eindex][index])
                Mon = int(DataList[1][Eindex][index])
                Norm = float(DataList[2][Eindex][index])
                NC = int(DataList[3][Eindex][index])
                printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
                print(printString)

        if not 'vmin' in kwargs or not 'vmax' in kwargs:
            minVal = np.min(np.concatenate(Int))
            maxVal = np.max(np.concatenate(Int))
            ax.set_clim(minVal,maxVal)
        ax.pmeshs = pmeshs
        ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,DataList))
        return ax,[intensity,monitorCount,Normalization,NormCount],qbins

    @_tools.KwargChecker()
    @_tools.overWritingFunctionDecorator(RLUAxes.createRLUAxes)
    def createRLUAxes(*args,**kwargs): # pragma: no cover
        raise RuntimeError('This code is not meant to be run but rather is to be overwritten by decorator. Something is wrong!! Should run {}'.format(RLUAxes.createRLUAxes))

    @_tools.KwargChecker()
    @_tools.overWritingFunctionDecorator(RLUAxes.createQEAxes)
    def createQEAxes(*args,**kwargs): # pragma: no cover
        raise RuntimeError('This code is not meant to be run but rather is to be overwritten by decorator. Something is wrong!! Should run {}'.format(RLUAxes.createQEAxes))
    
    #@_tools.KwargChecker(function=plt.pcolormesh,include=['vmin','vmax','colorbar','zorder'])
    def plotQPlane(self,EMin=None,EMax=None,EBins=None,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=False,ax=None,rlu=True,dataFiles=None,xScale=1.0,yScale=1.0,**kwargs):
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

        else: 
            DS = DataSet(convertedFiles = dataFiles)
            I,qx,qy,energy,Norm,Monitor, = DS.I,DS.qx,DS.qy,DS.energy,DS.Norm,DS.Monitor
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
            positions = np.array([qx,qy])
            qx,qy = np.einsum('ij,j...->i...',self.sample.RotMat,positions)
            
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
            ax.get_figure().colorbar(ax.pmeshs[0],pad=0.1)

        ax.set_clim(vmin,vmax)
        if _3D:
            minEBins = np.min(EBins)
            maxEBins = np.max(EBins)
            if not np.isclose(minEBins,maxEBins):
                ax.set_zlim(minEBins,maxEBins)
            else:
                ax.set_zlim(minEBins-0.1,maxEBins+0.1)
        else:
            def onclick(ax, event,Qx,Qy,data): # pragma: no cover
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

                print(printString)
            ax.cid = ax.figure.canvas.mpl_connect('button_press_event', lambda x: onclick(ax,x,Qx,Qy,[intensity,monitorCount,Normalization,NormCount]))
            
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
        
            - width (float): Width of the cut in 1/A (default 0.1).
            
            - minPixel (float): Minimal size of binning along the cutting directions. Points will be binned if they arecloser than minPixel (default=0.01)
        
            - rlu (bool): If True, provided QPoints are interpreted as (h,k,l) otherwise as (qx,qy), (default True).
        
            - dataFiles (list): List of dataFiles to cut. If none, the ones in the object will be used (default None).

            - constantBins (bool): If True only bins of size minPixel is used (default False)
        
        .. warning::
            The way the binning works is by extending the end points with 0.5*minPixel, but the method sorts away points not between the two Q points given and thus the start and end
            bins are only half filled. This might result in discrepancies between a single cut and the same cut split into different steps. Further, splitting lines into sub-cuts 
            forces a new binning to be done and the bin positions can then differ from the case where only one cut is performed.

        
        Returns: m = Q points, n = energy bins
                
            - Data list (m * n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
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
            sample =self.convertedFiles[0].sample
            positions = self.convertToQxQy(QPoints)
            
        elif rlu==False: # RLU is false
            positions = QPoints
            if QPoints.shape[1]!=2:
                raise AttributeError('Provide Q list is not 2 dimensional, should have shape (n,2) in QxQy mode but got shape {}.'.format(QPoints.shape))
        else:
            raise AttributeError('Given Q mode not understood. Got {} but must be either "RLU", "HKL" or "QxQy"')

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
        for pStart,pStop,w,mP,EB in zip(positions,positions[1:],width,minPixel,EnergyBins):
            _DataList,_BinList,_centerPosition,_binDistance = self.cutQE(q1=pStart,q2=pStop,width=w,minPixel=mP,EnergyBins=EB,rlu=False,
                                                                         dataFiles=dataFiles,extend=False,constantBins=constantBins)
            DataList.append(_DataList)
            if rlu:
                UB2D = self.sample.convertHKLINV # Matrix to calculate HKL from Qx,Qy
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
            
        return np.array(DataList),np.array(BinList),np.array(centerPosition),np.array(binDistance)


    @_tools.KwargChecker(include=np.concatenate([_tools.MPLKwargs,['vmin','vmax','log','ticks','seperatorWidth','tickRound','plotSeperator','cmap','colorbar','edgecolors']]))
    def plotCutQELine(self,QPoints,EnergyBins,width=0.1,minPixel=0.01,rlu=True,ax=None,dataFiles=None,constantBins=False,**kwargs):
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

        Return:  m = Q points, n = energy bins
            
            - ax: matplotlib axis in which the data is plotted
            
            - Data list (m * n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
                
            - Bin list (m * n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (m * n * 3D arrays): n instances of center positions for the bins.

            - binDistance (m * n arrays): n instances of arrays holding the distance in q to q1.

        .. note::
            
            The ax.set_clim function is created to change the colour scale. It takes inputs vmin,vmax. This function does however not work in 3D....

        """
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
            try:
                ax.name
            except:
                _3D = False
            else:
                if ax.name == '3d':
                    _3D = True
                else:
                    _3D = False
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
                
            if not 'log' in kwargs:
                log = False
            else:
                log = kwargs['log']
                del kwargs['log']
                
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
            idmax = len(DataList) # Number of segments
            edgeQDistance = []
            actualEnergy = []
            for BLT in BinListTotal:
                localE = []
                for BLTSub in BLT:
                    localE.append(BLTSub[2][0])
                
                localE.append(BLT[-1][2][1])
                
                actualEnergy.append(localE)
                
            for segID in range(idmax): # extract relevant data for current segment
                [intensityArray,monitorArray,normalizationArray,normcountArray] = DataList[segID]#[i] for i in range(4)]
                BinList = BinListTotal[segID]
                centerPosition = centerPositionTotal[segID]
                binDistance = binDistanceTotal[segID]
                q1 = positions[segID]
                q2 = positions[segID+1]
                if rlu:
                    q1 = np.dot(self.sample.convertHKLINV,q1)
                    q2 = np.dot(self.sample.convertHKLINV,q2)
                dirvec = np.array(q2)-np.array(q1)
                
                leftEdgeBins = np.array([np.dot(BL[0][0,:len(q1)],dirvec) for BL in BinList])
            
                leftEdgeIndex = np.argmin(leftEdgeBins)

                binedges = [np.linalg.norm(BL[0][:,:len(q1)]-BinList[leftEdgeIndex][0][0,:len(q1)],axis=1) for BL in BinList]
                
                binenergies = [BL[2] for BL in BinList]

                edgeQDistanceLocal = []
                for BL in BinList:
                    p = BL[0][:,:len(q1)] - q1
                    q = np.dot(p, dirvec)
                    if not (np.sort(q) == q).all():
                        raise RuntimeError("edgeQDistance[{}] is not sorted".format(iE))
                    edgeQDistanceLocal.append(q)
                edgeQDistance.append(edgeQDistanceLocal)

                # Calculate scattering intensity
                warnings.simplefilter('ignore')
                Int = [np.divide(iA*ncA,mA*nA) for iA,ncA,mA,nA in zip(intensityArray,normcountArray,monitorArray,normalizationArray)]
                warnings.simplefilter('once')
                IntTotal.append(Int)
                
                # plot intensity with provied kwargs
                for intensity,binEdge,binEnergy in zip(Int,binedges,binenergies):
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
                
                
                qstart= centerPosition[minimalDistanceIDEnergy][minimalDistanceID][:len(q1)]
                qstop = centerPosition[maximalDistanceIDEnergy][maximalDistanceID][:len(q1)]
                # Prepare the calculation of tick markers
                direction = (qstop-qstart)
                distanceChange = np.max(binDistance[maximalDistanceIDEnergy])-np.min(binDistance[minimalDistanceIDEnergy])
                if rlu:
                    qstartQ,qstopQ = [np.dot(self.sample.convertHKL,x) for x in [qstart,qstop]]
                    dirLen = np.linalg.norm(direction)
                    dirLenQ = np.linalg.norm(qstopQ-qstartQ)
                    distanceChange*=dirLen/dirLenQ
                
                binminmaxList = np.linspace(0,1,100)
                ticksCurrentSeg = ticksInSegments[segID]
                

                num=len(binminmaxList)
                
                if segID==idmax-1: ## Find best bins for tick marking
                    xvalues = np.round(np.linspace(0,num-1,ticksCurrentSeg+1)).astype(int) 
                else:
                    xvalues = np.round(np.linspace(0,num-1-int(num/ticksCurrentSeg),ticksCurrentSeg)).astype(int)

                my_xticks=[]
                for i in xvalues:
                    positionValues = binminmaxList[i]*direction+qstart
                    my_xticks.append('\n'.join([('{:.'+str(tickRound)+'f}').format(x+0.0) for x in positionValues]))
                
                xticks.append(binminmaxList[xvalues]*distanceChange +offset[-1]) # binDistanceAll[xvalues]
                xticklabels.append(my_xticks)
                offset.append(offset[-1]+np.max([np.max(binedge) for binedge in binedges]))

            if plotSeperator == True: # plot last seperator
                plt.plot([offset,offset],[np.min(EnergyBins[-1]),np.max(EnergyBins[-1])],'k',linewidth=seperatorWidth)
            
            ax.set_xticks(np.concatenate(xticks))
            ax.set_xticklabels(np.concatenate(xticklabels),fontsize=8, multialignment="center",ha="center")
            ax.EnergyBins = EnergyBins
            if rlu==True:
                ax.set_xlabel('$H$\n$K$\n$L$', fontsize=8)
            else:
                ax.set_xlabel('$Q_x/A$\n$Q_y/A$', fontsize=8)
            
            ax.xaxis.set_label_coords(1.15, -0.025) 
            ax.set_ylabel('E [meV]')
            if colorbar:
                ax.colorbar = ax.get_figure().colorbar(pmeshs[0],pad=0.1,format='%.2E')
            #plt.tight_layout()
            if 'pmeshs' in ax.__dict__:
                ax.pmeshs = np.concatenate([ax.pmeshs,pmeshs],axis=0)
            else:
                ax.pmeshs = pmeshs
        
        
            def set_clim(vmin,vmax):
                for pm in pmeshs:
                    pm.set_clim(vmin,vmax)
                    
            ax.set_clim = set_clim#lambda VMin,VMax: [pm.set_clim(VMin,VMax) for pm in pmeshs] 
            
            if not 'vmin' in kwargs:
                if log==True:
                    vmin = np.nanmin(np.log10(np.concatenate(Int)+1e-20))
                else:
                    vmin = np.nanmin(np.concatenate(Int))
            else:
                vmin = kwargs['vmin']
            if not 'vmax' in kwargs:
                if log==True:
                    vmax = np.nanmax(np.log10(np.concatenate(Int)+1e-20))
                else:
                    vmax = np.nanmax(np.concatenate(Int))
            else:
                vmax = kwargs['vmax']
            
            ax.set_clim(vmin,vmax)
            
            # Create mouse over function
            def format_coord(x,y,edgeQDistance,centerPos,EnergyBins,Int,rlu,offset,self):# pragma: no cover
                val = calculateIndex(x,y,offset,EnergyBins,edgeQDistance,textReturn=True)
                if type(val)==str:
                    return val
                segID,Eindex,index = val

                Intensity = Int[segID][Eindex][index][0]
                if rlu==False: 
                    qx, qy, E = centerPos[segID][Eindex][index]
                    return "qx = {0:.3f}, qy = {1:.3f}, E = {2:.3f}, I = {3:.3e}".format(qx+0.0,qy+0.0,E,Intensity)
                else:
                    H,K,L, E = centerPos[segID][Eindex][index]
                    return "h = {0:.3f}, h = {1:.3f}, l = {2:.3f}, E = {3:.3f}, I = {4:.3e}".format(H+0.0,K+0.0,L+0.0,E,Intensity)

            def calculateIndex(x,y,offset,EnergyBins,edgeQDistance,textReturn): # pragma: no cover
                if x<offset[0] or x>offset[-1]:
                    if textReturn == True:
                        return "x out of range: {:.3}".format(x)
                    else:
                        return -1,-1,-1
                segID = (np.arange(len(offset)-1)[[x>offStart and x<offStop for offStart,offStop in zip(offset,offset[1:])]])[0]
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
            ax.calculateIndex = lambda x,y: calculateIndex(x,y,offset,actualEnergy,edgeQDistance,textReturn=False)#EnergyBins
            def onclick(event,ax,DataList): # pragma: no cover
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
                    if not np.any([x==-1 for x in [segID,Eindex,index]]):
                        cts = int(DataList[segID][0][Eindex][index])
                        Mon = int(DataList[segID][1][Eindex][index][0])
                        Norm = float(DataList[segID][2][Eindex][index][0])
                        NC = int(DataList[segID][3][Eindex][index][0])
                        printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
                    print(printString)
            ax.format_coord = lambda x,y: format_coord(x,y,edgeQDistance,centerPositionTotal,actualEnergy,IntTotal,rlu,offset,self)#EnergyBins
            ax.figure.canvas.mpl_connect('button_press_event',lambda event:onclick(event,ax,DataList))
        else: # pragma: no cover
            # TODO: Make test!!!
            import matplotlib.colors

            if not 'vmin' in kwargs:
                if log==True:
                    vmin = np.nanmin(np.log10(np.concatenate(Int)+1e-20))
                else:
                    vmin = np.nanmin(np.concatenate(Int))
            else:
                vmin = kwargs['vmin']
            if not 'vmax' in kwargs:
                if log==True:
                    vmax = np.nanmax(np.log10(np.concatenate(Int)+1e-20))
                else:
                    vmax = np.nanmax(np.concatenate(Int))
            else:
                vmax = kwargs['vmax']


            ax.norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            if not 'cmap' in kwargs:
                from matplotlib.colors import ListedColormap
                cmap = plt.cm.coolwarm
            else:
                cmap = kwargs['cmap']
                kwargs = _tools.without_keys(dictionary=kwargs,keys='cmap')
            sfp = []
            for bins,datlist in zip(BinListTotal,DataList):
                energies = len(bins)
                
                energyEdges = np.array([bins[idx][2] for idx in range(energies)])
                ELength = np.array([len(x[0][:,0]) for x in bins])
                ELengthCummu = np.concatenate([[0],np.cumsum(ELength)],axis=0)
                H = np.concatenate([bins[idx][0][:,0] for idx in range(energies)],axis=0)
                K = np.concatenate([bins[idx][0][:,1] for idx in range(energies)],axis=0)
                L = np.concatenate([bins[idx][0][:,2] for idx in range(energies)],axis=0)
                P0,P1 = self.sample.calculateHKLToQxQy(H,K,L)
                P0,P1 = np.einsum('mj,j...->m...',self.sample.RotMat,[P0,P1])
                
                Data = np.array([np.concatenate(x,axis=0) for x in datlist]).squeeze()
                INT = np.divide(Data[0]*Data[3],Data[2]*Data[1])
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
                    normColor = INT[IntCommu[E]:IntCommu[E+1]].reshape(-1,1).repeat(2,axis=1)
                    normColor = np.concatenate([normColor,[normColor[-1]]])
                    if len(normColor)!=X.shape[0]:
                        continue
                    sf = ax.plot_surface(X,Y,Z,facecolors=cmap(ax.norm(normColor)),rstride=1,cstride=1,shade=False)
                    sf.value = normColor
                    # TODO: make colorlimits changeable after creation of axis
                    sfp.append(sf)

                    
            pmeshs = np.array(sfp).flatten()
            if 'pmeshs' in ax.__dict__:
                ax.pmeshs = np.concatenate([ax.pmeshs,pmeshs],axis=0)
            else:
                ax.pmeshs = pmeshs
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
    def cut1DE(self,E1,E2,q,rlu=True,width=0.02, minPixel = 0.1, dataFiles = None,constantBins=False):
        """Perform 1D cut through constant Q point returning binned intensity, monitor, normalization and normcount. The width of the cut is given by 
        the width attribute.
        
        .. note::
            Can only perform cuts for a constant energy plane of definable width.
        
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
            
        Returns:
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
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

            projectOrientation = np.dot(sample.orientationMatrix,q)
            Q = self.convertToQxQy(q)
            # Add width/2 to point in 4 direction and make width mean of the converted positions
            diffQ = projectOrientation.reshape(3,1)+np.array([[width*0.5,0,0],[-width*0.5,0,0],[0,width*0.5,0],[0,-width*0.5,0]]).T
            dist = [np.array(sample.tr(diffQ[0][i],diffQ[1][i])) for i in range(4)]-Q
            width = np.mean(np.linalg.norm(dist,axis=1))
        else: # Do nothing
            Q = np.array(q)


        return cut1DE(positions = positions, I=I, Norm=Norm,Monitor=Monitor,E1=E1,E2=E2,q=Q,width=width,minPixel=minPixel,constantBins=constantBins)

    @_tools.KwargChecker(function=createRLUAxes)
    def View3D(self,dQx,dQy,dE,rlu=True, log=False,grid=False,axis=2,**kwargs):
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
        warnings.simplefilter('ignore')
        Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
        warnings.simplefilter('once')

        Viewer = Viewer3D.Viewer3D(Data=Data,bins=bins,axis=axis,ax=axes,grid=grid,log=log)
        return Viewer

    def convertToQxQy(self,HKL):
        """Convert array or vector of HKL point(s) to corresponding Qx and QY

        Args:

            - HKL (array): array or vector of HKL point(s)

        Returns

            - Q (array): Converted HKL points in Qx QY of un-rotated coordinate system.
        """
        return convertToQxQy(self.sample,HKL)

    def convertToHKL(self,QxQy):
        """Convert array or vector of QxQy point(s) to corresponding HKL

        Args:

            - QxQy (array): array or vector of QxQy point(s)

        Returns

            - HKL (array): Converted QxQy points in HKL
        """
        return convertToHKL(self.sample,QxQy)

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
        
        - plotCoverage (bool): If True, generates plot of all points in the cutting plane and adds bounding box of cut (default False).

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
        return [np.array(np.array([])),np.array([]),np.array([]),np.array([])],[np.array([]),np.array([]),[Emin,Emax]]
        #raise AttributeError('No points are within the provided energy limits.')

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
   
    if plotCoverage: # pragma: no cover
        plt.figure()
        plt.scatter(positions2D[0],positions2D[1],s=0.5)
        plt.plot([binpositions[0][0]+orthopos[0][0],binpositions[-1][0]+orthopos[0][0]],[binpositions[0][1]+orthopos[0][1],binpositions[-1][1]+orthopos[0][1]],c='k')
        plt.plot([binpositions[0][0]+orthopos[1][0],binpositions[-1][0]+orthopos[1][0]],[binpositions[0][1]+orthopos[1][1],binpositions[-1][1]+orthopos[1][1]],c='k')
        for i in [0,-1]:
            plt.plot([binpositions[i][0]+orthopos[0][0],binpositions[i][0]+orthopos[1][0]],[binpositions[i][1]+orthopos[0][1],binpositions[i][1]+orthopos[1][1]],c='k')
        for binPos in binpositions:#i in range(len(binpositions)):
            plt.plot([binPos[0]+orthopos[0][0],binPos[0]+orthopos[1][0]],[binPos[1]+orthopos[0][1],binPos[1]+orthopos[1][1]],c='k',linewidth=0.5)
        if extend==False:
            plt.scatter(positions2D[0][insideQ][insideWidth],positions2D[1][insideQ][insideWidth],s=0.5)
        else:
            plt.scatter(positions2D[0][insideWidth],positions2D[1][insideWidth],s=0.5)
        ax = plt.gca()
        ax.set_aspect('equal', 'datalim')
        ax.set_xlabel('Qx [1/A]')
        ax.set_ylabel('Qy [1/A]')
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
    MonitorCount=  np.histogram(Energies,bins=bins,weights=np.array(Monitor[allInside].flatten(),dtype=np.int64))[0] # Need to change to int64 to avoid overflow
    Normalization= np.histogram(Energies,bins=bins,weights=Norm[allInside].flatten())[0]
    

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
        
        - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
        
        - qbins (n arrays): n arrays holding the bin edges along the length of q

    """
    qx,qy,energy = positions

    q = np.linalg.norm([qx,qy],axis=0)
    intensity = []
    monitorCount = []
    Normalization = []
    NormCount = []
    qbins = []
    
    #for i in range(len(EBinEdges)-1):
    for binStart,binEnd in zip(EBinEdges,EBinEdges[1:]):
        e_inside = np.logical_and(energy>binStart,energy<=binEnd)
        q_inside = q[e_inside]
        if constantBins==False:
            qbins.append(np.array(_tools.binEdges(q_inside,tolerance=qMinBin)))
        else:
            Min,Max = _tools.minMax(q_inside)
            qbins.append(np.arange(Min,Max+0.5*qMinBin,qMinBin))
            
        intensity.append(np.histogram(q_inside,bins=qbins[-1],weights=I[e_inside].flatten())[0].astype(I.dtype))
        monitorCount.append(np.histogram(q_inside,bins=qbins[-1],weights=Monitor[e_inside].flatten())[0].astype(Monitor.dtype))
        Normalization.append(np.histogram(q_inside,bins=qbins[-1],weights=Norm[e_inside].flatten())[0].astype(Norm.dtype))
        NormCount.append(np.histogram(q_inside,bins=qbins[-1],weights=np.ones_like(I[e_inside]).flatten())[0].astype(I.dtype))
    
    return [intensity,monitorCount,Normalization,NormCount],qbins


@_tools.KwargChecker()
def cutQE(positions,I,Norm,Monitor,q1,q2,width,minPixel,EnergyBins,extend=True,constantBins=False):
    """Cut data into maps of q and intensity between two q points and given energies. This is performed by doing consecutive constant energy planes.

    Args:

        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.
        
        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - q1 (2D array): Start position of cut in format (qx,qy).
        
        - q2 (2D array): End position of cut in format (qx,qy).
        
        - width (float): Full width of cut in q-plane.
        
        - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.

        - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

    Kwargs:

        - extend (bool): Whether or not the cut from q1 to q2 is to be extended throughout the data (default true)

        - constantBins (bool): If True only bins of size minPixel is used (default False)

    Returns:
        
        - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
        
        - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
        - center position (n * 3D arrays): n instances of center positions for the bins.

        - binDistance (n arrays): n instances of arrays holding the distance in q to q1.

    """
    intensityArray = []
    monitorArray = []
    normalizationArray = []
    normcountArray = []
    centerPos = []
    returnpositions = []
    binDistance = []

    dirvec = np.array(q2) - np.array(q1)
    dirvec /= np.linalg.norm(dirvec)

    for i in np.arange(len(EnergyBins)-1):
        [intensity,MonitorCount,Normalization,normcounts],position = cut1D(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,
                                                                            width=width,minPixel=minPixel,Emin=EnergyBins[i],Emax=EnergyBins[i+1],
                                                                            plotCoverage=False,extend=extend,constantBins=constantBins)
        if len(intensity)==0:
            continue
        returnpositions.append(position)
        intensityArray.append(intensity)
        monitorArray.append(MonitorCount)
        normalizationArray.append(Normalization)
        normcountArray.append(normcounts)

        thisCenterPos = 0.5*(position[0][:-1]+position[0][1:])
        centerPos.append(thisCenterPos)
        thisBinDistance = np.dot(thisCenterPos[:,:2] - q1, dirvec)
        binDistance.append(thisBinDistance)

    return [intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance

@_tools.KwargChecker(function=plt.errorbar,include=[_tools.MPLKwargs,'vmin','vmax'])
def plotCutQE(positions,I,Norm,Monitor,q1,q2,width,minPixel,EnergyBins,rlu=True,ax = None,constantBins=False,**kwargs):
    """Plotting wrapper for the cutQE method. Generates a 2D intensity map with the data cut by cutQE. 
    
    .. warning::
        Deprecated! Instead use the plotCutQELine tool with only two q points

    .. note::
        Positions shown in tool tip reflect the closes bin center and are thus limited to the area where data is present.

    Args:

        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.

        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array    

        - q1 (2D array): Start position of cut in format (qx,qy).
        
        - q2 (2D array): End position of cut in format (qx,qy).
        
        - width (float): Full width of cut in q-plane.
        
        - minPixel (float): Minimal size of binning along the cutting direction. Points will be binned if they are closer than minPixel.

        - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

    Kwargs:
        
        - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.

        - rlu (bool): If True, data is plotted in RLU. 
        
        - kwargs: All other keywords will be passed on to the ax.errorbar method.

        - constantBins (bool): If True only bins of size minPixel is used (default False)
    
    Returns:
        
        - ax (matplotlib axis): Matplotlib axis into which the plot was put.
        
        - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
        
        - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
        - center position (n * 3D arrays): n instances of center positions for the bins.

        - binDistance (n arrays): n instances of arrays holding the distance in q to q1.
    """
    warnings.warn("The plotCutQE is being deprecated. Use instead the plotCutQELine width only two QPoints",DeprecationWarning)

    [intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance = cutQE(positions=positions,I=I,Norm=Norm,Monitor=Monitor,q1=q1,q2=q2,width=width,minPixel=minPixel,
                                                                                                                    EnergyBins=EnergyBins,constantBins=constantBins)

    if len(returnpositions) < 1:
        raise RuntimeError("Expect at least one slice in energy dimension")

    q1 = np.array(q1)
    q2 = np.array(q2)

    dirvec = q2 - q1
    dirvec /= np.linalg.norm(dirvec)

    edgeQDistance = []
    for iE in range(len(EnergyBins)-1):
        p = returnpositions[iE][0][:,:2] - q1
        q = np.dot(p, dirvec)
        if not (np.sort(q) == q).all():
            raise RuntimeError("edgeQDistance[{}] is not sorted".format(iE))
        edgeQDistance.append(q)

    binEnergies = [x[2] for x in returnpositions]
    warnings.simplefilter('ignore')
    Int = [ np.divide( iA * ncA, mA * nA ) for iA,ncA,mA,nA in zip(intensityArray,normcountArray,monitorArray,normalizationArray) ]
    warnings.simplefilter('once')

    if ax is None:
        plt.figure()
        ax = plt.gca()

    pmeshs = []
    for intensity,edgeQ,binE in zip(Int,edgeQDistance,binEnergies):# in range(len(Int)):
        pmeshs.append(ax.pcolormesh(edgeQ, binE, intensity.T, **kwargs))

    def set_clim(VMin,VMax,pmeshs):
        [pm.set_clim(VMin,VMax) for pm in pmeshs]
    ax.set_clim = lambda VMin,VMax: set_clim(VMin,VMax,pmeshs)
    
    if not 'vmin' in kwargs or not 'vmax' in kwargs:
        minVal = np.nanmin(np.concatenate(Int))
        maxVal = np.nanmax(np.concatenate(Int))
        ax.set_clim(minVal,maxVal)

    def format_coord(x,y,edgeQDistance,centerPos,Int):# pragma: no cover
        if len(EnergyBins) < 2:
            return "len(EnergyBins) < 2"
        if y < EnergyBins[0] or y >= EnergyBins[-1]:
            return "E out of range {:.3}  {}..{}".format(y, EnergyBins[0], EnergyBins[-1])
        Eindex = EnergyBins.searchsorted(y) - 1
        if len(edgeQDistance[Eindex]) < 2:
            raise RuntimeError("len(edgeQDistance[{}]) < 2".format(Eindex))
        if x < edgeQDistance[Eindex][0] or x >= edgeQDistance[Eindex][-1]:
            return "x out of range: {:.3}".format(x)
        index = edgeQDistance[Eindex].searchsorted(x) - 1
        qx, qy, E = centerPos[Eindex][index]
        Intensity = Int[Eindex][index][0]
        return "qx = {0:.3f}, qy = {1:.3f}, E = {2:.3f}, I = {3:.3e}".format(qx, qy, E, Intensity)

    xtick_positions = []
    xtick_labels = []
    iE = 0
    m = len(centerPos[iE])
    for n in np.linspace(0, m-1, 4):
        i = int(round(n))
        xtick_positions.append(binDistance[iE][i])
        q = centerPos[iE][i][:2]
        xtick_labels.append("{0:.3f}\n{1:.3f}".format(q[0], q[1]))

    ax.format_coord = lambda x,y: format_coord(x,y,edgeQDistance,centerPos,Int)
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=8, multialignment="center", ha="center")
    ax.set_xlabel('$Q_h/A$\n$Q_k/A$', fontsize=8)
    ax.xaxis.set_label_coords(1.15, -0.025)
    ax.get_figure().colorbar(pmeshs[0],pad=0.1)
    ax.set_ylabel('E [meV]')
    plt.tight_layout()
    ax.pmeshs = pmeshs
    return ax,[intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance





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
            
    
    try:
        numFiles = len(files)
    except:
        numFiles = 1
        files = [files]

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
            ax[counter].get_figure().colorbar(ax[counter].collections[0], ax=ax[counter],format=ticker.FuncFormatter(fmt))
            
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

@_tools.KwargChecker()
def voronoiTessellation(points,plot=False,Boundary=False,numGroups=False):
    """Generate individual pixels around the given datapoints.

    Args:

        - points (list of list of points): Data points to generate pixels in shape [files,XY,N] i.e. [1,2,N] for one file with N points

    Kwargs:

        - plot (bool): If True, method plots pixels created with green as edge bins and red as internal (default False)

        - Boundary (list of Polygons): List of Shapely polygons constituting the boundaries (Default False)


    """

    if numGroups == False:
        numGroups = len(points)

    if Boundary==False:
        BoundPoly= [convexHullPoints(np.array(points[i][0]).flatten(),np.array(points[i][1]).flatten()) for i in range(numGroups)]
    else:
        BoundPoly = Boundary#[PolygonS(x.T) for x in Boundary]

    if numGroups == 1:

        combiPoly = BoundPoly[0]
        pointsX = np.array([points[0][0].flatten()])[0]
        pointsY = np.array([points[0][1].flatten()])[0]
    else: # Combine all files
        combiPoly = BoundPoly[0].union(BoundPoly[1])
        for i in range(len(BoundPoly)-2):
            combiPoly = combiPoly.union(BoundPoly[i+2])
        if Boundary==False:
            pointsX = np.concatenate([points[i][0].flatten() for i in range(numGroups)])
            pointsY = np.concatenate([points[i][1].flatten() for i in range(numGroups)])
        else:
            pointsX = points[0]
            pointsY = points[1]
        
    containsAllPoints=np.all([combiPoly.contains(PointS(pointsX[i],pointsY[i])) for i in range(len(pointsX))])
    if not containsAllPoints:
        plt.figure()
        plt.scatter(pointsX,pointsY,c='b')
        boundaryXY = np.array(combiPoly.boundary.coords)
        plt.plot(boundaryXY[:,0],boundaryXY[:,1],c='r')
        raise AttributeError('The provided boundary does not contain all points')
    # Add extra points to ensure that area is finite
    extraPoints = np.array([[np.mean(pointsX),np.max(pointsY)+50],[np.mean(pointsX),np.min(pointsY)-50],\
                             [np.min(pointsX)-50,np.mean(pointsY)],[np.max(pointsX)+50,np.mean(pointsY)],\
                             [np.min(pointsX)-50,np.max(pointsY)+50],[np.min(pointsX)-50,np.min(pointsY)-50],\
                             [np.max(pointsX)+50,np.max(pointsY)+50],[np.max(pointsX)+50,np.min(pointsY)-50]])

    AllPoints = np.array([np.concatenate([pointsX,extraPoints[:,0]]),np.concatenate([pointsY,extraPoints[:,1]])])
    

    vor = Voronoi(AllPoints.T)
    regions = np.array([reg for reg in vor.regions])
    boolval = np.array([len(x)>2 and not -1 in x for x in regions]) # Check if region has at least 3 points and is not connected to infinity (-1))
        
    PolyPoints = np.array([vor.vertices[reg,:] for reg in regions[boolval]])
    polygons = np.array([PolygonS(X) for X in PolyPoints])

    insidePolygonsBool = np.array([combiPoly.contains(P) for P in polygons])

    edgePolygonsBool = np.logical_not(insidePolygonsBool)
    
    intersectionPolygon = []
    for poly in polygons[edgePolygonsBool]:
        inter = poly.intersection(combiPoly)
        if not isinstance(inter,PolygonS): # Not a simple polygon
            inter = inter[np.argmax([x.area for x in inter])] # Return the polygon with biggest area inside boundary
        intersectionPolygon.append(inter)
    
    Polygons = np.concatenate([polygons[np.logical_not(edgePolygonsBool)],intersectionPolygon])
    
    
    if plot or len(pointsX)!=len(Polygons): # pragma: no cover
        plt.figure()
        insiders = np.logical_not(edgePolygonsBool)
        
        [plt.plot(np.array(inter.boundary.coords)[:,0],np.array(inter.boundary.coords)[:,1],c='r') for inter in polygons[insiders]]
        [plt.plot(np.array(inter.boundary.coords)[:,0],np.array(inter.boundary.coords)[:,1],c='g') for inter in intersectionPolygon]
        [plt.plot(np.array(bound.boundary.coords)[:,0],np.array(bound.boundary.coords)[:,1],'-.',c='r') for bound in BoundPoly]
        plt.scatter(extraPoints[:,0],extraPoints[:,1])

        from scipy.spatial import voronoi_plot_2d
        voronoi_plot_2d(vor)
    if not len(pointsX)==len(Polygons):
        raise AttributeError('The number of points given({}) is not the same as the number of polygons created({}). This can be due to many reasons, mainly:\n - Points overlap exactly\n - Points coinsides with the calculated edge\n - ??'.format(len(pointsX),len(Polygons)))


    return Polygons,np.array([np.array(P.boundary.coords[:-1]) for P in Polygons])
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
        raise AttributeError('Provided array has dimension(s) {} of size <= 1'.format(np.arange(xshape)[xshape<=1]))
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
        NormCount =    np.histogramdd(np.array(pos).T,bins=HistBins,weights=np.ones_like(data).flatten())[0].astype(int)
        returndata.append(Normalization)
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
    try: 
        t1 = first.dtype
    except:
        t1 = type(first)
    try:
        t2 = second.dtype
    except:
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

#________________________________________________TESTS_____________________________________________

def test_DataSet_Creation():

    try:
        dataset = DataSet(OtherSetting=10.0,YetAnotherWrongSetting=20.0)
        assert False
    except AttributeError:
        assert True
    dataset = DataSet(Author='Jakob Lass')
    if(dataset.settings['Author']!='Jakob Lass'):
        assert False


def test_Dataset_Initialization():

    emptyDataset = DataSet()
    del emptyDataset
    MJOLNIR.Data.DataFile.assertFile('Data/camea2018n000136.nxs')
    dataset = DataSet(dataFiles=['Data/camea2018n000136.hdf'],convertedFiles='Data/camea2018n000137.nxs',calibrationfiles=[])
    
    assert(dataset.dataFiles[0].name=='camea2018n000136.hdf')
    assert(dataset.convertedFiles[0].name=='camea2018n000137.nxs')
    assert(dataset.normalizationfiles == [])
    Str = str(dataset)

                                                                                                                 
def test_DataSet_Error():
    
    try:
        ds = DataSet(normalizationfiles=[10,11])
        assert False
    except:
        assert True

    ds = DataSet()
    
    try: # No data files
        ds.convertDataFile()
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.binData3D(0.1,0.1,0.1)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.cut1D([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.plotCut1D([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.cutQE([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.plotCutQE([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.cutPowder(np.linspace(0,4,5))
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.plotCutPowder(np.linspace(0,4,5))
        assert False
    except AttributeError:
        assert True
           
    

    try: # Wrong data file type
        ds.dataFiles = 100
        assert False
    except AttributeError:
        assert True


    try: # Can't overwrite settings
        ds.settings={}
        assert False
    except NotImplementedError:
        assert True

    try:# Wrong data file type
        ds.convertedFiles = 10
        assert False
    except AttributeError:
        assert True


    ds.dataFiles = 'Data/camea2018n000136.hdf'

def test_DataSet_Pythonic():
    dataFiles = ['Data/camea2018n000136.hdf','Data/camea2018n000137.hdf']
    dataset = DataSet(dataFiles=dataFiles)
    assert(len(dataset)==2)
    for df in dataset:
        print(df)
    initShape = dataset.I.shape
    names = [dataset[i].name for i in range(len(dataset))]
    names.reverse()
    for i,df in enumerate(list(reversed(dataset))):
        names[i]==df.name

    dataset.append(dataFiles)
    assert(len(dataset)==4)
    secondShape = dataset.Monitor.shape
    assert(np.all(secondShape!=initShape))
    del dataset[3]
    del dataset[2]
    try:
        dataset[10]
        assert False
    except IndexError:
        assert True
    
    try:
        del dataset[10]
        assert False
    except IndexError:
        assert True

    try:
        dataset.append('NoFile')
    except:
        assert True
    
    dataset.append(MJOLNIR.Data.DataFile.DataFile(dataFiles[0]))
    assert(len(dataset)==3)
    assert(dataset.I.shape!=secondShape)



def test_DataSet_Equality():
    D1 = DataSet(dataFiles='Data/camea2018n000136.hdf')#,convertedFiles=['TestData/VanNormalization.nxs'])
    assert(D1==D1)

def test_DataSet_SaveLoad():
    
    D1 = DataSet(dataFiles='Data/camea2018n000136.hdf')#,convertedFiles = 'TestData/VanNormalization.nxs')

    temp = 'temporary.bin'

    D1.save(temp)
    D2 = load(temp)
    os.remove(temp)
    assert(D1==D2) 

def test_DataSet_str():
    D1 = DataSet(dataFiles='Data/camea2018n000136.hdf')#,normalizationfiles = 'TestData/VanNormalization.hdf')
    string = str(D1)
    print(string)


def test_DataSet_Convert_Data():
    dataFiles = 'Data/camea2018n000136.hdf'
    dataset = DataSet(dataFiles=dataFiles)
    

    try:
        dataset.convertDataFile(dataFiles=dataFiles,binning=100)
        assert False
    except AttributeError: # Cant find normalization table
        assert True

    try:
        dataset.convertDataFile(dataFiles='FileDoesNotExist',binning=1)
        assert False
    except AttributeError: # FileDoesNotExist
        assert True

    try:
        os.remove('Data/camea2018n000136.nxs')
    except:
        pass
    dataset.convertDataFile(dataFiles=dataFiles,binning=8,saveLocation='Data/',saveFile=True)
    convertedFile = dataset.convertedFiles[0]
    
    otherFile = MJOLNIR.Data.DataFile.DataFile(dataFiles.replace('.hdf','.nxs'))
    assert(convertedFile==otherFile)
    os.remove('Data/camea2018n000136.nxs')
    


def test_DataSet_3DMesh():
    
    x = np.linspace(0,1,2)
    y = np.linspace(0,1,5)
    z = np.linspace(1,2,5)

    X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
    XT1,YT1,ZT1 = calculateGrid3D(X,Y,Z)

    assert(XT1.shape==(3,6,6))
    assert(np.all(XT1[:,0,0]==np.array([-0.5,0.5,1.5])))
    assert(np.all(YT1[0,:,0]==np.array([-0.125,0.125,0.375,0.625,0.875,1.125])))
    assert(np.all(YT1[0,:,0]==ZT1[0,0,:]-1.0))



def test_DataSet_BinData():
    I = np.random.randint(0,100,(10,20,30))
    Norm = np.random.rand(10,20,30)
    Posx = np.linspace(0,1,10)
    Posy = np.linspace(0,1,20)
    Posz = np.linspace(1,2,30)
    PosX,PosY,PosZ = np.meshgrid(Posx,Posy,Posz,indexing='ij')



    pos = [PosX.flatten(),PosY.flatten(),PosZ.flatten()]
    Data,bins = binData3D(0.5,0.25,0.25,pos,I,norm=Norm)

    ReBinnedI = Data[0]
    RebinnedNorm = Data[1]
    RebinnedNormCount = Data[2]


    assert(ReBinnedI.shape==(3,5,5))
    #assert(np.all(bins[0].shape=(4,6,6)))
    assert(RebinnedNorm.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.dtype==int)
    assert(RebinnedNorm.dtype==Norm.dtype)
    assert(ReBinnedI.dtype==I.dtype)


def test_DataSet_full_test():
    import MJOLNIR.Data.Viewer3D
    
    import matplotlib.pyplot as plt
    import os
    plt.ioff()
    DataFile = ['Data/camea2018n000136.hdf']

    dataset = DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='Data/')

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    
    warnings.simplefilter('ignore')
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')
    viewer = MJOLNIR.Data.Viewer3D.Viewer3D(Intensity,bins)
    viewer = dataset.View3D(0.08,0.08,0.25)
    
    if pythonVersion == 3: # Only possible in python 3
        viewer.ax.set_xticks_number(5)
        viewer.ax.set_yticks_number(10)

        viewer.ax.set_xticks_base(0.5)
        viewer.ax.set_yticks_base(0.5)

    viewer.setProjection(0)
    viewer.setPlane(4)
    del viewer 
    viewer = dataset.View3D(0.08,0.08,0.25,rlu=False)
    os.remove('Data/camea2018n000137.nxs')
    del viewer
    plt.close('all')

def test_DataSet_Visualization():
    import warnings
    from MJOLNIR.Data import Viewer3D,DataFile
    DataFiles = ['Data/camea2018n000136.hdf']

    dataset = DataSet(dataFiles=DataFiles)
    dataset.convertDataFile(saveLocation='Data')

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    Data,bins = dataset.binData3D(0.08,0.08,0.25,dataFiles = [MJOLNIR.Data.DataFile.DataFile('Data/camea2018n000136.nxs')])
    
    plt.ioff()
    warnings.simplefilter('ignore')
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')
    viewer = Viewer3D.Viewer3D(Intensity,bins)
    viewer.caxis = (0,100)
    try:
        viewer.caxis = 'Wrong type'
        assert False
    except AttributeError:
        assert True
    
    try:
        viewer.caxis = [0,1,2,3,4] # Too long input
        assert False
    except AttributeError:
        assert True
    
    try:
        viewer.setAxis(20) # Must bee 0,1, or 2
        assert False
    except AttributeError:
        assert True

    plt.plot()
    plt.close('all')
    os.remove('Data/camea2018n000136.nxs')

def test_DataSet_binEdges():
    X = np.random.rand(100)*3 # array between 0 and 3 -ish
    X.sort()
    tolerance = 0.01
    Bins = _tools.binEdges(X,tolerance=tolerance)

    assert(Bins[0]==X[0]-0.1*tolerance)
    assert(np.isclose(Bins[-1],X[-1],atol=5) or Bins[-1]>X[-1])
    assert(len(Bins)<=3.0/tolerance)
    assert(np.all(np.diff(Bins[:-1])>tolerance*0.99))

def test_DataSet_1Dcut():
    q1 =  np.array([1.23,-1.51])
    q2 =  np.array([1.54, -1.25])
    width = 0.1

    plt.ioff()
    convertFiles = ['Data/camea2018n000136.hdf','Data/camea2018n000137.hdf']
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile(saveFile=False)

    ax,D,P,binCenter,binDistance = ds.plotCut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,fmt='.',ticks=5,tickRound=2)
    D2,P2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5)
    assert(np.all([np.all(D[i]==D2[i]) for i in range(len(D))]))
    assert(np.all([np.all(P[i]==P2[i]) for i in range(len(P))]))

    [intensity,MonitorCount,Normalization,normcounts],bins = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,extend=False)
    assert(np.all(np.logical_and(bins[0][:,0]>=q1[0]-0.1,bins[0][:,0]<=q2[0]+0.1))) 
    # x-values should be between 1.1 and 2.0 correpsonding to q points given (add some extra space due to way bins are created (binEdges))

    #q3 = np.array([1.1,1.1])
    #q4 = np.array([2.0,2.0])
    [intensity,MonitorCount,Normalization,normcounts],bins = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=0.0,Emax=1.5,extend=False)
    assert(np.all(bins[0][:,0]>=q1[0]-0.1))
    assert(np.all(bins[0][:,0]<=q2[0]+0.1))
    assert(np.all(bins[0][:,1]>=q1[1]-0.1))
    assert(np.all(bins[0][:,1]<=q2[1]+0.1))
    # x and y-values should be between 1.1 and 2.0 correpsonding to q points given (add some extra space due to way bins are created (binEdges))

    Q1 = np.array([1,0,0])
    Q2 = np.array([0.5,1,0])

    ax,D,P,binCenter,binDistance = ds.plotCut1D(Q1,Q2,width,rlu=True,minPixel=0.01,Emin=2.0,Emax=2.5,fmt='.')
    D2,P2 = ds.cut1D(Q1,Q2,width,rlu=True,minPixel=0.01,Emin=2.0,Emax=2.5)
    assert(np.all([np.all(D[i]==D2[i]) for i in range(len(D))]))
    assert(np.all(np.array([np.all(np.isclose(P[i],P2[i])) for i in range(len(P))]).flatten()))

    q1,q2 = ds.convertToQxQy([Q1,Q2])
    D1,P1 = ds.cut1D(Q1,Q2,width,rlu=True,minPixel=0.01,Emin=2.0,Emax=2.5)
    D2,P2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5)

    BinPos,OrthoPos,E = P1
    BinPos = np.concatenate([ds.convertToQxQy(BinPos[:,:3]),BinPos[:,-1].reshape(-1,1)],axis=1)
    OrthoPos = ds.convertToQxQy(OrthoPos)
    P1 = [BinPos,OrthoPos,E]

    assert(np.all(np.array([np.all(np.isclose(D[i],D2[i])) for i in range(len(D))]).flatten()))
    assert(np.all(np.array([np.all(np.isclose(P1[i],P2[i])) for i in range(len(P))]).flatten()))
    

def test_DataSet_1DcutE():
    q =  np.array([1.23,-1.25]).reshape(2,1)
    width = 0.1
    Emin = 1.5
    Emax = 2.5
    plt.ioff()
    convertFiles = ['Data/camea2018n000137.hdf']
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    Datset._getData()
    I,qx,qy,energy,Norm,Monitor = Datset.I.flatten(),Datset.qx.flatten(),Datset.qy.flatten(),Datset.energy.flatten(),Datset.Norm.flatten(),Datset.Monitor.flatten()

    [intensity,MonitorCount,Normalization,normcounts],[bins] = cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=Emin,E2=Emax,q=q,width=width,minPixel=0.01)
    Q = Datset.convertToHKL(q.reshape(2))
    
    [intensity,MonitorCount,Normalization,normcounts],[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=Q,width=width,minPixel=0.01)
    assert(np.min(bins)>=Emin-0.01) # Check that bins do not include data outside of cut
    assert(np.max(bins)<=Emax+0.01)
    assert(len(bins)==len(intensity)+1)# Bins denotes edges and must then be 1 more than intensity

    assert(intensity.shape==MonitorCount.shape) # Check that all matrices are cut equally
    assert(intensity.shape==Normalization.shape)
    assert(intensity.shape==normcounts.shape)

    [intensity,MonitorCount,Normalization,normcounts],[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=q,width=width,minPixel=0.01,rlu=False)
    
    Data,[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    
    assert(np.all(np.isclose(np.diff(bins),0.01)))
    assert(bins.min()>=Emin)
    assert(bins.max()<=Emax)

    try: # no points inside energy interval
        cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=500,E2=700,q=q,width=width,minPixel=0.01)
        assert False
    except AttributeError:
        assert True

    try: # no points inside q
        cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=5,E2=7,q=np.array([20.0,0]).reshape(2,1),width=width,minPixel=0.01)
        assert False
    except AttributeError:
        assert True

def test_DataSet_2Dcut():
    q1 =  np.array([1.23,-1.25])
    q2 =  np.array([1.54, -1.51])
    width = 0.1
    minPixel=0.02
    EnergyBins = np.linspace(2,3,4)
    plt.ioff()
    convertFiles = ['Data/camea2018n000137.hdf']
    
    Datset = DataSet(dataFiles = convertFiles)
    try:
        os.remove('Data/camea2018n000137.nxs')
    except:
        pass
    Datset.convertDataFile(saveFile=False)
    ax,Data,pos,cpos,distance = Datset.plotCutQE(q1,q2,width,minPixel,EnergyBins,rlu=False)# Remove to improve test coverage ,vmin=0.0 , vmax= 5e-06)
    Data2,pos2,cpos2,distance2 = Datset.cutQE(q1,q2,width,minPixel,EnergyBins,rlu=False)
    
    for i in range(len(Data)):
        for j in range(len(Data[i])):
            assert(np.all(Data[i][j]==Data2[i][j]))

    for i in range(len(pos)):
        for j in range(len(pos[i])):
            assert(np.all(pos[i][j]==pos2[i][j]))
    
    for i in range(len(cpos)):
        for j in range(len(cpos[i])):
            assert(np.all(cpos2[i][j]==cpos[i][j]))
        
    for i in range(len(distance)):
        for j in range(len(distance[i])):
            assert(np.all(distance2[i][j]==distance[i][j]))

    Q1 = np.array([1,0,0])
    Q2 = np.array([0.5,1,0])

    q1,q2 = Datset.convertToQxQy([Q1,Q2])

    Data1,pos1,cpos1,distance1 = Datset.cutQE(Q1,Q2,width,minPixel,EnergyBins,rlu=True)
    Data2,pos2,cpos2,distance2 = Datset.cutQE(q1,q2,width,minPixel,EnergyBins,rlu=False)
    for i in range(len(Data)):
        for j in range(len(Data[i])):
            assert(np.all(Data1[i][j]==Data2[i][j]))

    for i in range(len(pos)):
        for j in range(len(pos[i])):
            assert(np.all(pos1[i][j]==pos2[i][j]))
    
    for i in range(len(cpos)):
        for j in range(len(cpos[i])):
            assert(np.all(cpos2[i][j]==cpos1[i][j]))
        
    for i in range(len(distance)):
        for j in range(len(distance[i])):
            assert(np.all(distance2[i][j]==distance1[i][j]))

def test_DataSet_cutPowder():
    Tolerance = 0.01

    plt.ioff()
    convertFiles = ['Data/camea2018n000136.hdf']
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    mask = np.ones_like(Datset.I.data)

    Datset.mask = mask
    Datset.mask = np.logical_not(mask)
    

    eBins = _tools.binEdges(Datset.energy,0.25)

    ax,D,q = Datset.plotCutPowder(eBins,Tolerance)# Remove to improve test ,vmin=0,vmax=1e-6)
    D2,q2 = Datset.cutPowder(eBins,Tolerance)

    for i in range(len(D)):
        for j in range(len(D[i])):
            print(D[i][j],D2[i][j])
            assert(np.all(D[i][j]==D2[i][j]))

    for i in range(len(q)):
        for j in range(len(q[i])):
            assert(np.all(q[i][j]==q2[i][j]))

def test_DataSet_createRLUAxes():
    plt.ioff()
    fig = plt.figure()
    convertFiles = ['Data/camea2018n000136.hdf']
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile()

    ax = ds.createRLUAxes()
    ax = ds.createRLUAxes(nbinsx=5)
    ax = ds.createRLUAxes(nbinsy=5)
    ax = ds.createRLUAxes(basex=0.5,figure=fig)
    ax = ds.createRLUAxes(basey=0.5)

    if pythonVersion == 3: # Only possible in python 3
        ax.set_xticks_number(5)
        ax.set_yticks_number(8)

        ax.set_xticks_base(0.2)
        ax.set_yticks_base(0.5)

        ax.set_xticks_number(5)
        ax.set_yticks_number(8)

    V1,V2,V3 = [2,0,0],[-2,3,0],[2,-3,0]
    ax.set_axis(V1,V2)
    ax.set_axis(V1,V2,V3)

    plt.close('all')


def test_DataSet_createQEAxes():
    plt.ioff()
    convertFiles = ['Data/camea2018n000136.hdf']
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile()

    ax = ds.createQEAxes(projectionVector1=ds.sample.projectionVector1,projectionVector2=ds.sample.projectionVector2)

    try:
        ax = ds.createQEAxes(axis=2) # Axis only allowed to be 0 or 1
    except AttributeError:
        assert True
    
    try:
        ax = ds.createQEAxes(projectionVector1=[1,0,0],projectionVector2=[1,2,3,4,5]) # Wrong shape of vector
    except AttributeError:
        assert True
    plt.close('all')



def test_DataSet_plotQPlane():
    plt.ioff()
    convertFiles = ['Data/camea2018n000137.hdf']#'TestData/ManuallyChangedData/A3.hdf']
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()

    EmptyDS = DataSet()
    try:
        Datset.plotQPlane() # No Bins, Emin or Emax
        assert False
    except AttributeError:
        assert True
    try:
        Datset.plotQPlane(EBins=[10]) # Length of bins is 1
        assert False
    except AttributeError:
        assert True
    
    try:
        Datset.plotQPlane(EMin=20,EMax=10) # EMin>EMax
        assert False
    except AttributeError:
        assert True
    
    try:
        EmptyDS.plotQPlane(EMin=2,EMax=3) # Empty DataSet
        assert False
    except AttributeError:
        assert True


    EMin = np.min(Datset.energy)
    EMax = EMin+0.5
    Data,[Qx,Qy],ax1 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=True,log=False,rlu=True)
    Data,[Qx,Qy],ax2 = Datset.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=True,rlu=True)
    fig,AX = plt.subplots()
    Data,[Qx,Qy],ax3 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,ax=AX,colorbar=True,vmin=0,vmax=1e-6,zorder=10)
    
    ax1.set_clim(-20,-15)
    ax2.set_clim(0,1e-6)
    Data,[Qx,Qy],ax3 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05)
    
    cmap = plt.cm.coolwarm

    Dataset = DataSet(dataFiles=convertFiles)
    for d in Dataset.dataFiles:
        d.A3Off +=90 # rotate data to fall into problem of arctan2
    Data,[Qx,Qy],ax2 = Datset.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=True,rlu=True,cmap=cmap)
    QxShape = np.array(Qx[0]).shape
    QyShape = np.array(Qy[0]).shape
    assert(QxShape==QyShape)
    assert(np.all(np.array(Data[0][0]).shape == np.array(QxShape)-np.array([1,1])))
    try:
        Datset.plotQPlane(EMin,EMax,binning='notABinningMethod')
        assert False
    except:
        assert True

    # 3D
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.colors import ListedColormap
    cmap = plt.cm.coolwarm
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)

    fig = plt.figure(figsize=(10,11))
    ax = fig.add_subplot(111, projection='3d')

    Energies = np.concatenate(Datset.energy,axis=0)
    E = np.arange(Energies.min()+0.35,Energies.max(),0.35)


    [I,Monitor,Norm,NormCount],[xBins,yBins],ax = \
    Datset.plotQPlane(EBins=E,ax = ax,xBinTolerance=0.03,yBinTolerance=0.03,
            binning='polar',vmin=7.5e-7,vmax=7e-6,antialiased=True,cmap=cmap,rlu=True,extend='max')
    plt.close('all')

@pytest.mark.unit
def test_DataSet_plotA3A4(quick):
    plt.ioff()

    File1 = 'Data/camea2018n000136.hdf'
    File2 = 'Data/camea2018n000137.hdf'

    DS = DataSet(dataFiles=[File1,File2])
    DS.convertDataFile()

    F1 = DS.convertedFiles[0]
    F2 = DS.convertedFiles[1]

    files = [F1,F2]
    axes = [plt.figure().gca(),plt.figure().gca()]
    try:
        plotA3A4(files,planes=[],ax=axes) # 64 planes and only 2 axes
        assert False
    except AttributeError:
        assert True

    try:
        plotA3A4(files,planes=None,ax=[]) # 64 planes and only 2 axes
        assert False
    except AttributeError:
        assert True 

    try:
        plotA3A4(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
        assert False
    except AttributeError:
        assert True
    
    try:
        ei = F1.Ei
        F2.Ei = F1.Ei*10
        plotA3A4(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
        assert False
    except AttributeError:
        F2.Ei = ei
        assert True

    try:
        plotA3A4(files,planes=[10,[22]],ax=axes,singleFigure=True) # 2 axes and singleFigure true
        assert False
    except AttributeError:
        assert True
    if not quick==True:
        print('___________')
        plotA3A4(files,planes=[10,[22,23]],ax=axes) # Plot plane 10 and 22+23 in the provided axes
        print('___________')
        DS.plotA3A4(planes=[19,[22,25]]) # Plot planes in new axes
        print('___________')
        DS.plotA3A4([F1,F1],planes=[19,[22,25]]) # Plot planes in new axes
        print('___________')
        patches,energies=DS.plotA3A4([F1],planes=[10,25],returnPatches=True)
        print('___________')
        assert(len(patches)==2)
        assert(len(energies)==2)
    plt.close('all')

@pytest.mark.unit
def test_DataSet_plotQPatches(quick):
    assert True
    #     plt.ioff()

#     File1 = 'TestData/T0Phonon10meV.nxs'
#     File2 = 'TestData/T0Phonon10meV93_5A4.nxs'

#     DS = DataSet(convertedFiles=[File1,File2])

#     F1 = DS.convertedFiles[0]
#     F2 = DS.convertedFiles[1]

#     files = [F1,F2]
#     axes = [plt.figure().gca(),plt.figure().gca()]
#     try:
#         plotQPatches(files,planes=[],ax=axes) # 64 planes and only 2 axes
#         assert False
#     except AttributeError:
#         assert True
        
#     try:
#         plotQPatches(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
#         assert False
#     except AttributeError:
#         assert True

#     try:
#         plotQPatches(files,planes=[10,[22]],ax=axes,singleFigure=True) # 2 axes and singleFigure true
#         assert False
#     except AttributeError:
#         assert True

#     if not quick==True:
#         plotQPatches(files,planes=[10,[22,23]],ax=axes) # Plot plane 10 and 22+23 in the provided axes
#         DS.plotQPatches(planes=[19,[22,25]],A4Extend=0.5,A3Extend=1) # Plot planes in new axes
#         DS.plotQPatches(dataFiles=[files[0],files[0]],planes=[19,[22,25]],A4Extend=0.5,A3Extend=1) # Plot planes in new axes and only one file
#     plt.close('all')
    

def test_DataSet_fmt():
    assert('$1.00 \\times 10^{1}$' == fmt(10,'Unused'))
    assert('$1.00 \\times 10^{-10}$' == fmt(1e-10,'Unused'))
    assert('$2.55 \\times 10^{-2}$' == fmt(0.0255,'Unused'))
    assert('$2.56 \\times 10^{-2}$' == fmt(0.02556,'Unused'))
    

def test_DataSet_figureRowColumns():
    assert(np.all(np.array([3,4])==np.array(figureRowColumns(10)))) # 10 -> 3,4
    assert(np.all(np.array([3,3])==np.array(figureRowColumns(9)))) # 9 -> 3,3
    assert(np.all(np.array([1,1])==np.array(figureRowColumns(1)))) # 1 -> 1,1
    try:
        figureRowColumns(0) # No figures
        assert False
    except AttributeError:
        assert True

    assert(np.all(np.array([8,8])==np.array(figureRowColumns(63)))) # 63 -> 8,8
    

def test_DataSet_centeroidnp():
    pos = np.array([[0,0],[1,0],[0,1],[1,1]],dtype=float)
    assert(np.all(np.isclose(np.array([0.5,0.5]),centeroidnp(pos))))

    pos2 = np.array([[1.2,2.2],[7.5,1.0],[11.0,0.0],[4.0,-1.0],[2.0,2.0]],dtype=float)
    assert(np.all(np.isclose(np.array([5.14,0.84]),centeroidnp(pos2))))
    
def test_DataSet_compareNones():
    assert(compareNones(np.array([None]),np.array([None]),0.1))
    assert(not compareNones(np.array([None]),np.array([0.5]),0.1))
    assert(not compareNones(np.array([0.5]),np.array([None]),0.1))
    assert(compareNones(np.array([0.4]),np.array([0.5]),0.2))
    assert(not compareNones(np.array([0.4]),np.array([0.5]),0.001))

    assert(not np.all(compareNones(np.array([0.4,10.2,10.0]),np.array([0.5]),0.001)))
    assert(np.all(compareNones(np.array([0.4,10.2,10.0]),np.array([0.4,10.2,10.0]),0.001)))


def test_DataSet_cutQELine():
    QPoints = np.array([[0.3,-1],[0.7,-1.4],[1.6,-0.9],[0.3,-0.9]],dtype=float)
    QPointsHKL=np.array([[1.0,0.0,0.0],
                        [0.5,1.5,0.0],
                        [1.7,-0.1,0.0],
                        [1.0,1.0,0.0]])


    EnergyBins = np.linspace(1.7,2.7,5)
    minPixel = 0.001
    width=0.1
    DataFile = ['Data/camea2018n000137.hdf']

    dataset = DataSet(convertedFiles=DataFile)
    dataset.convertDataFile(saveFile=False)
    
    try: # No Q-points
        dataset.cutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu=True)
        assert False
    except AttributeError:
        assert True
    try: # Wrong RLU-input
        dataset.cutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu='42') # Wrong RLU-input
        assert False
    except AttributeError:
        assert True

    DataList,BinList,centerPosition,binDistance=dataset.cutQELine(QPoints,EnergyBins,width=width,minPixel=minPixel,rlu=False)
    DataList,BinList,centerPosition,binDistance=dataset.cutQELine(QPointsHKL,EnergyBins,width=width,minPixel=minPixel,rlu=True)
    for i in range(len(DataList)): # Check each segment
        print("i="+str(i))
        assert(DataList[i].shape==(4,len(EnergyBins)-1)) # Must be of size (4,len(E)-1)
        for j in range(DataList[i].shape[1]): # Loop through energies
            print("j="+str(j))   
            assert(np.all([DataList[i][k,j].shape == DataList[i][-1,j].shape for k in range(len(DataList[i])-1)])) # All list are of same size
            assert(centerPosition[i,j].shape == (DataList[i][0,j].shape[0],4))
            assert(centerPosition[i][j].shape[0] == binDistance[i][j].shape[0])
            
    assert(BinList.shape==(len(QPoints)-1,len(EnergyBins)-1,3))

def test_DataSet_plotCutQELine():
    Points = np.array([[0.7140393034102988,-0.4959224853328328],
                        [1.128363301356428,-1.6520150761601147],
                        [1.9002545852012716,-0.9393552598967219],
                        [1.0432282332853056,-0.12375569239528339]],dtype=float)
    QPoints = np.zeros((Points.shape[0],3))
    QPoints[:,:2]=Points
    EnergyBins = np.linspace(1.7,2.7,11)
    minPixel = 0.001
    width=0.1
    
    DataFile = ['Data/camea2018n000136.hdf','Data/camea2018n000137.hdf']
    dataset = DataSet(convertedFiles=DataFile)
    dataset.convertDataFile(saveFile=False)
    
    try: # No Q-points
        dataset.plotCutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu=False)
        assert False
    except AttributeError:
        assert True

    try: # No points in E range
        dataset.plotCutQELine(QPoints,EnergyBins+100,width=width,minPixel=minPixel,rlu=True,vmin=0.0,vmax=1.5e-6,ticks=10)
        assert False
    except AttributeError:
        assert True

    try: # No wrong dim of QPonts
        dataset.plotCutQELine(QPoints,EnergyBins,width=width,minPixel=minPixel,rlu=False)
        assert False
    except AttributeError:
        assert True

    try: # No wrong dim of QPonts
        dataset.plotCutQELine(QPoints[:,:2],EnergyBins,width=width,minPixel=minPixel,rlu=True)
        assert False
    except AttributeError:
        assert True


    fig = plt.figure()
    ax = fig.gca()

    ax,DataList,BinListTotal,centerPositionTotal,binDistanceTotal = dataset.plotCutQELine(
        QPoints[:,:2],EnergyBins,width=width,minPixel=minPixel,rlu=False,ax=ax,vmin=0.0,vmax=1.5e-6,log=True,seperatorWidth=3)


    HKLPoints = np.array([[1.0,0.0,0.0],
                        [0.5,1.5,0.0],
                        [1.7,-0.1,0.0],
                        [1.0,1.0,0.0]])



    ax,DataList,BinListTotal,centerPositionTotal,binDistanceTotal = dataset.plotCutQELine(
        HKLPoints,EnergyBins,width=width,minPixel=minPixel,rlu=True,plotSeperator = False,ticks=1,tickRound=1,colorbar=True,log=True)


    # 3D
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.colors import ListedColormap
    cmap = plt.cm.coolwarm
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)

    fig = plt.figure(figsize=(10,11))
    ax = fig.add_subplot(111, projection='3d')

    Energies = np.concatenate(dataset.energy,axis=0)
    E = np.arange(Energies.min()+0.35,Energies.max(),0.35)
    

    ax,DataList,BinListTotal,centerPositionTotal,binDistanceTotal = \
    dataset.plotCutQELine(QPoints=HKLPoints,EnergyBins=E,ax = ax,width=0.05,minPixel=0.01,
            vmin=7.5e-7,vmax=7e-6,cmap=cmap,rlu=True)

    plt.close('all')


def test_DataSet_extractDetectorData():
    DataFile = ['Data/camea2018n000136.hdf','Data/camea2018n000137.hdf']#['TestData/ManuallyChangedData/A3.nxs','TestData/ManuallyChangedData/A3.nxs']
    dataset = DataSet(DataFile)

    binning = 1
    dataset.convertDataFile(binning=binning)

    
    try:
        dataset.extractDetectorData(A4=10000.0) # A4 outside of detector
        assert False
    except AttributeError:
        assert True

    try:
        dataset.extractDetectorData(Ef=10000.0) # Ef outside of detector
        assert False
    except AttributeError:
        assert True
    

    Efs = dataset.convertedFiles[0].instrumentCalibrationEf[:,1].reshape(104,8*binning)
    AnalyserSelection = 5
    Ef = np.mean(Efs[:,AnalyserSelection])

    A4s = dataset.convertedFiles[0].instrumentCalibrationA4.reshape(104,8*binning)
    DetectorSelection = 19
    A4 = np.mean(A4s[DetectorSelection])-dataset.convertedFiles[0].A4Off


    DatBoth = dataset.extractData(A4=A4,Ef=Ef)
    DatBothId = dataset.extractData(A4=A4,EfId=AnalyserSelection)
    DatOne = dataset.extractData(A4=A4)
    DatOne2= dataset.extractData(Ef=Ef)
    DatAll = dataset.extractData()
    DatAllRaw = dataset.extractData(raw=True)


    # Independent of number of files:
    assert(len(DatAllRaw)==3) # Check that 3 lists are returned
    assert(len(DatAllRaw[0])==len(DatAllRaw[1]) and len(DatAllRaw[0])==len(DatAllRaw[2])) # Check that 3 list have same number of files

    assert(np.all(DatBothId[0]==DatBoth[0])) # Check that ID and value gives the same.

    # The shape of raw is the same as non-raw
    assert(len(DatAllRaw[0])==len(DatAll)) # Have same number of files

    for i in range(len(dataset.convertedFiles)):
        assert(DatAllRaw[0][i].shape==DatAllRaw[1][i].shape and DatAllRaw[0][i].shape==DatAllRaw[2][i].shape) # Check that 3 list have same shape
        assert(DatAllRaw[0][i].shape==DatAll[i].shape) 
        

def test_DataSet_OxfordList():
    l = ['Apples','Pears']
    S = OxfordList(l)
    assert(S=='Apples and Pears')

    l.append('Oranges')
    S = OxfordList(l)
    assert(S=='Apples, Pears, and Oranges')

    assert(OxfordList([]) is None)
    assert(OxfordList(['Apples'])=='Apples')
