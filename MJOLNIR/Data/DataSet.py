import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
#import scipy


import numpy as np
import pickle as pickle
import h5py as hdf

from scipy.ndimage import filters
import scipy.optimize
from scipy.spatial import Voronoi,ConvexHull,KDTree

from shapely.geometry import Polygon as PolygonS
from shapely.geometry import Point as PointS

import matplotlib.pyplot as plt
from mpl_toolkits.axisartist.grid_helper_curvelinear import \
    GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker

import datetime
import warnings
from MJOLNIR.Data import DataFile

import time

dataLocation = 'entry/data/intensity'#'entry/Detectors/Detectors'
EiLocation = 'entry/data/incident_energy' # 'entry/Ei'
monLocation = 'entry/control/data'#'entry/Monitor'

class DataSet(object):
    def __init__(self, dataFiles=None, normalizationfiles=None, calibrationfiles=None, convertedFiles=None, **kwargs):
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


        if dataFiles is not None:
            self.dataFiles = dataFiles

        if normalizationfiles is not None:
            self.normalizationfiles = normalizationfiles
        
        if convertedFiles is not None:
            self.convertedFiles = convertedFiles
            self._getData()

        if calibrationfiles is not None:
            self.calibrationfiles = calibrationfiles


        self._settings = {}
            
        
        
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
            self._dataFiles = isListOfDataFiles(dataFiles)
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
            self._convertedFiles = isListOfDataFiles(convertedFiles)
        except Exception as e:
            raise(e)


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
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')    

    def load(self,filename):
        """Method to load an object from a pickled file."""
        try:                                # Opening the given file with an error catch
            fileObject = open(filename, 'rb')
        except IOError as e:                        # Catch all IO-errors
            print("Error in opening file:\n{}".format(e))
        else:
            tmp_dict = pickle.load(fileObject)
            
            fileObject.close()
            # TODO: Make checks that the object loaded is of correct format?
            self=tmp_dict

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

    def convertDataFile(self,dataFiles=None,binning=8,saveLocation=None):
        """Conversion method for converting scan file(s) to hkl file. Converts the given h5 file into NXsqom format and saves in a file with same name, but of type .nxs.
        Copies all of the old data file into the new to ensure complete reduncency. Determins the binning wanted from the file name of normalization file.

        Kwargs:

            - dataFiles (DataFile, string or list of): File path(s), file must be of hdf format (default self.dataFiles).

            - binning (int): Binning to be used when converting files (default 8).

            - saveLocation (string): File path to save location of data file(s) (defaults to same as raw file).

        Raises:

            - IOError

            - AttributeError
            
        """


        if dataFiles is None:
            if len(self.dataFiles)==0:
                raise AttributeError('No data files file provided either through input of in the DataSet object.')
            dataFiles = self.dataFiles

        dataFiles = isListOfDataFiles(dataFiles)
        #if not isinstance(dataFiles,list):
        #    dataFiles=[dataFiles]
        for datafile in dataFiles:
            if not os.path.isfile(datafile.fileLocation):
                raise AttributeError('Provided file does not exist, '+str(datafile))
            file = hdf.File(datafile.fileLocation,mode='r+')             
            instrument = getInstrument(file)
           
            if instrument.name.split('/')[-1] == 'CAMEA':
                EPrDetector = 8 
            elif instrument.name.split('/')[-1] in ['MULTIFLEXX','FLATCONE']:
                EPrDetector = 1
            
            
            
            normalization = np.array(file.get('entry/calibration/{}_pixels'.format(binning)))
            if(normalization.shape==()):
                raise AttributeError('Binning not found in data file ({})'.format(binning))

            
            
            Data = np.array(instrument.get('detector/data'))

            detectors = Data.shape[1]
            
            liveDetectors = np.array(instrument.get('detector/online'),dtype=bool)
            if liveDetectors is None:
                liveDetectors = np.ones((detectors),dtype=bool)

            if instrument.name.split('/')[-1] in ['MULTIFLEXX','FLATCONE']:
                Data.shape = (Data.shape[0],Data.shape[1],-1)

            A4Zero = file.get('entry/zeros/A4')
            
            if A4Zero is None:
                A4Zero=0.0
            else:
                A4Zero = np.deg2rad(np.array(A4Zero))

            
            A3Zero = file.get('entry/zeros/A3')
            if A3Zero is None:
                A3Zero=0.0
            else:
                A3Zero = np.deg2rad(np.array(A3Zero))

            A4 = np.deg2rad(np.array(normalization[:,9]))+A4Zero
            A4=A4.reshape(detectors,binning*EPrDetector,order='C')
            Ef = np.array(normalization[:,4])
            Ef=Ef.reshape(detectors,binning*EPrDetector,order='C')

            PixelEdge = normalization[:,7:9].reshape(detectors,EPrDetector,binning,2).astype(int)

            factorsqrtEK = 0.694692
            
            A4File = np.array(instrument.get('detector/polar_angle'))
            A4File = A4File.reshape((-1))
            A4Shape = A4.shape
            A4Mean = A4.reshape((1,A4Shape[0],A4Shape[1]))-np.deg2rad(A4File).reshape((A4File.shape[0],1,1))

            DataMean=np.zeros((Data.shape[0],Data.shape[1],EPrDetector*binning),dtype=int)
            for i in range(detectors): # for each detector
                for j in range(EPrDetector):
                    for k in range(binning):
                        DataMean[:,i,j*binning+k] = np.sum(Data[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)

            EfMean = normalization[:,4].reshape(A4.shape[0],EPrDetector*binning)
            Normalization = (normalization[:,3]*np.sqrt(2*np.pi)*normalization[:,5]).reshape(A4.shape[0],EPrDetector*binning)

            kf = factorsqrtEK*np.sqrt(EfMean)
            Ei = np.array(instrument.get('monochromator/energy'))
            
            ki = factorsqrtEK*np.sqrt(Ei)
            
            A3 = np.deg2rad(np.array(file.get('/entry/sample/rotation_angle/')))+A3Zero
            
            #### Qx = ki-kf*cos(A4), Qy = -kf*sin(A4)

            Qx = ki.reshape((ki.shape[0],1,1,1))-(kf.reshape((1,1,EfMean.shape[0],EfMean.shape[1]))*np.cos(A4Mean)).reshape((1,A4Mean.shape[0],A4Mean.shape[1],A4Mean.shape[2]))
            Qy = np.zeros((ki.shape[0],1,1,1))-kf.reshape((1,1,EfMean.shape[0],EfMean.shape[1]))*np.sin(A4Mean.reshape((1,A4Mean.shape[0],A4Mean.shape[1],A4Mean.shape[2])))
            
            QX = Qx.reshape((1,Qx.shape[0],Qx.shape[1],Qx.shape[2],Qx.shape[3]))*np.cos(A3.reshape((A3.shape[0],1,1,1,1)))-Qy.reshape((1,Qy.shape[0],Qy.shape[1],Qy.shape[2],Qy.shape[3]))*np.sin(A3.reshape((A3.shape[0],1,1,1,1)))
            QY = Qx.reshape((1,Qx.shape[0],Qx.shape[1],Qx.shape[2],Qx.shape[3]))*np.sin(A3.reshape((A3.shape[0],1,1,1,1)))+Qy.reshape((1,Qy.shape[0],Qy.shape[1],Qy.shape[2],Qy.shape[3]))*np.cos(A3.reshape((A3.shape[0],1,1,1,1)))

            EnergyShape = (1,len(Ei),1,EfMean.shape[0],EfMean.shape[1])
            DeltaE = (Ei.reshape((Ei.shape[0],1,1))-EfMean.reshape((1,EfMean.shape[0],EfMean.shape[1]))).reshape(EnergyShape
            Intensity = DataMean.reshape((QX.shape[0],QX.shape[1],QX.shape[2],QX.shape[3],QX.shape[4]))
        
            DeltaE=DeltaE.repeat(QX.shape[0],axis=0)
            DeltaE=DeltaE.repeat(QX.shape[2],axis=2)
            
            Monitor = np.array(file.get('/entry/control/data'),dtype=int)
            Monitor.shape = (len(Monitor),1,1)
            Monitor = np.repeat(Monitor,Intensity.shape[3],axis=1)
            Monitor = np.repeat(Monitor,Intensity.shape[4],axis=2)
            Monitor.shape = Intensity.shape

            Normalization.shape = (1,1,1,Normalization.shape[0],Normalization.shape[1])
            Normalization = np.repeat(Normalization,Intensity.shape[0],axis=0)
            Normalization = np.repeat(Normalization,Intensity.shape[1],axis=1)
            Normalization = np.repeat(Normalization,Intensity.shape[2],axis=2)
            ## TODO: Don't let all things vary at the same time!!

            # Filter out offline detectors
            #for attribute in [Intensity,Monitor,Normalization,QX,QY,DeltaE]:
            #    attribute = attribute[:,:,:,liveDetectors,:]
            
            if not saveLocation is None:
                if saveLocation[-1]!='/':
                    saveLocation+='/'
                saveloc = saveLocation+datafile.fileLocation.replace('.h5','.nxs').split('/')[-1]
            else:
                saveloc = datafile.fileLocation.replace('.h5','.nxs')
            saveNXsqom(datafile,file,saveloc,Intensity,Monitor,QX,QY,DeltaE,binning,Normalization)
            
            file.close()
            convFil = DataFile.DataFile(saveloc)
            self.convertedFiles.append(convFil)
        self._getData()
            
    def _getData(self): # Internal method to populate I,qx,qy,energy,Norm and Monitor
        self.I,self.qx,self.qy,self.energy,self.Norm,self.Monitor,self.a3,self.a4,self.instrumentCalibration,self.Ei = DataFile.extractData(self.convertedFiles)
        


    def binData3D(self,dx,dy,dz,dataFiles=None):
        """Bin a converted data file into voxels with sizes dx*dy*dz. Wrapper for the binData3D functionality.

        Args:

            - dx (float): step sizes along the x direction (required).

            - dy (float): step sizes along the y direction (required).

            - dz (float): step sizes along the z direction (required).

        Kwargs:

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
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor = DataFile.extractData()
            
        

        pos=[qx.flatten(),qy.flatten(),energy.flatten()]

        returnData,bins = binData3D(dx,dy,dz,pos,I,norm=Norm,mon=Monitor)

        return returnData,bins

    def cut1D(self,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False,dataFiles=None):
        """Wrapper for 1D cut through constant energy plane from q1 to q2 function returning binned intensity, monitor, normalization and normcount. The full width of the line is width while height is given by Emin and Emax. 
        the minimum step sizes is given by minPixel.
        
        .. note::
            Can only perform cuts for a constant energy plane of definable width.
        
        Args:
            
            - q1 (2D array): Start position of cut in format (qx,qy).
            
            - q2 (2D array): End position of cut in format (qx,qy).
            
            - width (float): Full width of cut in q-plane.
            
            - minPixel (float): Minimal size of binning aling the cutting direction. Points will be binned if they are closer than minPixel.
            
            - Emin (float): Minimal energy to include in cut.
            
            - Emax (float): Maximal energy to include in cut
            
        Kwargs:
            
            - plotCoverage (bool): If True, generates plot of all points in the cutting plane and adds bounding box of cut (default False).

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.
        
        
        Returns:
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
        """
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor = DataFile.extractData()
            
        positions = [qx,qy,energy]
        positions = [qx,qy,energy]
        
        return cut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False)

    def plotCut1D(self,q1,q2,width,minPixel,Emin,Emax,ax=None,plotCoverage=False,dataFiles=None,**kwargs):  
        """Plotting wrapper for the cut1D method. Generates a 1D plot with bins at positions corresponding to the distance from the start point. 
        Adds the 3D position on the x axis with ticks.
        
        .. note::
            Can only perform cuts for a constant energy plane of definable width.
        
        Kwargs:
            
            - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
            
            - kwargs: All other keywords will be passed on to the ax.errorbar method.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - binCenter (3D array): Array containing the position of the bin centers of size (n,3)
            
            - binDistance (array): Distance from centre of bins to start position.
        """
        
        
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor = DataFile.extractData()
            
        positions = [qx,qy,energy]

        return plotCut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,ax,plotCoverage,**kwargs)


    def cutQE(self,q1,q2,width,minPixel,EnergyBins,dataFiles=None):
        """Wrapper for cut data into maps of q and intensity between two q points and given energies. This is performed by doing consecutive constant energy planes.

        Args:

            - q1 (2D array): Start position of cut in format (qx,qy).
            
            - q2 (2D array): End position of cut in format (qx,qy).
            
            - width (float): Full width of cut in q-plane.
            
            - minPixel (float): Minimal size of binning aling the cutting direction. Points will be binned if they are closer than minPixel.

            - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

        Kwargs:

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.
    

        Returns:
            
            - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
            - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (n * 3D arrays): n instances of center positions for the bins.

            - binDistance (n arrays): n isntances of arrays holding the distance in q to q1.

        """
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor = DataFile.extractData()
            
        positions = [qx,qy,energy]
        
        return cutQE(positions,I,Norm,Monitor,q1,q2,width,minPixel,EnergyBins)

    def plotCutQE(self,q1,q2,width,minPixel,EnergyBins,ax=None,dataFiles=None,**kwargs): 
        """Plotting wrapper for the cutQE method. Generates a 2D intensity map with the data cut by cutQE. 
    
        .. note::
            Positions shown in tool tip reflect the closes bin center and are thus limited to the area where data is present.
        
        Args:

            - q1 (2D array): Start position of cut in format (qx,qy).
            
            - q2 (2D array): End position of cut in format (qx,qy).
            
            - width (float): Full width of cut in q-plane.
            
            - minPixel (float): Minimal size of binning aling the cutting direction. Points will be binned if they are closer than minPixel.

            - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

        Kwargs:
            
            - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.
        
            - kwargs: All other keywords will be passed on to the ax.errorbar method.
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
            - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
            - center position (n * 3D arrays): n instances of center positions for the bins.

            - binDistance (n arrays): n isntances of arrays holding the distance in q to q1.
        """
        
        
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor = DataFile.extractData()
            
        positions = [qx,qy,energy]

        return plotCutQE(positions,I,Norm,Monitor,q1,q2,width,minPixel,EnergyBins,ax = None,**kwargs)


    def cutPowder(self,EBinEdges,qMinBin=0.01,dataFiles=None):
        """Cut data powder map with intensity as function of the length of q and energy. 

        Args:
            
            - EBinEdges (list): Bin edges between which the cuts are performed.

        Kwargs:

            - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.

            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.


        Returns:
            
            - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
            
            - qbins (n arrays): n arrays holding the bin edges along the lenght of q

        """
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor = DataFile.extractData()
            
        positions = [qx,qy,energy]

        return cutPowder(positions,I,Norm,Monitor,EBinEdges,qMinBin)

    def plotCutPowder(self,EBinEdges,qMinBin=0.01,ax=None,dataFiles=None,**kwargs):
        """Plotting wrapper for the cutPowder method. Generates a 2D plot of powder map with intensity as function of the length of q and energy.  
        
        .. note::
            Can only perform cuts for a constant energy plane of definable width.
        
        Args:

            - EBinEdges (list): Bin edges between which the cuts are performed.

        Kwargs:
            
            - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.
            
            - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
            
            - dataFiles (list): List of dataFiles to cut (default None). If none, the ones in the object will be used.

            - kwargs: All other keywords will be passed on to the ax.pcolormesh method.
        
        Returns:
            
            - ax (matplotlib axis): Matplotlib axis into which the plot was put.
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).

        """
        if dataFiles is None:
            if len(self.convertedFiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                self._getData()#datafiles = self.convertedFiles
                I = self.I
                qx = self.qx
                qy = self.qy
                energy = self.energy
                Norm = self.Norm
                Monitor = self.Monitor

        else: 
            dataFiles = isListOfDataFiles(dataFiles)
            I,qx,qy,energy,Norm,Monitor,Ei = DataFile.extractData()
            
        positions = [qx,qy,energy]

        return plotCutPowder(positions,I,Norm,Monitor,EBinEdges,qMinBin,ax,**kwargs)

    def createRLUAxes(self): # pragma: no cover
        """Wrapper for the createRLUAxes method.

        Returns:

            - ax (Matplotlib axes): Created reciprocal lattice axes.

        .. note::
           Uses sample from the first converted data file. However, this should be taken care of by the comparison of datafiles to ensure same sample and settings.

        """
        return createRLUAxes(self)


    def plotQPlane(self,EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=False,ax=None,RLUPlot=True,**kwargs):
        """Wrapper for plotting tool to show binned intensities in the Q plane between provided energies.
        
        Args:
            
            - EMin (float): Lower energy limit.
            
            - EMax (float): Upper energy limit.
            
        Kwargs:
            
            - binning (str): Binning scheme, either 'xy' or 'polar' (default 'xy').
            
            - xBinTolerance (float): bin sizes along x direction (default 0.05). If enlargen is true, this is the minimum bin size.

            - yBinTolerance (float): bin sizes along y direction (default 0.05). If enlargen is true, this is the minimum bin size.
            
            - enlargen (bool): If the bin sizes should be adaptive (default False). If set true, bin tolereces are used as minimum bin sizes.

            - log (bool): Plot intensities as the logarithm (defautl False).
            
            - ax (matplotlib axes): Axes in which the data is plotted (default None). If None, the function creates a new axes object.

            - RLUPlot (bool): If true and axis is None, a new reciprocal lattice axis is created and used for plotting (default True).
            
            - other: Other key word arguments are passed to the pcolormesh plotting algorithm.
            
        Returns:
            
            - ax (matplotlib axes)
            
        .. note::
            The axes object gets a new method denoted 'set_clim' taking two parameters (VMin and VMax) used to change axes coloring.
            
            
        """
        I = self.I
        qx = self.qx
        qy = self.qy
        energy = self.energy
        Norm = self.Norm
        Monitor = self.Monitor

        if ax is None and RLUPlot is True:
            ax = self.createRLUAxes()
        
        pos = [qx,qy,energy]
        return plotQPlane(I,Monitor,Norm,pos,EMin,EMax,binning=binning,xBinTolerance=xBinTolerance,yBinTolerance=yBinTolerance,enlargen=enlargen,log=log,ax=ax,**kwargs)

    def plotA3A4(self,files=None,ax=None,dimension='2D',planes=[],binningDecimals=3,singleFigure=False,plotTesselation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
        """Plot data files together with pixels created around each point in A3-A4 space. Data is binned in the specified planes through their A3 and A4 values. 
        This can result in distordet binning when binning across large energy regions. Data is plotted using the pixels calulated for average plane value, i.e. 
        binning 7,8,9,10, and 11 patches for plane 9 are used for plotting.

        Kwargs:
            - files (DataFiles): single file or list of files to be binned together (Default self.convertedFiles)

            - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)

            - dimension ('2D' or '3D'): Plot data in 2 or 3 dimensions (default '2D')

            - planes (list (of lists)): Planes to be plotted and binned (default [])

            - binningDecimals (int): Number of decimal places Q positions are rounded before binning (default 3)

            - singleFigure (bool): If true, all planes are plotted in same figure (default False)

            - plotTesselation (bool): Plot tesselation of points (default False)

            - Ei_err (float): Tolerence of E_i for which the values are equal (default = 0.05)

            - temperature_err (float): Tolerence of temperature for which the values are equal (default = 0.2)
            
            - magneticField_err (float): Tolerence of magnetic field for which the values are equal (default = 0.2)
            
            - electricField_err (float): Tolerence of electric field for which the values are equal (default = 0.2)

        Returns:
            
            - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.

        Raises:

            - NotImplimentedError

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
        if files is None:
            files = self.convertedFiles
        
        return plotA3A4(files,ax=ax,dimension=dimension,planes=planes,binningDecimals=binningDecimals,
        singleFigure=singleFigure,plotTesselation=plotTesselation,Ei_err=Ei_err,temperature_err=temperature_err,\
        magneticField_err=magneticField_err,electricField_err=electricField_err)

    def plotQPatches(self,files=None,ax=None,dimension='2D',planes=[],binningDecimals=3,A4Extend=0.2,A3Extend=0.5,singleFigure=False,plotTesselation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
        """Plot data files together with pixels created around each point in Q space. 

        .. warning::
           This method plots all measurement points unless they are literaly on top of each other and is thus really slow! Binning 8 planes for two files takes approximately
           3.5 minutes. Alternatively use binning, i.e. plotQPlane.

        Kwargs:

            - files (DataFiles): single file or list of files to be binned together (Default self.convertedFiles)

            - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)

            - dimension ('2D' or '3D'): Plot data in 2 or 3 dimensions (default '2D')

            - planes (list (of lists)): Planes to be plotted and binned (default [])

            - binningDecimals (int): Number of decimal places Q positions are rounded before binning (default 3)

            - A4Extend (float): Angle value with which the boundary is extended away from points in A4 direction (default 0.2)
            
            - A3Extend (float): Angle value with which the boundary is extended away from points in A3 direction (default 0.5)

            - singleFigure (bool): If true, all planes are plotted in same figure (default False)

            - plotTesselation (bool): Plot tesselation of points (default False)

            - Ei_err (float): Tolerence of E_i for which the values are equal (default = 0.05)

            - temperature_err (float): Tolerence of temperature for which the values are equal (default = 0.2)
            
            - magneticField_err (float): Tolerence of magnetic field for which the values are equal (default = 0.2)
            
            - electricField_err (float): Tolerence of electric field for which the values are equal (default = 0.2)

        Returns:
            
            - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.

        Raises:


            - AttributeError

        Examples: REDO!

        The following example will combine the two files and plot all of the available planes in different figures.

        >>> DS = DataSet.DataSet(convertedFiles=[--.nxs,---.nxs])
        >>> plt.figure()
        >>> ax = plt.gca()
        >>>
        >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax)

        If only a subset of planes or different planes are to be combined the following will achieve this:

        >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax,planes=[0,1,2,3,[4,5,6],[8,9]])

        Here planes 0 through 3 are plotted separately while 4,5, and 6 as well as 8 and 9 are binned.

        .. note::
            Binning planes from different analysers might result in nonsensible binnings.

        """
        if files is None:
            files = self.convertedFiles
        
        return plotQPatches(files,ax=ax,dimension=dimension,planes=planes,binningDecimals=binningDecimals,A4Extend=A4Extend,A3Extend=A3Extend,singleFigure=singleFigure,\
        plotTesselation=plotTesselation,Ei_err=Ei_err,temperature_err=temperature_err,\
        magneticField_err=magneticField_err,electricField_err=electricField_err)


def cut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False):
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
        
        - minPixel (float): Minimal size of binning aling the cutting direction. Points will be binned if they are closer than minPixel.
        
        - Emin (float): Minimal energy to include in cut.
        
        - Emax (float): Maximal energy to include in cut
        
    Kwargs:
        
        - plotCoverage (bool): If True, generates plot of all points in the cutting plane and adds bounding box of cut (default False).
    
    Returns:
        
        - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
        
        - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
    """
    dirvec = np.array(q2,dtype=float)-np.array(q1,dtype=float)
    dirvec/=np.linalg.norm(dirvec)
    orthovec=np.array([dirvec[1],-dirvec[0]])
    
    ProjectMatrix = np.array([dirvec,orthovec])
    
    insideEnergy = np.logical_and(positions[2]<=Emax,positions[2]>=Emin)
    if(len(insideEnergy)==0):
        raise AttributeError('No points are within the provided energy limits.')

    positions2D = np.array([positions[0][insideEnergy], positions[1][insideEnergy]])
    propos = np.dot(ProjectMatrix,positions2D-q1.reshape(2,1))
    
    orthobins = [-width/2.0,width/2.0]
    insideWidth = np.logical_and(propos[1]<orthobins[1],propos[1]>orthobins[0])
    
    lenbins = np.array(binEdges(propos[0][insideWidth],minPixel))
    orthopos = np.outer(orthobins,orthovec)
    binpositions = np.outer(lenbins,dirvec)+q1
    
    if len(lenbins)==0:
        return [np.array(np.array([])),np.array([]),np.array([]),np.array([])],[np.array([]),orthopos,[Emin,Emax]]
    
    normcounts = np.histogramdd(propos.T,bins=[lenbins,orthobins],weights=np.ones((propos.shape[1])).flatten())[0]
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
         for i in range(len(binpositions)):
             plt.plot([binpositions[i][0]+orthopos[0][0],binpositions[i][0]+orthopos[1][0]],[binpositions[i][1]+orthopos[0][1],binpositions[i][1]+orthopos[1][1]],c='k',linewidth=0.5)
         plt.scatter(positions2D[0][insideWidth],positions2D[1][insideWidth],s=0.5)
         ax = plt.gca()
         ax.set_aspect('equal', 'datalim')
         ax.set_xlabel('Qx [1/A]')
         ax.set_ylabel('Qy [1/A]')
    return [intensity,MonitorCount,Normalization,normcounts],[binpositionsTotal,orthopos,np.array([Emin,Emax])]



        
def binEdges(values,tolerance):
    """Generate binning of values array with minimum bin size of tolerance. Binning starts at values[0]-tolerance/2.0 and ends at values[-1]+tolerance/2.0.
    
    Args:
        
        - values (array): 1D array to be binned.
        
        - tolerance (float): Minimum length of bin sizes.
        
    Returns:
        
        - bins (array)
    
    """
    values_array = np.array(values).ravel()
    unique_values = np.asarray(list(set(values_array)))
    unique_values.sort()
    if len(unique_values)==0:
        return []
    #    bin_edges = [unique_values[0] - tolerance / 2]
    #    for i in range(len(unique_values) - 1):
    #        if unique_values[i+1] - unique_values[i] > tolerance:
    #            bin_edges.append((unique_values[i] + unique_values[i+1]) / 2)
    #        else:
    #            pass
    #    
    #    bin_edges.append(unique_values[-1] + tolerance / 2)
    bin_edges = [unique_values[0] - tolerance / 2.0]
    add = 1
    current = 0
    while current<len(unique_values) - 2:
        add=1
        while unique_values[current+add] - unique_values[current] < tolerance:
            if current+add < len(unique_values) - 2:
                add+=1
            else:
                current=len(unique_values)-add-1
                break
        bin_edges.append((unique_values[current] + unique_values[current+add]) / 2)
        current+=add+1
    bin_edges.append(unique_values[-1] + tolerance / 2)
    return np.array(bin_edges)

def plotCut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,ax=None,plotCoverage=False,**kwargs):
    """Plotting wrapper for the cut1D method. Generates a 1D plot with bins at positions corresponding to the distance from the start point. 
    Adds the 3D position on the x axis with ticks.
    
    .. note::
        Can only perform cuts for a constant energy plane of definable width.
    
    Kwargs:
        
        - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
        
        
        - kwargs: All other keywords will be passed on to the ax.errorbar method.
    
    Returns:
        
        - ax (matplotlib axis): Matplotlib axis into which the plot was put.
        
        - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
        
        - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
        - binCenter (3D array): Array containing the position of the bin centers of size (n,3)
        
        - binDistance (array): Distance from centre of bins to start position.
    """
    
    
    q1 = np.array(q1)
    q2 = np.array(q2)
    
    D,P = cut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,plotCoverage)
    INT = np.divide(D[0]*D[3],D[1]*D[2])
    INT_err = np.divide(np.sqrt(D[0])*D[3],D[1]*D[2])
    
    
    binCenter = 0.5*(P[0][:-1]+P[0][1:])
    num=len(binCenter)
    
    xvalues = np.round(np.linspace(0,num-1,5)).astype(int)
    my_xticks=[]
    for i in xvalues:
        my_xticks.append('\n'.join(map(str,[np.round(binCenter[i,0],2),np.round(binCenter[i,1],2),np.round(binCenter[i,2],2)])))
    
    binDistance = np.linalg.norm(binCenter-P[0][0],axis=1)
    
    if ax is None:
        plt.figure()
        ax = plt.gca()
    
    ax.errorbar(binDistance,INT,yerr=INT_err,**kwargs)

    ax.set_xticks(binDistance[xvalues])
    ax.set_xticklabels(my_xticks,fontsize=10, multialignment="center",ha="center")
    ax.set_xlabel('$Q_h/A$\n$Q_k/A$\nE/meV', fontsize=10)
    ax.xaxis.set_label_coords(1.15, -0.025)
    ax.set_ylabel('Int [arb]')
    plt.tight_layout()
    # TODO: Incorporate RLU figure/axis
    def format_coord(x,y,binDistance,binCenter):# pragma: no cover
        index = np.argmin(np.abs(binDistance-x))
        qx,qy,E = binCenter[index]
        return  "qx = {0:.3f}, qy = {1:.3f}, E = {2:.3f}, I = {3:0.4e}".format(qx,qy,E,y)
    
    ax.format_coord = lambda x,y: format_coord(x,y,binDistance,binCenter)

    return ax,D,P,binCenter,binDistance



def cutPowder(positions,I,Norm,Monitor,EBinEdges,qMinBin=0.01):
    """Cut data powder map with intensity as function of the length of q and energy. 

    Args:

        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.

        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - EBinEdges (list): Bin edges between which the cuts are performed.

    Kwargs:

        - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.

    Returns:
        
        - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
        
        - qbins (n arrays): n arrays holding the bin edges along the lenght of q

    """
    qx,qy,energy = positions
    q = np.linalg.norm([qx,qy],axis=0)
    intensity = []
    monitorCount = []
    Normalization = []
    NormCount = []
    qbins = []
    
    for i in range(len(EBinEdges)-1):
        e_inside = np.logical_and(energy>EBinEdges[i],energy<=EBinEdges[i+1])
        q_inside = q[e_inside]
        qbins.append(np.array(binEdges(q_inside,tolerance=qMinBin)))
            
        intensity.append(np.histogram(q_inside,bins=qbins[-1],weights=I[e_inside].flatten())[0].astype(I.dtype))
        monitorCount.append(np.histogram(q_inside,bins=qbins[-1],weights=Monitor[e_inside].flatten())[0].astype(Monitor.dtype))
        Normalization.append(np.histogram(q_inside,bins=qbins[-1],weights=Norm[e_inside].flatten())[0].astype(Norm.dtype))
        NormCount.append(np.histogram(q_inside,bins=qbins[-1],weights=np.ones_like(I[e_inside]).flatten())[0].astype(I.dtype))
    
    return [intensity,monitorCount,Normalization,NormCount],qbins



def plotCutPowder(positions, I,Norm,Monitor,EBinEdges,qMinBin=0.01,ax=None,**kwargs):
    """Plotting wrapper for the cutPowder method. Generates a 2D plot of powder map with intensity as function of the length of q and energy.  
    
    .. note::
       Can only perform cuts for a constant energy plane of definable width.
    
    Args:

        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.

        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - EBinEdges (list): Bin edges between which the cuts are performed.

    Kwargs:
        
        - qMinBin (float): Minimal size of binning along q (default 0.01). Points will be binned if they are closer than qMinBin.
        
        - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
        
        - kwargs: All other keywords will be passed on to the ax.pcolormesh method.
    
    Returns:
        
        - ax (matplotlib axis): Matplotlib axis into which the plot was put.
        
        - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
        
        - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).

    """
    
    [intensity,monitorCount,Normalization,NormCount],qbins = cutPowder(positions,I,Norm,Monitor,EBinEdges,qMinBin)
    Int = [np.divide(intensity[i]*NormCount[i],monitorCount[i]*Normalization[i]) for i in range(len(EBinEdges)-1)]
    
    eMean = 0.5*(EBinEdges[:-1]+EBinEdges[1:])
    
    if ax is None:
        plt.figure()
        ax = plt.gca()
    pmeshs = []
    
    for i in range(len(EBinEdges)-1):
        pmeshs.append(ax.pcolormesh(qbins[i],[EBinEdges[i],EBinEdges[i+1]],Int[i].reshape((len(qbins[i])-1,1)).T,**kwargs))
    
    
    def format_coord(x,y,qBin,eMean,Int):# pragma: no cover
            EIndex = np.argmin(np.abs(y-eMean))
            qIndex = np.argmin(np.abs(x-0.5*(qBin[EIndex][:-1]+qBin[EIndex][1:])))
            Intensity = Int[EIndex][qIndex] 
            return  "|q| = {0:.3f}, E = {1:.3f}, I = {2:0.4e}".format(qBin[EIndex][qIndex],eMean[EIndex],Intensity)
        
    ax.format_coord = lambda x,y: format_coord(x,y,qbins,eMean,Int)
    ax.set_xlabel('|q| [1/A]')
    ax.set_ylabel('E [meV]')
    
    ax.set_clim = lambda VMin,VMax: [pm.set_clim(VMin,VMax) for pm in pmeshs]
    
    if not 'vmin' in kwargs or not 'vmax' in kwargs:
        minVal = np.min(np.concatenate(Int))
        maxVal = np.max(np.concatenate(Int))
        ax.set_clim(minVal,maxVal)
    ax.pmeshs = pmeshs
    return ax,[intensity,monitorCount,Normalization,NormCount],qbins
 




def cutQE(positions,I,Norm,Monitor,q1,q2,width,minPix,EnergyBins):
    """Cut data into maps of q and intensity between two q points and given energies. This is performed by doing consecutive constant energy planes.

    Args:

        - positions (3 arrays): position in Qx, Qy, and E in flattend arrays.
        
        - I (array): Flatten intensity array
        
        - Norm (array): Flatten normalization array
        
        - Monitor (array): Flatten monitor array
        
        - q1 (2D array): Start position of cut in format (qx,qy).
        
        - q2 (2D array): End position of cut in format (qx,qy).
        
        - width (float): Full width of cut in q-plane.
        
        - minPixel (float): Minimal size of binning aling the cutting direction. Points will be binned if they are closer than minPixel.

        - EnergyBins (list): Bin edges between which the 1D constant energy cuts are performed.

    Returns:
        
        - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
        
        - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
        - center position (n * 3D arrays): n instances of center positions for the bins.

        - binDistance (n arrays): n isntances of arrays holding the distance in q to q1.

    """
    intensityArray = []
    monitorArray = []
    normalizationArray = []
    normcountArray = []
    centerPos = []
    returnpositions = []
    binDistance = []
    
    for i in np.arange(len(EnergyBins)-1):
        [intensity,MonitorCount,Normalization,normcounts],position = cut1D(positions,I,Norm,Monitor,q1,q2,width,minPix,EnergyBins[i],EnergyBins[i+1],plotCoverage=False)
        if len(intensity)==0:
            continue
        returnpositions.append(position)
        intensityArray.append(intensity)
        monitorArray.append(MonitorCount)
        normalizationArray.append(Normalization)
        normcountArray.append(normcounts)
        centerPos.append(0.5*(position[0][:-1]+position[0][1:]))
        binDistance.append(np.linalg.norm(centerPos[-1][:,:2]-position[0][0][:2],axis=1))#np.linalg.norm(centerPos[-1][:,:2]-q1,axis=1))
    
    
    
    return [intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance

def plotCutQE(positions,I,Norm,Monitor,q1,q2,width,minPix,EnergyBins,ax = None,**kwargs):
    """Plotting wrapper for the cutQE method. Generates a 2D intensity map with the data cut by cutQE. 
    
    .. note::
        Positions shown in tool tip reflect the closes bin center and are thus limited to the area where data is present.
    

    Kwargs:
        
        - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
        
        
        - kwargs: All other keywords will be passed on to the ax.errorbar method.
    
    Returns:
        
        - ax (matplotlib axis): Matplotlib axis into which the plot was put.
        
        - Data list (n * 4 arrays): n instances of [Intensity, monitor count, normalization and normalization counts].
        
        - Bin list (n * 3 arrays): n instances of bin edge positions in plane of size (m+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
        - center position (n * 3D arrays): n instances of center positions for the bins.

        - binDistance (n arrays): n isntances of arrays holding the distance in q to q1.
    """

    [intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance = cutQE(positions,I,Norm,Monitor,q1,q2,width,minPix,EnergyBins)
    
    if ax is None:
        plt.figure()
        ax = plt.gca()
    
    dirvec = np.array(q2)-np.array(q1)
    leftEdgeBins = np.array([np.dot(returnpositions[i][0][0,:2],dirvec) for i in range(len(returnpositions))])
    rightEdgeBins = np.array([np.dot(returnpositions[i][0][-1,:2],dirvec) for i in range(len(returnpositions))])

    leftEdgeIndex = np.argmin(leftEdgeBins)
    rightEdgeIndex= np.argmax(rightEdgeBins)

    binedges = [np.linalg.norm(returnpositions[i][0][:,:2]-returnpositions[leftEdgeIndex][0][0,:2],axis=1) for i in range(len(returnpositions))]
    
    #print(binedges)

    binenergies = [returnpositions[i][2] for i in range(len(returnpositions))]
    warnings.simplefilter('ignore')
    Int = [np.divide(intensityArray[i]*normcountArray[i],monitorArray[i]*normalizationArray[i]) for i in range(len(intensityArray))]
    warnings.simplefilter('once')
    pmeshs = []
    for i in range(len(Int)):
        pmeshs.append(ax.pcolormesh(binedges[i],binenergies[i],Int[i].T,**kwargs))
    
    ax.set_clim = lambda VMin,VMax: [pm.set_clim(VMin,VMax) for pm in pmeshs]
    
    if not 'vmin' in kwargs or not 'vmax' in kwargs:
        minVal = np.min(np.concatenate(Int))
        maxVal = np.max(np.concatenate(Int))
        ax.set_clim(minVal,maxVal)

    binCenter = centerPos[0]#np.concatenate(centerPos,axis=0)
    #binCenter.sort(axis=0)
    binDistanceAll = binDistance[0]#np.linalg.norm(binCenter[:,:2]-q1,axis=1)
    #print(binDistanceAll)
    
    num=len(binDistanceAll)
    
    # leftPoint = returnpositions[leftEdgeIndex][0][0,:2]
    # rightPoint = returnpositions[rightEdgeIndex][0][-1,:2]
    # diffVector = rightPoint-leftPoint

    xvalues = np.round(np.linspace(0,num-1,5)).astype(int)

    my_xticks=[]
    for i in xvalues:
        my_xticks.append('\n'.join(map(str,[np.round(binCenter[i,0],2),np.round(binCenter[i,1],2)])))
    
    
    def format_coord(x,y,binDistance,centerPos,Int):# pragma: no cover
        Eindex = np.argmin(np.abs([x[0][2] for x in centerPos]-y))
        index = np.argmin(np.abs(binDistance[Eindex]-x))
        qx,qy,E = centerPos[Eindex][index]
        Intensity = Int[Eindex][index][0]
        return  "qx = {0:.3f}, qy = {1:.3f}, E = {2:.3f}, I = {3:.3e}".format(qx,qy,E,Intensity)
    
    ax.format_coord = lambda x,y: format_coord(x,y,binDistance,centerPos,Int)


    ax.set_xticks(binDistanceAll[xvalues])
    ax.set_xticklabels(my_xticks,fontsize=10, multialignment="center",ha="center")
    ax.set_xlabel('$Q_h/A$\n$Q_k/A$', fontsize=10)
    ax.xaxis.set_label_coords(1.15, -0.025)
    ax.set_ylabel('E [meV]')
    plt.tight_layout()
    ax.pmeshs = pmeshs
    return ax,[intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance

def createRLUAxes(Dataset): # pragma: no cover
    """Create a reciprocal lattice plot for a given DataSet object.
    
    Args:
        
        - Dataset (DataSet): DataSet object for which the RLU plot is to be made.
        
    Returns:
        
        - ax (Matplotlib axes): Axes containing the RLU plot.
    
    """
    fileObject = Dataset.convertedFiles[0]
    
    fig = plt.figure(figsize=(7, 4))
    fig.clf()
    grid_helper = GridHelperCurveLinear((fileObject.sample.tr, fileObject.sample.inv_tr))
    
    ax = Subplot(fig, 1, 1, 1, grid_helper=grid_helper)
  
    fig.add_subplot(ax)
    ax.set_aspect(1.)
    ax.grid(True, zorder=0)
    
    ax.format_coord = fileObject.sample.format_coord
    projV1 = fileObject.sample.orientationMatrix[0].astype(int)
    projV2 = fileObject.sample.orientationMatrix[1].astype(int)
    ax.set_xlabel('hkl = [{0:d},{1:d},{2:d}]'.format(projV1[0],projV1[1],projV1[2]))
    ax.set_ylabel('hkl = [{0:d},{1:d},{2:d}]'.format(projV2[0],projV2[1],projV2[2]))
    return ax

def plotQPlane(I,Monitor,Norm,pos,EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=False,ax=None,**kwargs):
    """Plotting tool to show binned intensities in the Q plane between provided energies.
    
    Args:
        
        - I (array): Intensity of data.
        
        - Monitor (array): Monitor of data.
        
        - Norm (array): Nornmalization of data.
        
        - pos (3 array): Position of data in qx, qy, and energy.
        
        - EMin (float): Lower energy limit.
        
        - EMax (float): Upper energy limit.
        
    Kwargs:
        
        - binning (str): Binning scheme, either 'xy' or 'polar' (default 'xy').
        
        - xBinTolerance (float): bin sizes along x direction (default 0.05). If enlargen is true, this is the minimum bin size.

        - yBinTolerance (float): bin sizes along y direction (default 0.05). If enlargen is true, this is the minimum bin size.
        
        - enlargen (bool): If the bin sizes should be adaptive (default False). If set true, bin tolereces are used as minimum bin sizes.

        - log (bool): Plot intensities as the logarithm (defautl False).
        
        - ax (matplotlib axes): Axes in which the data is plotted (default None). If None, the function creates a new axes object.
        
        - other: Other key word arguments are passed to the pcolormesh plotting algorithm.
        
    Returns:
        
        - ax (matplotlib axes)
        
    .. note::
        The axes object gets a new method denoted 'set_clim' taking two parameters (VMin and VMax) used to change axes coloring.
        
        
    """
    qx,qy,energy=pos


    
    if ax is None:
    #        if RLUPlot:
    #            ax = self.createRLUAxes()
    #        else:
        plt.figure()
        ax = plt.gca()
            
    
    
    binnings = ['xy','polar']#,'rlu']
    if not binning in binnings:
        raise AttributeError('The provided binning is not understood, should be {}'.format(', '.join(binnings)))
    if binning == 'polar':
        x = np.arctan2(qy,qx)
        y = np.linalg.norm([qx,qy],axis=0)  
    
    elif binning == 'xy':
        x = qx
        y = qy
        
    #elif binning == 'rlu':
    #    raise NotImplementedError('Currently the RLU binning is not implimented')
    
    
    EBinEdges = [EMin,EMax]

    intensity = []
    monitorCount = []
    Normalization = []
    NormCount = []
    bins = []
    
    
    for i in range(len(EBinEdges)-1):
        e_inside = np.logical_and(energy>EBinEdges[i],energy<=EBinEdges[i+1])
        if enlargen:
            yBins = binEdges(y[e_inside],yBinTolerance)
        else:
            yBins = np.arange(np.min(y[e_inside]),np.max(y[e_inside]),yBinTolerance)
        for j in range(len(yBins)-1):
            ey_inside = np.logical_and(np.logical_and(e_inside,np.logical_and(y>yBins[j],y<yBins[j+1])),(1-np.isnan(Norm)).astype(bool))
            
            x_inside = x[ey_inside]
            #y_inside = y[ey_inside]
            
            if enlargen:
                xbins = binEdges(x_inside,tolerance=xBinTolerance)
            else:
                xbins = np.arange(np.min(x),np.max(x),xBinTolerance)
                
            if len(xbins)==0:
                continue
            bins.append(np.array([xbins,np.array([yBins[j],yBins[j+1]])]))
            
            intensity.append(np.histogram(x_inside,bins=bins[-1][0],weights=I[ey_inside].flatten())[0].astype(I.dtype))
            monitorCount.append(np.histogram(x_inside,bins=bins[-1][0],weights=Monitor[ey_inside].flatten())[0].astype(Monitor.dtype))
            Normalization.append(np.histogram(x_inside,bins=bins[-1][0],weights=Norm[ey_inside].flatten())[0].astype(Norm.dtype))
            NormCount.append(np.histogram(x_inside,bins=bins[-1][0],weights=np.ones_like(I[ey_inside]).flatten())[0].astype(I.dtype))

    warnings.simplefilter('ignore')
    Int = [np.divide(intensity[i]*NormCount[i],monitorCount[i]*Normalization[i]) for i in range(len(intensity))]
    warnings.simplefilter('once')
    
    if binning == 'polar':
        Qx = [np.outer(bins[i][1],np.cos(bins[i][0])).T for i in range(len(intensity))]
        Qy = [np.outer(bins[i][1],np.sin(bins[i][0])).T for i in range(len(intensity))]
    
    elif binning == 'xy':
        Qx = [np.outer(bins[i][0],np.ones_like(bins[i][1])) for i in range(len(intensity))]
        Qy = [np.outer(np.ones_like(bins[i][0]),bins[i][1]) for i in range(len(intensity))]
        
   
    pmeshs = []
    if log:
        Int = [np.log(1e-20+np.array(Int[i])) for i in range(len(Int))]
    for i in range(len(intensity)):
        pmeshs.append(ax.pcolormesh(Qx[i],Qy[i],Int[i].reshape((len(Int[i]),1)),zorder=10,**kwargs))
    ax.set_aspect('equal', 'datalim')
    ax.grid(True, zorder=0)
    ax.set_clim = lambda VMin,VMax: [pm.set_clim(VMin,VMax) for pm in pmeshs]
    ax.pmeshs = pmeshs
    return ax

def plotA3A4(files,ax=None,dimension='2D',planes=[],binningDecimals=3,singleFigure=False,plotTesselation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
    """Plot data files together with pixels created around each point in A3-A4 space. Data is binned in the specified planes through their A3 and A4 values. 
    This can result in distordet binning when binning across large energy regions. Data is plotted using the pixels calulated for average plane value, i.e. 
    binning 7,8,9,10, and 11 patches for plane 9 are used for plotting.

    Args:
        
        - files (DataFiles): single file or list of files to be binned together

    Kwargs:

        - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)

        - dimension ('2D' or '3D'): Plot data in 2 or 3 dimensions (default '2D')

        - planes (list (of lists)): Planes to be plotted and binned (default [])

        - binningDecimals (int): Number of decimal places A3-A4 positions are rounded before binning (default 3)

        - singleFigure (bool): If true, all planes are plotted in same figure (default False)

        - plotTesselation (bool): Plot tesselation of points (default False)

        - Ei_err (float): Tolerence of E_i for which the values are equal (default = 0.05)

        - temperature_err (float): Tolerence of temperature for which the values are equal (default = 0.2)
        
        - magneticField_err (float): Tolerence of magnetic field for which the values are equal (default = 0.2)
        
        - electricField_err (float): Tolerence of electric field for which the values are equal (default = 0.2)

    Returns:
        
        - ax (matplotlib axis or list of): axis (list of) containing figures for plotted planes.

    Raises:

        - NotImplimentedError

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
    if dimension!='2D':
        raise NotImplementedError('Only 2D plotting is currently supported')
    
    if not isinstance(ax, (list,)) and ax is not None:
        ax = np.array([ax])
    
    if not isinstance(planes, (list,)):
        planes = np.array([planes])
        
    if not ax is None:
        if singleFigure and np.array([ax]).size != 1:
            raise AttributeError('Single figure chosen but multiple axes given ({}).'.format(np.array([ax]).size))
        
        elif not singleFigure and len(ax) != len(planes) and not len(planes)==0:
            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,len(planes)))
            
    
    files = np.asarray(files)
    numFiles = len(files)

    
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
    

    A4All = np.array([files[i].A4 for i in range(numFiles)]) # +files[i].A4Off Already taken care of in the file
    A3All = np.array([files[i].A3 for i in range(numFiles)]) # +files[i].A3Off 
    
    Ishape = files[0].I.shape
    IAll = np.array([files[i].I[:,0,0,:,:].reshape((A3All[i].size,Ishape[3],Ishape[4])) for i in range(numFiles)]) # into shape sum(A3),104,64 for CAMEA
    NormAll = np.array([files[i].Norm[:,0,0,:,:].reshape((A3All[i].size,Ishape[3],Ishape[4])) for i in range(numFiles)])
    MonitorAll = np.array([files[i].Monitor[:,0,0,:,:].reshape((A3All[i].size,Ishape[3],Ishape[4])) for i in range(numFiles)])
    
    if not ax is None:
        if not singleFigure and len(ax) != Ishape[4] and len(planes) == 0: # Plot all planes in provided axes
            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,Ishape[4]))

    I = np.concatenate(IAll,axis=0)
    Norm = np.concatenate(NormAll,axis=0)
    Mon = np.concatenate(MonitorAll,axis=0)

    A4InstrAll = np.array([files[i].instrumentCalibration[:,9]-A4All[i] for i in range(numFiles)])
    
    # Find binning (All are equal through testing)
    binning = files[0].binning
    if binning==1:
        if A4InstrAll.shape[1]==155: #MULTIFLEXX
            A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,5,binning))
        else: # FLATCONE
            A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,1,binning))  
    else:
        A4InstrAll = np.reshape(A4InstrAll,(numFiles,-1,8,binning))

        A4InstrAll[:,92]-=0.1 # TODO: Change this!!
    
    ####################################################################### Assume that all energies have same A4
    A4InstrAll = A4InstrAll.reshape(numFiles,A4InstrAll[0].shape[0],-1)[:,:,0]

    # Generate measured points in A3-A4 space
    points = []

    for i in range(numFiles):
        X,Y = [x.flatten() for x in np.meshgrid(A3All[i],A4InstrAll[i],indexing='ij')]
        points.append([X,Y])
    
    PosAll = np.concatenate(points,axis=1)
    unique,uindex,count = np.unique(PosAll,axis=1,return_index=True,return_counts=True)
    
    if np.sum(count>1)>0: # If there is any duplicate points
        #print('Binning!')
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
        index = np.lexsort((PosAll[1], PosAll[0]))

        shape = I.shape[2] #(64 or 8 depending on instrument and binning)
        Isorted = I.reshape(-1,shape)[index,:]
        Normsorted = Norm.reshape(-1,shape)[index,:]
        Monsorted = Mon.reshape(-1,shape)[index,:]
    

    polygons,GoodPolyPoints = voronoiTesselation(points,plot = plotTesselation,Boundary = BoundPoly)

    # Sort centroids (i.e. polygons) like measurement points
    centroids = np.array([centeroidnp(x) for x in GoodPolyPoints]).T

    if isinstance(points,list):
        X = np.concatenate(points,axis=1).T
    else:
        X = points.T
    Y = centroids.T
    def closest_node(node, nodes):
        nodes = np.asarray(nodes)
        deltas = nodes - node
        dist_2 = np.einsum('ij,ij->i', deltas, deltas)
        return np.argmin(dist_2)
    
    A = [closest_node(X,y) for y in Y]
        
        
    _,SortUindex,SortCount = np.unique(A,return_index=True,return_counts=True)
    if np.sum(SortCount>1)!=0:
        raise AttributeError('The number of points connecting the centroids from tesselation and points are not equal...')
    centInd = SortUindex
    
    sortedPolyPoints = GoodPolyPoints[centInd]
    factorsqrtEK = 0.694692
    
    # Calcualte k vectors
    Ei = files[0].Ei
    ki = np.sqrt(Ei)*factorsqrtEK
    kf = np.sqrt(Ei-files[0].energy[0,0,0,:,:].mean(axis=0))*factorsqrtEK
    
    
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
    E = np.mean(files[0].energy,axis=(0,1,2,3))
    
    if len(planes)==0:
        planes = range(len(E))
        
    plots = len(planes)
    if ax is None: # Create needed axes
        if singleFigure: # create only one
            rows,cols = figureRowColumns(plots)
            fig,ax = plt.subplots(nrows=rows, ncols=cols)
            ax = np.array(ax).flatten()
    if singleFigure:
        if ax is None:
            ax = plt.figure().gca()
    else:
        if ax is None:
            ax = [plt.figure().gca() for _ in range(plots)]
            
    counter = 0
    for plane in planes:
        
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
        
        
         # Generate polygons in Qspace
        patches = [Polygon(np.array([QRX[i][:,plotPlane],QRY[i][:,plotPlane]]).T) for i in range(len(QRX))]
        pcollection = PatchCollection(patches)
        
        currentInt = IntensityBin#
        
        pcollection.set_array(currentInt)
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
    return ax


def plotQPatches(files,ax=None,dimension='2D',planes=[],binningDecimals=3,A4Extend=0.2,A3Extend=0.5,singleFigure=False,plotTesselation=False,Ei_err = 0.05,temperature_err=0.2,magneticField_err=0.2,electricField_err=0.2):
    """Plot data files together with pixels created around each point in Q space. 

    .. warning::
        This method plots all measurement points unless they are literaly on top of each other and is thus really slow! Binning 8 planes for two files takes approximately
        3.5 minutes. Alternatively use binning, i.e. plotQPlane.


    Args:
        
        - files (DataFiles): single file or list of files to be binned together

    Kwargs:

        - ax (matplotlib axis): Axis into which the planes are to be plotted (Default None, i.e. new)

        - dimension ('2D' or '3D'): Plot data in 2 or 3 dimensions (default '2D')

        - planes (list (of lists)): Planes to be plotted and binned (default [])

        - binningDecimals (int): Number of decimal places Q positions are rounded before binning (default 3)

        - A4Extend (float): Angle value with which the boundary is extended away from points in A4 direction (default 0.2)
        
        - A3Extend (float): Angle value with which the boundary is extended away from points in A3 direction (default 0.5)

        - singleFigure (bool): If true, all planes are plotted in same figure (default False)

        - plotTesselation (bool): Plot tesselation of points (default False)

        - Ei_err (float): Tolerence of E_i for which the values are equal (default = 0.05)

        - temperature_err (float): Tolerence of temperature for which the values are equal (default = 0.2)
        
        - magneticField_err (float): Tolerence of magnetic field for which the values are equal (default = 0.2)
        
        - electricField_err (float): Tolerence of electric field for which the values are equal (default = 0.2)

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
    >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax)

    If only a subset of planes or different planes are to be combined the following will achieve this:

    >>> DataSet.plotQPatches(DS.convertedFiles,ax=ax,planes=[0,1,2,3,[4,5,6],[8,9]])

    Here planes 0 through 3 are plotted separately while 4,5, and 6 as well as 8 and 9 are binned.

    .. note::
        Binning planes from different analysers might result in nonsensible binnings.

    """
    t1 = time.time()
    if dimension!='2D':
        raise NotImplementedError('Only 2D plotting is currently supported')
    
    if not isinstance(ax, (list,)) and ax is not None:
        ax = np.array([ax])
    
    if not isinstance(planes, (list,)):
        planes = np.array([planes])
        
    if not ax is None:
        if singleFigure and np.array([ax]).size != 1:
            raise AttributeError('Single figure chosen but multiple axes given ({}).'.format(np.array([ax]).size))
        
        elif not singleFigure and len(ax) != len(planes) and not len(planes)==0:
            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,len(planes)))
            
    
    files = np.asarray(files)
    numFiles = len(files)

    
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
    
    Ishape = files[0].I.shape
    if not ax is None:
        if not singleFigure and len(ax) != Ishape[4] and len(planes) == 0: # Plot all planes in provided axes
            raise AttributeError('Number of axes ({}) provided does not match number of planes ({}).'.format(np.array([ax]).size,Ishape[4]))

    
    IAll = np.array([files[i].I[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)]) # into shape sum(A3),104,64 for CAMEA
    NormAll = np.array([files[i].Norm[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
    MonitorAll = np.array([files[i].Monitor[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
  
    I = np.concatenate(IAll,axis=0)
    Norm = np.concatenate(NormAll,axis=0)
    Mon = np.concatenate(MonitorAll,axis=0)
    
    QxAll = np.array([files[i].qx[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
    QyAll = np.array([files[i].qy[:,0,0,:,:].reshape((-1,Ishape[3],Ishape[4])) for i in range(numFiles)])
    Qx = np.concatenate(QxAll,axis=0)
    Qy = np.concatenate(QyAll,axis=0)

    
    if len(planes)==0:
        planes = range(len(I.shape[-1]))
    
    plots = len(planes)
    if ax is None: # Create needed axes
        if singleFigure: # create only one
            rows,cols = figureRowColumns(plots)
            fig,ax = plt.subplots(nrows=rows, ncols=cols)
            ax = np.array(ax).flatten()
    if singleFigure:
        if ax is None:
            ax = plt.figure().gca()
    else:
        if ax is None:
            ax = [plt.figure().gca() for _ in range(plots)]
    t2 = time.time()
    counter = 0
    for plane in planes:
        mp = []
        for i in range(len(files)):
            xx = boundaryQ(files[i],plane,A4Extend=A4Extend,A3Extend=A3Extend)
            polygons = [PolygonS(x.T) for x in xx.transpose(1,0,2)]
            if isinstance(plane,list):
                if len(plane)>1:
                    mplocal = polygons[0]
                    for j in range(len(polygons)-1):
                        mplocal = mplocal.union(polygons[j+1])
                    mp.append(mplocal)
                else:
                    mp.append(polygons[0])
            else:
                mp.append(polygons[0])
        
        
        if len(mp)>1:
            boundary = mp[0]
            for i in range(len(mp)-1):
                boundary = boundary.union(mp[i+1])
            boundary = [boundary]
        else:
            boundary = mp
        
        t3 = time.time()
        if isinstance(plane,list) or isinstance(plane,np.ndarray):         
            IAlive = []
            NormAlive = []
            MonAlive = []
            QxAlive = []
            QyAlive = []
            for i in range(len(plane)):
                alive = np.logical_not(np.isnan(Norm[:,:,plane[i]]))
                IAlive.append(I[alive,plane[i]])
                NormAlive.append(Norm[alive,plane[i]])
                MonAlive.append(Mon[alive,plane[i]])
                
                QxAlive.append(Qx[alive,plane[i]])
                QyAlive.append(Qy[alive,plane[i]])
            IAlive = np.concatenate(IAlive)
            NormAlive = np.concatenate(NormAlive)
            MonAlive = np.concatenate(MonAlive)
            QxAlive = np.concatenate(QxAlive)
            QyAlive = np.concatenate(QyAlive)
        else:
            alive = np.logical_not(np.isnan(Norm[:,:,plane]))
            IAlive = I[alive,plane]
            NormAlive = Norm[alive,plane]
            MonAlive = Mon[alive,plane]
            QxAlive = Qx[alive,plane]
            QyAlive = Qy[alive,plane]
            

        t4 = time.time()
        points = np.array([QxAlive,QyAlive])
        unique,uindex = np.unique(np.round(points,binningDecimals),axis=1,return_index=True)
        t41 = time.time()
        if unique.shape[1]!=points.shape[1]:
            #print('BINNING!')
            mask = np.ones(points.shape[1],dtype=bool)
            mask[uindex] = False
            doublePoints = points[:,mask]
            t42 = time.time()
            kdtree = KDTree(unique.T)
            t425= time.time()
            doubleIndex = kdtree.query(np.round(doublePoints,binningDecimals).T,distance_upper_bound=np.power(10,-binningDecimals*1.0)*1.1)[1]

            t43 = time.time()
            points = unique
           
            Isorted = IAlive[uindex]
            Normsorted = NormAlive[uindex]
            Monsorted = MonAlive[uindex]

            t44 = time.time()

            IAliveDouble = IAlive[mask]
            NormAliveDouble = NormAlive[mask]
            MonAliveDouble = MonAlive[mask]

            t45 = time.time()

            Isorted[doubleIndex]+=IAliveDouble
            Normsorted[doubleIndex]=np.mean([Normsorted[doubleIndex],NormAliveDouble],axis=0)
            Monsorted[doubleIndex]+=MonAliveDouble

            t46 = time.time()
            currentInt = np.divide(Isorted,Normsorted*Monsorted)
            t47 = time.time()
        else:

            pointIndex = np.lexsort((points[1], points[0]))
            currentInt = np.divide(IAlive[pointIndex],NormAlive[pointIndex]*MonAlive[pointIndex])
            
        t5 = time.time()    
        polygons,GoodPolyPoints = voronoiTesselation([points],plot = plotTesselation,Boundary = boundary)
        t6 = time.time()
        centroids = np.array([np.array(x.centroid.coords).reshape(2) for x in polygons]).T
        t7 = time.time()

        X = unique.T
        Y = Y = centroids.T
        

        kdtree = KDTree(X)

        A = kdtree.query(Y,distance_upper_bound=0.2)[1]

        t8 = time.time()
        SortUnique,SortUindex,SortCount = np.unique(A,return_index=True,return_counts=True)
        if np.sum(SortCount>1)!=0:
            raise AttributeError('The number of points tieing the centroids and Q poinst together are not equal, difference is {}...'.format(np.sum(SortCount>1)))
        patchIndex = SortUindex
        t9 = time.time()
        E = files[0].energy
        patches = [Polygon(np.array([np.array(x.boundary.coords)[:,0],np.array(x.boundary.coords)[:,1]]).T) for x in polygons[patchIndex]]
        pcollection = PatchCollection(patches)
        
        try:
            bpoints = np.array(boundary[0].boundary.coords)
        except:
            bpoints = np.concatenate([np.array(x.boundary.coords) for x in boundary[0]])
            
        qxmin = np.min(bpoints[:,0])
        qymin = np.min(bpoints[:,1])
        qxmax = np.max(bpoints[:,0])
        qymax = np.max(bpoints[:,1])
        
        QXlim = np.max(np.abs([qxmin,qxmax]))
        QYlim = np.max(np.abs([qymin,qymax]))
        
        pcollection.set_array(currentInt)
        pcollection.set_edgecolor('face')
        currIntMin = np.max([np.nanmin(currentInt),0.0])
        pcollection.set_clim(currIntMin,np.nanmax(currentInt))
        
        ax[counter].add_collection(pcollection)
        ax[counter].set_xlim(-QXlim,QXlim)
        ax[counter].set_ylim(-QYlim,QYlim)
        ax[counter].axes.grid('on')
        ax[counter].get_figure().colorbar(ax[counter].collections[0], ax=ax[counter],format=ticker.FuncFormatter(fmt))
        
        ax[counter].collections[0].set_clim(currIntMin,np.max(currentInt))
        if not isinstance(plane,list):
            ax[counter].set_title('Energy {0:.3f} meV - plane {1}'.format(np.mean(E[:,:,:,:,plane]),plane))
        else:
            if len(plane) == 1:
                ax[counter].set_title('Energy {0:.3f} meV - plane {1}'.format(np.mean(E[:,:,:,:,plane]),plane))
            else:
                ax[counter].set_title('Energy {0:.3f} meV - planes '.format(np.mean(E[:,:,:,:,plane]))+\
                  ','.join([str(x) for x in plane]))
        counter +=1
        t10 = time.time()
    if False:
        print('Initialization: {}s'.format(t2-t1))
        print('Boundary: {}s'.format(t3-t2))
        print('Alive sorting: {}s'.format(t4-t3))
        print('Binning: {}s'.format(t5-t4))
        print('\tUnique: {}\n\tMask: {}\n\tGenerate KDTree: {}\n\tDouble index: {}\n\tUnique Mat: {}\n\tApply mask: {}\n\tAdd together: {}\n\tCalculate I: {}'.format(t41-t4,t42-t41,t425-t42,t43-t425,t44-t43,t45-t44,t46-t45,t47-t46))
        print('Voronoi: {}s'.format(t6-t5))
        print('Centroid Calc: {}s'.format(t7-t6))
        print('Sorting: {}s'.format(t8-t7))
        print('Unique: {}s'.format(t9-t8))
        print('Plot: {}s'.format(t10-t9))
    return ax



def voronoiTesselation(points,plot=False,Boundary=False):
    """Generate individual pixels around the given datapoints.

    Args:

        - points (list of list of points): Data points to generate pixels in shape [files,XY,N] i.e. [1,2,N] for one file with N points

    Kwargs:

        - plot (bool): If True, method plots pixels created with green as edge bins and red as internal (default False)

        - Boundary (lost of Polygons): List of Shapely polygons constituting the boundaries (Default False)


    """
    numGroups = len(points)
    
    if Boundary==False:
        BoundPoly= [convexHullPoints(points[i][0].flatten(),points[i][1].flatten()) for i in range(numGroups)]
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

    # Sort vertecies for the polygon generation
    #GoodPolyPoints = PolyPoints[np.array([not np.any(x==-1) for x in PolyPoints],dtype=bool)]
    
       
    polygons = np.array([PolygonS(X) for X in PolyPoints])

    insidePolygonsBool = np.array([combiPoly.contains(P) for P in polygons])
    #crosses = np.array([combiPoly.crosses(P) for P in polygons])

    edgePolygonsBool = np.logical_not(insidePolygonsBool)
    
    intersectionPolygon = []#[poly.intersection(combiPoly) for poly in polygons[edgePolygonsBool]]
    for poly in polygons[edgePolygonsBool]:
        #try: # If polygon side coincides with boundary intersection fails and following is needed
        inter = poly.intersection(combiPoly)
        if not isinstance(inter,PolygonS): # Not a simple polygon
            #print(inter)
            #print(type(inter))
            #[x.xy for x in inter]
            inter = inter[np.argmax([x.area for x in inter])] # Return the polygon with biggest area inside boundary
            #inter.boundary.coords
       # except:

        #       diff = poly.difference(combiPoly)
        #    inter = poly.difference(diff)
    
        intersectionPolygon.append(inter)
    
    Polygons = np.concatenate([polygons[np.logical_not(edgePolygonsBool)],intersectionPolygon])
    
    
    if plot or len(pointsX)!=len(Polygons):
        plt.figure()
        insiders = np.logical_not(edgePolygonsBool)
        
        [plt.plot(np.array(inter.boundary.coords)[:,0],np.array(inter.boundary.coords)[:,1],c='r') for inter in polygons[insiders]]
        [plt.plot(np.array(inter.boundary.coords)[:,0],np.array(inter.boundary.coords)[:,1],c='g') for inter in intersectionPolygon]
        [plt.plot(np.array(bound.boundary.coords)[:,0],np.array(bound.boundary.coords)[:,1],'-.',c='r') for bound in BoundPoly]
        plt.scatter(extraPoints[:,0],extraPoints[:,1])

        from scipy.spatial import voronoi_plot_2d
        voronoi_plot_2d(vor)
    if not len(pointsX)==len(Polygons):
        raise AttributeError('The number of points given({}) is not the same as the number of polygons created({}). This can be due to many reasons, mainly:\n - Points overlap exactly\n - Points coinsides with the calulated edge\n - ??'.format(len(pointsX),len(Polygons)))


    return Polygons,np.array([np.array(P.boundary.coords[:-1]) for P in Polygons])


def boundaryQ(file,plane,A4Extend=0.0,A3Extend=0.0):
    """Calculate the boundary of a given scan in Q space
    A4Extend: in degrees
    A3Extend: in degrees
    """
    energy = file.energy[:,0,0,:,:]
    A3 = file.A3+file.A3Off
    
    A4 = file.A4-file.A4Off
    Ei = file.Ei
    IC = file.instrumentCalibration
        
    InstrumentA4 = IC[:,-1].reshape(energy.shape[1],-1)[:,plane]
    
    factorsqrtEK = 0.694692
    InstA4 = (InstrumentA4-A4)*np.pi/180.0 
    A4Min = np.min(InstA4,axis=0)
    A4Max = np.max(InstA4,axis=0)
    
    InstrumentEnergy = IC[:,4].reshape(energy.shape[1],-1)[:,plane]
    
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
    if isinstance(inputFiles,list):
        for file in inputFiles:
            if isinstance(file,DataFile.DataFile):
                returnList.append(file)
            elif isinstance(file,str):
                returnList.append(DataFile.DataFile(file))
    elif isinstance(inputFiles,DataFile.DataFile):
        returnList.append(inputFiles)
    elif isinstance(inputFiles,str):
        returnList.append(DataFile.DataFile(inputFiles))
    else:
        raise AttributeError('File provided is not of type string, list, or DataFile')
    if len(returnList)>1:
        sameSample = [returnList[0].sample==file.sample for file in returnList]
        if not np.all(sameSample):
            raise AttributeError('Files does not have the same sample! Compared to first entry: {}'.format(sameSample))
    return returnList


def saveNXsqom(datafile,fs,savefilename,Intensity,Monitor,QX,QY,DeltaE,binning,Normalization):
    
    fd = hdf.File(savefilename,'w')
    group_path = fs['/entry'].parent.name
    
    group_id = fd.require_group(group_path)
    
    
    fs.copy('/entry', group_id, name="/entry")
    
    definition = fd.create_dataset('entry/definition',(1,),dtype='S70',data=np.string_('NXsqom'))
    definition.attrs['NX_class'] = 'NX_CHAR'
    
    process = fd.create_group('entry/reduction')
    process.attrs['NX_class']=b'NXprocess'
    proc = process.create_group('MJOLNIR_algorithm_convert')
    proc.attrs['NX_class']=b'NXprocess'
    author= proc.create_dataset('author',shape=(1,),dtype='S70',data=np.string_('Jakob Lass'))
    author.attrs['NX_class']=b'NX_CHAR'
    
    date= proc.create_dataset('date',shape=(1,),dtype='S70',data=np.string_(datetime.datetime.now()))
    date.attrs['NX_class']=b'NX_CHAR'
    
    description = proc.create_dataset('description',shape=(1,),dtype='S70',data=np.string_('Conversion from pixel to Qx,Qy,E in reference system of instrument.'))
    description.attrs['NX_class']=b'NX_CHAR'
    
    rawdata = proc.create_dataset('rawdata',shape=(1,),dtype='S200',data=np.string_(datafile.name))
    rawdata.attrs['NX_class']=b'NX_CHAR'
    
    normalizationString = proc.create_dataset('binning',shape=(1,),dtype='int32',data=binning)
    normalizationString.attrs['NX_class']=b'NX_INT'
    
    data = fd.get('entry/data')
    #data['rawdata']=data['data']
    #del data['data']
    
    
    fileLength = Intensity.shape
    #print(Intensity.shape)
    Int = data.create_dataset('intensity',shape=(fileLength),dtype='int32',data=Intensity)
    Int.attrs['NX_class']='NX_INT'
    
    monitor = data.create_dataset('monitor',shape=(fileLength),dtype='int32',data=Monitor)
    monitor.attrs['NX_class']=b'NX_INT'

    normalization = data.create_dataset('normalization',shape=(fileLength),dtype='float32',data=Normalization)
    normalization.attrs['NX_class']=b'NX_FLOAT'
    
    qx = data.create_dataset('qx',shape=(fileLength),dtype='float32',data=QX)
    qx.attrs['NX_class']=b'NX_FLOAT'
    
    qy = data.create_dataset('qy',shape=(fileLength),dtype='float32',data=QY)
    qy.attrs['NX_class']=b'NX_FLOAT'
    
    #qz = data.create_dataset('qz',shape=(fileLength,),dtype='float32',data=np.zeros((fileLength,)))
    #qz.attrs['NX_class']=b'NX_FLOAT'
    
    en = data.create_dataset('en',shape=(fileLength),dtype='float32',data=DeltaE)
    en.attrs['NX_class']=b'NX_FLOAT'

    fd.close()


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

    Now XX is a 21x11x67 array containing all x coordinates of the edges exactly midway bewteen the points. Same goes for YY and ZZ with y and z coordinates respectively.
    """

    xshape = X.shape
    
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

        Rebinned intensity (and if provided Normalization, Monitor, and Normalization Count) and X, Y, and Z bins in 3 3D arrays.


    Example:

    >>> pos = [Qx,Qy,E]
    >>> Data,bins = DataSet.binData3D(0.05,0.05,0.2,pos,I,norm=Norm,mon=Monitor)

    """

    if bins is None:
        bins = calculateBins(dx,dy,dz,pos)
    
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
    
    XX,YY,ZZ = calculateGrid3D(X,Y,Z)
    
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
    location = file.visititems(lambda x,y: getNX_class(x,y,b'NXinstrument'))
    return file.get(location)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def figureRowColumns(subplots):
    if subplots<1:
        raise AttributeError('Negative or zero number of subplots requiested.')
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


def centeroidnp(arr): # Calcualted centroid
        length = arr.shape[0]
        Totsum = np.sum(arr,axis=0)
        return Totsum/length

def compareNones(first,second,margin): # Function to compare
        if first.dtype == type(None) and second.dtype == type(None):
            return True
        elif first.dtype == second.dtype:
            return np.isclose(first,second,atol=margin)
        else:
            return False

#________________________________________________TESTS_____________________________________________

def test_DataSet_Creation():

    dataset = DataSet(OtherSetting=10.0)
    
    if(dataset.settings['OtherSetting']!=10.0):
        assert False


def test_Dataset_Initialization():

    emptyDataset = DataSet()
    del emptyDataset
    dataset = DataSet(OhterSetting=10.0,dataFiles='TestData/cameasim2018n000001.h5',convertedFiles='TestData/cameasim2018n000001.nxs')
    assert(dataset.dataFiles[0].name=='cameasim2018n000001.h5')
    assert(dataset.convertedFiles[0].name=='cameasim2018n000001.nxs')
                                                                                                                 
def test_DataSet_Error():
    

    ds = DataSet()
    
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


    ds.dataFiles = 'TestData/VanNormalization.h5'

def test_DataSet_Equality():
    D1 = DataSet(dataFiles='TestData/VanNormalization.h5')#,convertedFiles=['TestData/VanNormalization.nxs'])
    assert(D1==D1)

def test_DataSet_SaveLoad():
    
    D1 = DataSet(dataFiles='TestData/VanNormalization.h5')#,convertedFiles = 'TestData/VanNormalization.nxs')

    temp = 'temporary.bin'

    D1.save(temp)
    D2 = DataSet()
    D2.load(temp)
    os.remove(temp)
    assert(D1==D2) 

def test_DataSet_str():
    D1 = DataSet(dataFiles='TestData/cameasim2018n000001.h5',normalizationfiles = 'TestData/VanNormalization.h5')
    string = str(D1)
    print(string)


def test_DataSet_Convert_Data():
    dataFiles = 'TestData/cameasim2018n000001.h5'
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
    dataset.convertDataFile(dataFiles=dataFiles,binning=8,saveLocation='TestData/')
    convertedFile = dataset.convertedFiles[0]
    otherFile = DataFile.DataFile(dataFiles.replace('.h5','.nxs'))
    assert(convertedFile==otherFile)

    os.remove('TestData/cameasim2018n000001.nxs')
    


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
    DataFile = ['TestData/cameasim2018n000001.h5']

    dataset = DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='TestData/')

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    
    warnings.simplefilter('ignore')
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')
    viewer = MJOLNIR.Data.Viewer3D.Viewer3D(Intensity,bins)
    
    os.remove('TestData/cameasim2018n000001.nxs')
    del viewer

def test_DataSet_Visualization():
    import warnings
    from MJOLNIR.Data import Viewer3D
    DataFile = ['TestData/cameasim2018n000001.h5']

    dataset = DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='TestData/')

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    plt.ioff()
    warnings.simplefilter('ignore')
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')
    viewer = Viewer3D.Viewer3D(Intensity,bins)
    viewer.caxis = (0,100)
    plt.plot()
    plt.close()

def test_DataSet_binEdges():
    X = np.random.rand(100)*3 # array between 0 and 3 -ish
    X.sort()
    tolerance = 0.01
    Bins = binEdges(X,tolerance=tolerance)

    assert(Bins[0]==X[0]-0.5*tolerance)
    assert(Bins[-1]==X[-1]+0.5*tolerance)
    assert(len(Bins)<=3.0/tolerance)
    assert(np.all(np.diff(Bins[:-1])>tolerance))

def test_DataSet_1Dcut():
    q1 =  np.array([0,0.0])
    q2 =  np.array([3.0, 0.0])
    width = 0.1

    plt.ioff()
    convertFiles = ['TestData/cameasim2018n000011.h5']
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    ax,D,P,binCenter,binDistance = Datset.plotCut1D(q1,q2,width,minPixel=0.01,Emin=5.5,Emax=6.0,fmt='.')
    D2,P2 = Datset.cut1D(q1,q2,width,minPixel=0.01,Emin=5.5,Emax=6.0)
    assert(np.all([np.all(D[i]==D2[i]) for i in range(len(D))]))
    assert(np.all([np.all(P[i]==P2[i]) for i in range(len(P))]))

def test_DataSet_2Dcut():
    q1 =  np.array([0,0.0])
    q2 =  np.array([3.0, 0.0])
    width = 0.1
    minPixel=0.02
    EnergyBins = np.linspace(4,7,4)
    plt.ioff()
    convertFiles = ['TestData/cameasim2018n000011.h5']
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    ax,Data,pos,cpos,distance = Datset.plotCutQE(q1,q2,width,minPixel,EnergyBins,vmin=0.0 , vmax= 5e-06)
    Data2,pos2,cpos2,distance2 = Datset.cutQE(q1,q2,width,minPixel,EnergyBins)
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

def test_DataSet_cutPowder():
    Tolerance = 0.01

    plt.ioff()
    convertFiles = ['TestData/cameasim2018n000011.h5']
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    eBins = binEdges(Datset.energy,0.25)

    ax,D,q = Datset.plotCutPowder(eBins,Tolerance,vmin=0,vmax=1e-6)
    D2,q2 = Datset.cutPowder(eBins,Tolerance)
    for i in range(len(D)):
        for j in range(len(D[i])):
            assert(np.all(D[i][j]==D2[i][j]))

    for i in range(len(q)):
        for j in range(len(q[i])):
            assert(np.all(q[i][j]==q2[i][j]))


def test_DataSet_plotQPlane():
    plt.ioff()
    convertFiles = ['TestData/cameasim2018n000011.h5']
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    EMin = np.min(Datset.energy)
    EMax = EMin+0.5
    ax1 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=True,log=False,ax=None,RLUPlot=True)
    ax2 = Datset.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=True,ax=None,RLUPlot=True)

    ax1.set_clim(-20,-15)

    try:
        Datset.plotQPlane(EMin,EMax,binning='notABinningMethod')
        assert False
    except:
        assert True

def test_DataSet_plotA3A4():
    plt.ioff()

    File1 = 'TestData/T0Phonon10meV.nxs'
    File2 = 'TestData/T0Phonon10meV93_5A4.nxs'

    DS = DataSet(convertedFiles=[File1,File2])

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
        plotA3A4(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
        assert False
    except AttributeError:
        assert True

    try:
        plotA3A4(files,planes=[10,[22]],ax=axes,singleFigure=True) # 2 axes and singleFigure true
        assert False
    except AttributeError:
        assert True

    plotA3A4(files,planes=[10,[22,23]],ax=axes) # Plot plane 10 and 22+23 in the provided axes
    DS.plotA3A4(planes=[19,[22,25]]) # Plot planes in new axes
    DS.plotA3A4([F1,F1],planes=[19,[22,25]]) # Plot planes in new axes
    plt.close('all')


def test_DataSet_plotQPatches():
    plt.ioff()

    File1 = 'TestData/T0Phonon10meV.nxs'
    File2 = 'TestData/T0Phonon10meV93_5A4.nxs'

    DS = DataSet(convertedFiles=[File1,File2])

    F1 = DS.convertedFiles[0]
    F2 = DS.convertedFiles[1]

    files = [F1,F2]
    axes = [plt.figure().gca(),plt.figure().gca()]
    try:
        plotQPatches(files,planes=[],ax=axes) # 64 planes and only 2 axes
        assert False
    except AttributeError:
        assert True
        
    try:
        plotQPatches(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
        assert False
    except AttributeError:
        assert True

    try:
        plotQPatches(files,planes=[10,[22]],ax=axes,singleFigure=True) # 2 axes and singleFigure true
        assert False
    except AttributeError:
        assert True

    plotQPatches(files,planes=[10,[22,23]],ax=axes) # Plot plane 10 and 22+23 in the provided axes
    DS.plotQPatches(planes=[19,[22,25]],A4Extend=0.5,A3Extend=1) # Plot planes in new axes
    DS.plotQPatches(files=[files[0],files[0]],planes=[19,[22,25]],A4Extend=0.5,A3Extend=1) # Plot planes in new axes and only one file
    plt.close('all')

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
=======
    import matplotlib.pyplot as plt
    import os
    plt.ioff()
    DataFile = ['TestData/cameasim2018n000001.h5']

    dataset = DataSet(datafiles=DataFile)
    dataset.ConvertDatafile(savelocation='TestData/')

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    
    warnings.simplefilter('ignore')
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')
    viewer = MJOLNIR.Data.Viewer3D.Viewer3D(Intensity,bins)
    
    os.remove('TestData/cameasim2018n000001.nxs')
    del viewer

def test_DataSet_Visualization():
    import warnings
    from MJOLNIR.Data import Viewer3D
    DataFile = ['TestData/cameasim2018n000001.h5']

    dataset = DataSet(datafiles=DataFile)
    dataset.ConvertDatafile(savelocation='TestData/')

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    plt.ioff()
    warnings.simplefilter('ignore')
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')
    viewer = Viewer3D.Viewer3D(Intensity,bins)
    viewer.caxis = (0,100)
    plt.plot()
    plt.close()
