import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import scipy
from scipy.ndimage import filters
import matplotlib.pyplot as plt
import numpy as np
import pickle as pickle
import h5py as hdf
import scipy.optimize
import datetime
import warnings

dataLocation = 'entry/data/intensity'#'entry/Detectors/Detectors'
EiLocation = 'entry/data/incident_energy' # 'entry/Ei'
monLocation = 'entry/control/data'#'entry/Monitor'





class DataSet(object):
    def __init__(self, datafiles=None, normalizationfiles=None, calibrationfiles=None, convertedfiles=None, **kwargs):
        """DataSet object to hold all informations about data.
        
        Kwargs:
            
            - datafiles (list of strings): List of datafiles to be used in conversion (default None).

            - normalizationfiles (string or list of strings): Location of Vanadium normalization file(s) (default None).

            - calibrationfiles (string or list of strings): Location of calibration normalization file(s) (default None).

            - convertedfiles (string or list of strings): Location of converted data files (default None).

        Raises:

            - ValueError
            
            - NotImplementedError
        
        """
        
        self._datafiles = []
        self._normalizationfiles = []
        self._convertedfiles = []
        self._calibrationfiles = []


        if datafiles is not None:
            self.datafiles = datafiles

        if normalizationfiles is not None:
            self.normalizationfiles = normalizationfiles
        
        if convertedfiles is not None:
            self.convertedfiles = convertedfiles

        if calibrationfiles is not None:
            self.calibrationfiles = calibrationfiles


        self._settings = {}
            
        
        
        # Add all other kwargs to settings
        for key in kwargs:
            self.settings[key]=kwargs[key]
        



    @property
    def datafiles(self):
        return self._datafiles

    @datafiles.getter
    def datafiles(self):
        return self._datafiles

    @datafiles.setter
    def datafiles(self,datafiles):
        try:
            self._datafiles = isListOfStrings(datafiles)
        except Exception as e:
            #print('Error {} while parsing input.!'.format(e))
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
            #print('Error {} while parsing input.!'.format(e))
            raise(e)


    @property
    def convertedfiles(self):
        return self._convertedfiles

    @convertedfiles.getter
    def convertedfiles(self):
        return self._convertedfiles

    @convertedfiles.setter
    def convertedfiles(self,convertedfiles):
        try:
            self._convertedfiles = isListOfStrings(convertedfiles)
        except Exception as e:
            #print('Error {} while parsing input.!'.format(e))
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
            #print('Error {} while parsing input.!'.format(e))
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




    def convertDatafile(self,datafiles=None,binning=8,savelocation=None):
        """Conversion method for converting scan file(s) to hkl file. Converts the given h5 file into NXsqom format and saves in a file with same name, but of type .nxs.
        Copies all of the old data file into the new to ensure complete reduncency. Determins the binning wanted from the file name of normalization file.

        Kwargs:

            - datafiles (string or list of): File path(s), file must be of hdf format (default DataSet.datafiles).

            - binning (int): Binning to be used when converting files (default 8).

            - savelocation (string): File path to save location of data file(s) (defaults to same as raw file).

        Raises:

            - IOError

            - AttributeError
            
        """


        if datafiles is None:
            if len(self.datafiles)==0:
                raise AttributeError('No data files file provided either through input of in the DataSet object.')
            datafiles = self.datafiles


        if not isinstance(datafiles,list):
            datafiles=[datafiles]
        for datafile in datafiles:
            
            file = hdf.File(datafile,mode='r+')             
            instrument = getInstrument(file)
            
            
            if instrument.name.split('/')[-1] == 'CAMEA':
                EPrDetector = 8 # <----------_ REDO
            elif instrument.name.split('/')[-1] == 'MULTIFLEXX':
                EPrDetector = 1
            
            
            
            normalization = np.array(file.get('entry/calibration/{}_pixels'.format(binning)))
            if(normalization.shape==()):
                raise AttributeError('Binning not found in data file ({})'.format(binning))

            
            
            Data = np.array(instrument.get('detector/data'))

            detectors = Data.shape[1]
            
            if instrument.name.split('/')[-1] == 'MULTIFLEXX':
                Data.shape = (Data.shape[0],Data.shape[1],1)

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
            #print(QX.shape)
            #print(QY)
            #print(QX.shape.count(1))
            
            #if QX.shape.count(1)!=2:
            #    raise ValueError('At least two parameters changed simulatneously!')
            
            EnergyShape = (1,len(Ei),1,EfMean.shape[0],EfMean.shape[1])
            DeltaE = (Ei.reshape((Ei.shape[0],1,1))-EfMean.reshape((1,EfMean.shape[0],EfMean.shape[1]))).reshape(EnergyShape)
            
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
            
            if not savelocation is None:
                if savelocation[-1]!='/':
                    savelocation+='/'
                saveloc = savelocation+datafile.replace('.h5','.nxs').split('/')[-1]
            else:
                saveloc = datafile.replace('.h5','.nxs')
            saveNXsqom(datafile,file,saveloc,Intensity,Monitor,QX,QY,DeltaE,binning,Normalization)
            
            file.close()
            self.convertedfiles.append(saveloc)

    def binData3D(self,dx,dy,dz,datafiles=None):
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
        
        if datafiles is None:
            if len(self.convertedfiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                datafiles = self.convertedfiles
        elif not isinstance(datafiles,list):
            datafiles = [datafiles]
        #else:
            #raise AttributeError('datafiles attribute not understood.')

        #returnData = []
        I = []
        posx = []
        posy = []
        energy = []
        Norm = []
        Monitor = []

        for data in datafiles:
            
            file = hdf.File(data,'r')

            I.append(np.array(file.get('entry/data/intensity')))
            posx.append(np.array(file.get('entry/data/qx')))
            posy.append(np.array(file.get('entry/data/qy')))
            energy.append(np.array(file.get('entry/data/en')))
            Norm.append(np.array(file.get('entry/data/normalization')))
            Monitor.append(np.array(file.get('entry/data/monitor')))
            file.close()

        I = np.concatenate(I)
        posx = np.concatenate(posx)
        posy = np.concatenate(posy)
        energy = np.concatenate(energy)
        Norm = np.concatenate(Norm)
        Monitor = np.concatenate(Monitor)

        pos=[posx,posy,energy]

        returnData,bins = binData3D(dx,dy,dz,pos,I,norm=Norm,mon=Monitor)

        return returnData,bins

    def cut1D(self,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False,datafiles=None):
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

            - datafiles (list): List of datafiles to cut (default None). If none, the ones in the object will be used.
        
        
        Returns:
            
            - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
            
            - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
            
        """
        if datafiles is None:
            if len(self.convertedfiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                datafiles = self.convertedfiles
        elif not isinstance(datafiles,list):
            datafiles = [datafiles]
        I = []
        posx = []
        posy = []
        energy = []
        Norm = []
        Monitor = []

        for data in datafiles:
            
            file = hdf.File(data,'r')

            I.append(np.array(file.get('entry/data/intensity')))
            posx.append(np.array(file.get('entry/data/qx')))
            posy.append(np.array(file.get('entry/data/qy')))
            energy.append(np.array(file.get('entry/data/en')))
            Norm.append(np.array(file.get('entry/data/normalization')))
            Monitor.append(np.array(file.get('entry/data/monitor')))
            file.close()

        I = np.concatenate(I)
        posx = np.concatenate(posx)
        posy = np.concatenate(posy)
        energy = np.concatenate(energy)
        Norm = np.concatenate(Norm)
        Monitor = np.concatenate(Monitor)

        goodPixels = np.logical_and((1- np.isnan(Norm)).astype(bool),Norm!=0)

        I = I[goodPixels]
        qx=posx[goodPixels]
        qy=posy[goodPixels]
        energy=energy[goodPixels]
        Norm = Norm[goodPixels]
        Monitor = Monitor[goodPixels]
        positions = [qx,qy,energy]
        
        return cut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False)

    def plotCut1D(self,q1,q2,width,minPixel,Emin,Emax,ax=None,plotCoverage=False,datafiles=None,**kwargs):  # pragma: no cover
        """Plotting wrapper for the cut1D method. Generates a 1D plot with bins at positions corresponding to the distance from the start point. 
    Adds the 3D position on the x axis with ticks.
    
    .. note::
        Can only perform cuts for a constant energy plane of definable width.
    
    Kwargs:
        
        - ax (matplotlib axis): Figure axis into which the plots should be done (default None). If not provided, a new figure will be generated.
        
        - kwargs: All other keywords will be passed on to the ax.errorbar method.

        - datafiles (list): List of datafiles to cut (default None). If none, the ones in the object will be used.
    
    Returns:
        
        - ax (matplotlib axis): Matplotlib axis into which the plot was put.
        
        - Data list (4 arrays): Intensity, monitor count, normalization and normalization counts binned in the 1D cut.
        
        - Bin list (3 arrays): Bin edge positions in plane of size (n+1,3), orthogonal positions of bin edges in plane of size (2,2), and energy edges of size (2).
        
        - binCenter (3D array): Array containing the position of the bin centers of size (n,3)
        
        - binDistance (array): Distance from centre of bins to start position.
        """
        
        
        if datafiles is None:
            if len(self.convertedfiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                datafiles = self.convertedfiles
        elif not isinstance(datafiles,list):
            datafiles = [datafiles]
        I = []
        posx = []
        posy = []
        energy = []
        Norm = []
        Monitor = []

        for data in datafiles:
            
            file = hdf.File(data,'r')

            I.append(np.array(file.get('entry/data/intensity')))
            posx.append(np.array(file.get('entry/data/qx')))
            posy.append(np.array(file.get('entry/data/qy')))
            energy.append(np.array(file.get('entry/data/en')))
            Norm.append(np.array(file.get('entry/data/normalization')))
            Monitor.append(np.array(file.get('entry/data/monitor')))
            file.close()

        I = np.concatenate(I)
        posx = np.concatenate(posx)
        posy = np.concatenate(posy)
        energy = np.concatenate(energy)
        Norm = np.concatenate(Norm)
        Monitor = np.concatenate(Monitor)

        goodPixels = np.logical_and((1- np.isnan(Norm)).astype(bool),Norm!=0)

        I = I[goodPixels]
        qx=posx[goodPixels]
        qy=posy[goodPixels]
        energy=energy[goodPixels]
        Norm = Norm[goodPixels]
        Monitor = Monitor[goodPixels]
        positions = [qx,qy,energy]

        return plotCut1D(positions,I,Norm,Monitor,q1,q2,width,minPixel,Emin,Emax,ax,plotCoverage,**kwargs)


    def cutQE(self,q1,q2,width,minPixel,Emin,Emax,plotCoverage=False,datafiles=None):
        """Wrapper for cut data into maps of q and intensity between two q points and given energies. This is performed by doing consecutive constant energy planes.

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
        if datafiles is None:
            if len(self.convertedfiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                datafiles = self.convertedfiles
        elif not isinstance(datafiles,list):
            datafiles = [datafiles]
        I = []
        posx = []
        posy = []
        energy = []
        Norm = []
        Monitor = []

        for data in datafiles:
            
            file = hdf.File(data,'r')

            I.append(np.array(file.get('entry/data/intensity')))
            posx.append(np.array(file.get('entry/data/qx')))
            posy.append(np.array(file.get('entry/data/qy')))
            energy.append(np.array(file.get('entry/data/en')))
            Norm.append(np.array(file.get('entry/data/normalization')))
            Monitor.append(np.array(file.get('entry/data/monitor')))
            file.close()

        I = np.concatenate(I)
        posx = np.concatenate(posx)
        posy = np.concatenate(posy)
        energy = np.concatenate(energy)
        Norm = np.concatenate(Norm)
        Monitor = np.concatenate(Monitor)

        goodPixels = np.logical_and((1- np.isnan(Norm)).astype(bool),Norm>0)

        I = I[goodPixels]
        qx=posx[goodPixels]
        qy=posy[goodPixels]
        energy=energy[goodPixels]
        Norm = Norm[goodPixels]
        Monitor = Monitor[goodPixels]
        positions = [qx,qy,energy]
        
        return cutQE(positions,I,Norm,Monitor,q1,q2,width,minPix,EnergyBins)

    def plotCutQE(self,q1,q2,width,minPixel,EnergyBins,ax=None,plotCoverage=False,datafiles=None,**kwargs): # pragma: no cover
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
        
        
        if datafiles is None:
            if len(self.convertedfiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                datafiles = self.convertedfiles
        elif not isinstance(datafiles,list):
            datafiles = [datafiles]
        I = []
        posx = []
        posy = []
        energy = []
        Norm = []
        Monitor = []

        for data in datafiles:
            
            file = hdf.File(data,'r')

            I.append(np.array(file.get('entry/data/intensity')))
            posx.append(np.array(file.get('entry/data/qx')))
            posy.append(np.array(file.get('entry/data/qy')))
            energy.append(np.array(file.get('entry/data/en')))
            Norm.append(np.array(file.get('entry/data/normalization')))
            Monitor.append(np.array(file.get('entry/data/monitor')))
            file.close()

        I = np.concatenate(I)
        posx = np.concatenate(posx)
        posy = np.concatenate(posy)
        energy = np.concatenate(energy)
        Norm = np.concatenate(Norm)
        Monitor = np.concatenate(Monitor)

        goodPixels = np.logical_and((1- np.isnan(Norm)).astype(bool),Norm>0)

        I = I[goodPixels]
        qx=posx[goodPixels]
        qy=posy[goodPixels]
        energy=energy[goodPixels]
        Norm = Norm[goodPixels]
        Monitor = Monitor[goodPixels]
        positions = [qx,qy,energy]

        return plotCutQE(positions,I,Norm,Monitor,q1,q2,width,minPixel,EnergyBins,ax = None,**kwargs)


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
    dirvec = np.array(q2)-np.array(q1)
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
    return bin_edges

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
    
    
    def format_coord(x,y,binDistance,binCenter):
        index = np.argmin(np.abs(binDistance-x))
        qx,qy,E = binCenter[index]
        return  "qx = {0:.3f}, qy = {1:.3f}, E = {2:.3f}, I = {3:0.4e}".format(qx,qy,E,y)
    
    ax.format_coord = lambda x,y: format_coord(x,y,binDistance,binCenter)

    return ax,D,P,binCenter,binDistance



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
    for i in range(len(Int)):
        ax.pcolormesh(binedges[i],binenergies[i],Int[i].T,**kwargs)
    
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
    
    
    def format_coord(x,y,binDistance,centerPos,Int):
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

    return ax,[intensityArray,monitorArray,normalizationArray,normcountArray],returnpositions,centerPos,binDistance


        
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
    
    rawdata = proc.create_dataset('rawdata',shape=(1,),dtype='S200',data=np.string_(datafile))
    rawdata.attrs['NX_class']=b'NX_CHAR'
    
    normalizationString = proc.create_dataset('binning',shape=(1,),dtype='int32',data=binning)
    normalizationString.attrs['NX_class']=b'NX_INT'
    
    data = fd.get('entry/data')
    #data['rawdata']=data['data']
    #del data['data']
    
    
    fileLength = Intensity.size
    
    Int = data.create_dataset('intensity',shape=(fileLength,),dtype='int32',data=Intensity.flatten())
    Int.attrs['NX_class']='NX_INT'
    
    monitor = data.create_dataset('monitor',shape=(fileLength,),dtype='int32',data=Monitor.flatten())
    monitor.attrs['NX_class']=b'NX_INT'

    normalization = data.create_dataset('normalization',shape=(fileLength,),dtype='float32',data=Normalization.flatten())
    normalization.attrs['NX_class']=b'NX_FLOAT'
    
    qx = data.create_dataset('qx',shape=(fileLength,),dtype='float32',data=QX.flatten())
    qx.attrs['NX_class']=b'NX_FLOAT'
    
    qy = data.create_dataset('qy',shape=(fileLength,),dtype='float32',data=QY.flatten())
    qy.attrs['NX_class']=b'NX_FLOAT'
    
    qz = data.create_dataset('qz',shape=(fileLength,),dtype='float32',data=np.zeros((fileLength,)))
    qz.attrs['NX_class']=b'NX_FLOAT'
    
    en = data.create_dataset('en',shape=(fileLength,),dtype='float32',data=DeltaE.flatten())
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


# def binEdges(values,tolerance):
#     """Generate bins for a value array and bin points within tolerance. If tolerance is set too high, only the outer bin edges are returned.
#     Low tolerance results in all points having their own bin.

#     Args:

#         - values (array): Value array to be binned

#         - tolerance (float): Tolerance in with points are binned togehter.

#     returns:

#         - binEdges (array): Edges starting tolerance*0.5 below lowest value and ends tolerance*0.5 above largest value.

#     """
#     values_array = np.array(values).ravel()
#     unique_values = np.asarray(list(set(values_array)))
#     unique_values.sort()
#     if len(unique_values)==0: # If no values are given, return empty list
#         return []
#     bin_edges = [unique_values[0] - tolerance / 2]
#     for i in range(len(unique_values) - 1):
#         if unique_values[i+1] - unique_values[i] > tolerance:
#             bin_edges.append((unique_values[i] + unique_values[i+1]) / 2)
#         else:
#             pass
    
#     bin_edges.append(unique_values[-1] + tolerance / 2)
#     return bin_edges



#________________________________________________TESTS_____________________________________________

def test_DataSet_Creation():

    dataset = DataSet(OtherSetting=10.0)
    
    if(dataset.settings['OtherSetting']!=10.0):
        assert False


def test_Dataset_Initialization():

    emptyDataset = DataSet()
    del emptyDataset
    dataset = DataSet(OhterSetting=10.0,datafiles='SomeFile',convertedfiles='Converted.nxs')
    assert(dataset.datafiles[0]=='SomeFile')
    assert(dataset.convertedfiles[0]=='Converted.nxs')


def test_DataSet_Error():
    

    ds = DataSet()
    
    try: # Wrong data file type
        ds.datafiles = 100
        assert False
    except AttributeError:
        assert True


    try: # Can't overwrite settings
        ds.settings={}
        assert False
    except NotImplementedError:
        assert True

    try:# Wrong data file type
        ds.convertedfiles = 10
        assert False
    except AttributeError:
        assert True

    ds.datafiles = 'Here.h5'



def test_DataSet_Equality():
    D1 = DataSet(datafiles='Here',convertedfiles=['Van/1.nxs','Van/2.nxs'])
    assert(D1==D1)

def test_DataSet_SaveLoad():
    
    D1 = DataSet(datafiles='Here',convertedfiles = 'Van.nxs')

    temp = 'temporary.bin'

    D1.save(temp)
    D2 = DataSet()
    D2.load(temp)
    os.remove(temp)
    assert(D1==D2) 

def test_DataSet_str():
    D1 = DataSet(datafiles='Here',normalizationfiles = 'Van.nxs')
    string = str(D1)
    print(string)


def test_DataSet_Convert_Data():

    DataFiles = 'TestData/cameasim2018n000001.h5'
    dataset = DataSet(datafiles=DataFiles)
    

    try:
        dataset.convertDatafile(datafiles=DataFiles,binning=100)
        assert False
    except AttributeError: # Cant find normalization table
        assert True

    dataset.convertDatafile(datafiles=DataFiles,binning=8,savelocation='TestData/')
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

    dataset = DataSet(datafiles=DataFile)
    dataset.convertDatafile(savelocation='TestData/')

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
    dataset.convertDatafile(savelocation='TestData/')

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
    assert(np.all(np.diff(Bins)>tolerance))

def test_DataSet_1Dcut():
    q1 =  np.array([0,0.0])
    q2 =  np.array([3.0, 0.0])
    width = 0.1

    plt.ioff()
    convertFiles = ['TestData/cameasim2018n000011.h5']
    
    Datset = DataSet(datafiles = convertFiles)
    Datset.convertDatafile()
    ax,D,P,binCenter,binDistance = Datset.plotCut1D(q1,q2,width,minPixel=0.01,Emin=5.5,Emax=6.0,fmt='.')

def test_DataSet_2Dcut():
    q1 =  np.array([0,0.0])
    q2 =  np.array([3.0, 0.0])
    width = 0.1
    minPixel=0.02
    EnergyBins = np.linspace(4,7,4)
    plt.ioff()
    convertFiles = ['TestData/cameasim2018n000011.h5']
    
    Datset = DataSet(datafiles = convertFiles)
    Datset.convertDatafile()
    ax3,DD3,pos3,cpos3,distance3 = Datset.plotCutQE(q1,q2,width,minPixel,EnergyBins,vmin=0.0 , vmax= 5e-06)

