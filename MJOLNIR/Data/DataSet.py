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

dataLocation = 'entry/data/data'#'entry/Detectors/Detectors'
EiLocation = 'entry/data/incident_energy' # 'entry/Ei'
monLocation = 'entry/control/data'#'entry/Monitor'





class DataSet(object):
    def __init__(self, datafiles=None,normalizationfiles=None, calibrationfiles=None, convertedfiles=None, **kwargs):
        """DataSet object to hold all informations about data.
        
        Kwargs:
            
             - datafiles (list of strings): List of datafiles to be used in conversion.

            - normalizationfiles (string or list of strings): Location of Vanadium normalization file(s).

            - calibrationfiles (string or list of strings): Location of calibration normalization file(s).

            - convertedfiles (string or list of strings): Location of converted data files.

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
            self._datafiles = IsListOfStrings(datafiles)
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
            self._normalizationfiles = IsListOfStrings(normalizationfiles)
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
            self._convertedfiles = IsListOfStrings(convertedfiles)
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
            self._calibrationfiles = IsListOfStrings(calibrationfiles)
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




    def ConvertDatafile(self,datafiles=None,binning=8,savelocation=None):
        """Conversion method for converting scan file(s) to hkl file. Converts the given h5 file into NXsqom format and saves in a file with same name, but of type .nxs.
        Copies all of the old data file into the new to ensure complete reduncency. Determins the binning wanted from the file name of normalization file.

        Kwargs

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

            A4 = np.deg2rad(np.array(normalization[:,9]))#-A4Zero
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
            
            A3 = np.deg2rad(np.array(file.get('/entry/sample/rotation_angle/')))#+A3Zero
            
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

            I.append(np.array(file.get('entry/data/data')))
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
        






        
def IsListOfStrings(object):
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
    data['rawdata']=data['data']
    del data['data']
    
    
    fileLength = Intensity.size
    
    Int = data.create_dataset('data',shape=(fileLength,),dtype='int32',data=Intensity.flatten())
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

        Rebinned intensity (and if provided Normalization, Monitor, and Normalization Count) and X, Y, and Z bins in 3 1D arrays.


    Example:

    >>> pos = [Qx,Qy,E]
    >>> Data,bins = DataSet.binData3D(0.05,0.05,0.2,pos,I,norm=Norm,mon=Monitor)

    """

    if bins is None:
        bins = calculateBins(dx,dy,dz,pos)
    
    #NonNaNs = 1-np.isnan(data.flatten())

    #pos = [np.array(x[NonNaNs]) for x in pos]

    intensity =    np.histogramdd(np.array(pos).T,bins=bins,weights=data.flatten())[0].astype(data.dtype)

    returndata = [intensity]
    if mon is not None:
        MonitorCount=  np.histogramdd(np.array(pos).T,bins=bins,weights=mon.flatten())[0].astype(mon.dtype)
        returndata.append(MonitorCount)
    if norm is not None:
        Normalization= np.histogramdd(np.array(pos).T,bins=bins,weights=norm.flatten())[0].astype(norm.dtype)
        NormCount =    np.histogramdd(np.array(pos).T,bins=bins,weights=np.ones_like(data).flatten())[0].astype(int)
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
    
    bins=[XX[:,0,0],YY[0,:,0],ZZ[0,0,:]]
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

#________________________________________________TESTS_____________________________________________

def test_DataSet_Creation():

    dataset = DataSet(OtherSetting=10.0)
    
    if(dataset.settings['OtherSetting']!=10.0):
        assert False


def test_Dataset_Initialization():

    emptyDataset = DataSet()
        
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
        dataset.ConvertDatafile(datafiles=DataFiles,binning=100)
        assert False
    except AttributeError: # Cant find normalization table
        assert True

    dataset.ConvertDatafile(datafiles=DataFiles,binning=8,savelocation='TestData/')
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
    assert(np.all(bins[0]==np.linspace(-0.25,1.25,4)))
    assert(np.all(bins[1]==np.linspace(-0.125,1.125,6)))
    assert(np.all(bins[2]==np.linspace(1-0.125,2.125,6)))

    assert(RebinnedNorm.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.dtype==int)
    assert(RebinnedNorm.dtype==Norm.dtype)
    assert(ReBinnedI.dtype==I.dtype)


def test_DataSet_full_test():
    import MJOLNIR.Data.Viewer3D
    import warnings
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