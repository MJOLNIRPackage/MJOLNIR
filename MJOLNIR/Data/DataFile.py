import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import scipy
import matplotlib.pyplot as plt
import numpy as np
import h5py as hdf
import warnings
from MJOLNIR import _tools
import datetime
import math
import shapely
from shapely.geometry import Polygon as PolygonS, Point as PointS
from MJOLNIR import TasUBlibDEG as TasUBlib
from MJOLNIR._tools import Marray

factorsqrtEK = 0.694692

def cosd(x):
    return np.cos(np.deg2rad(x))

def sind(x):
    return np.sin(np.deg2rad(x))

class DataFile(object):
    """Object to load and keep track of HdF files and their conversions"""
    def __init__(self,fileLocation=None):
        # Check if file exists
        if isinstance(fileLocation,DataFile): # Copy everything in provided file
            self.updateProperty(fileLocation.__dict__)
        elif isinstance(fileLocation,str) :
            if not os.path.isfile(fileLocation):
                raise AttributeError('File location does not exist({}).'.format(fileLocation))
            if fileLocation.split('.')[-1]=='nxs':
    	        self.type='nxs'
        
            elif fileLocation.split('.')[-1]=='hdf':
	            self.type='hdf'
            else:
                raise AttributeError('File is not of type nxs or hdf.')
            self.name = fileLocation.split('/')[-1]
            self.fileLocation = os.path.abspath(fileLocation)		
			
            with hdf.File(fileLocation) as f:
			
                sample=f.get('/entry/sample')
                self.sample = Sample(sample=f.get('/entry/sample'))
                instr = getInstrument(f)
                if self.type == 'hdf':
                    self.I = Marray(instr.get('detector/counts')).swapaxes(1,2)
                else:
                    self.I=Marray(f.get('entry/data/intensity'))
                    self.counts = Marray(instr.get('detector/counts')).swapaxes(1,2)
                    self.qx=Marray(f.get('entry/data/qx'))
                    self.qy=Marray(f.get('entry/data/qy'))
                    self.h=Marray(f.get('entry/data/h'))
                    self.k=Marray(f.get('entry/data/k'))
                    self.l=Marray(f.get('entry/data/l'))
                    self.energy=Marray(f.get('entry/data/en'))
                    self.Norm=Marray(f.get('entry/data/normalization'))
                self.MonitorMode = np.array(f.get('entry/control/mode'))[0].decode()
                self.MonitorPreset=np.array(f.get('entry/control/preset'))                
                if self.type == 'hdf':
                    self.Monitor = Marray(f.get('entry/control/data'))
                    if not self.MonitorMode == 't' and len(self.Monitor)>1: # If not counting on time and more than one point saved
                        if self.Monitor.flatten()[0]!=self.MonitorPreset: # For all data in 2018 with wrong monitor saved
                            self.Monitor = np.ones_like(self.Monitor)*self.MonitorPreset ### TODO: Make Mark save the correct monitor!!
                else:
                    self.Monitor=Marray(f.get('entry/data/monitor'))
                self.Time = np.array(f.get('entry/control/time'))

                instr = getInstrument(f)
                self.instrument = instr.name.split('/')[-1]
                self.possibleBinnings = np.array([int(x[-1]) for x in np.array(instr) if x[:5]=='calib'])
                self.Ei = Marray(instr.get('monochromator/energy'))
                self.A3 = Marray(f.get('entry/sample/rotation_angle'))
                self.A4 = Marray(instr.get('analyzer/polar_angle')).reshape(-1)
                try:
                    self.A3Off = self.sample.A3Off#np.array(f.get('entry/sample/rotation_angle_zero'))  
                except:
                    self.A3Off = [0.0]
                self.A4Off = Marray(instr.get('analyzer/polar_angle_offset'))
                if self.type == 'hdf':
                    self.binning=1 # Choose standard binning 1
                else:
                    self.binning = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/binning'))[0]
                calibrations = []
                for binning in self.possibleBinnings:
                    Ef = np.array(instr.get('calib{}/final_energy'.format(str(binning))))
                    width = np.array(instr.get('calib{}/width'.format(str(binning))))
                    bg = np.array(instr.get('calib{}/background'.format(str(binning))))
                    amp = np.array(instr.get('calib{}/amplitude'.format(str(binning))))
                    EfTable = np.array([amp,Ef,width,bg]).T
                    A4 = np.array(instr.get('calib{}/a4offset'.format(str(binning))))
                    bound = np.array(instr.get('calib{}/boundaries'.format(str(binning))))
                    calibrations.append([EfTable,A4,bound])
                self.instrumentCalibrations = np.array(calibrations)
                self.loadBinning(self.binning)
                
                self.temperature = np.array(sample.get('temperature'))
                self.magneticField = np.array(sample.get('magnetic_field'))
                self.electricField = np.array(sample.get('electric_field'))
                self.scanParameters,self.scanValues,self.scanUnits = getScanParameter(f)
                self.scanCommand = np.array(f.get('entry/scancommand'))
                if self.type == 'nxs':
                    self.original_file = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/rawdata'))[0].decode()
                self.title = np.array(f.get('entry/title'))

                try:
                    self.scanParameters,self.scanValues,self.scanUnits = getScanParameter(f)
                except:
                    pass

                if self.type == 'hdf':
                    ###################
                    #self.I[:,:,:150]=0#
                    ###################
                    pass
                
            for key in ['magneticField','temperature','electricField']:
                if self.__dict__[key].dtype ==object: # Is np nan object
                    self.__dict__[key] = None
        else:
            pass
        

    @property
    def A3Off(self):
        return self._A3Off

    @A3Off.getter
    def A3Off(self):
        return self._A3Off

    @A3Off.setter
    def A3Off(self,A3Off):
        if A3Off is None:
            self._A3Off = 0.0
        else:
            self._A3Off = A3Off

    @property
    def A4Off(self):
        return self._A4Off

    @A4Off.getter
    def A4Off(self):
        return self._A4Off

    @A4Off.setter
    def A4Off(self,A4Off):
        if A4Off is None:
            self._A4Off = 0.0
        else:
            self._A4Off = A4Off

    @property
    def A3(self):
        return self._A3

    @A3.getter
    def A3(self):
        return self._A3

    @A3.setter
    def A3(self,A3):
        if A3.shape == ():
            self._A3 = [0.0]
        else:
            if A3[0] is None:
                self._A3 = [0.0]
            else:
                self._A3 = A3

    @property
    def A4(self):
        return self._A4

    @A4.getter
    def A4(self):
        return self._A4

    @A4.setter
    def A4(self,A4):
        if A4.shape == ():
            self._A4 = [0.0]
        else:
            self._A4 = A4

    def updateProperty(self,dictionary):
        if isinstance(dictionary,dict):
            for key in dictionary.keys():
                self.__setattr__(key,dictionary[key])

    def __eq__(self,other):
        print(self.difference(other))
        return len(self.difference(other))==0
    
    def difference(self,other,keys = set(['sample','instrument','Ei','I','_A3','_A4','binning','scanParameters'])):
        """Return the difference between two data files by keys"""
        dif = []
        if not set(self.__dict__.keys()) == set(other.__dict__.keys()): # Check if same generation and type (hdf or nxs)
            return list(set(self.__dict__.keys())-set(other.__dict__.keys()))

        comparisonKeys = keys
        for key in comparisonKeys:
            skey = self.__dict__[key]
            okey = other.__dict__[key]
            if isinstance(skey,np.ndarray):
                try:
                    if not np.all(np.isclose(skey,okey)):
                        if not np.all(np.isnan(skey),np.isnan(okey)):
                            dif.append(key)
                except (TypeError, AttributeError,ValueError):
                    if np.all(skey!=okey):
                        dif.append(key)
            elif not np.all(self.__dict__[key]==other.__dict__[key]):
                dif.append(key)
        return dif

    def __str__(self):
        returnStr = 'Data file {} from the MJOLNIR software package of type {}\n'.format(self.name,self.type)
        returnStr+= 'Ei: '+ str(self.Ei) + '\nA3: ' + ','.join([str(x) for x in self.A3])
        returnStr+= '\nA4: ' + ','.join([str(x) for x in self.A4]) + '\nSample: '+str(self.sample)
        return returnStr

    def __add__(self,other):
        raise NotImplementedError('Adding two data files is not yet supported.')

    def __hasattr__(self,s):
        return s in self.__dict__.keys()

    @_tools.KwargChecker()
    def convert(self,binning):
        if self.instrument == 'CAMEA':
            EPrDetector = 8 
        elif self.instrument in ['MULTIFLEXX','FLATCONE']:
            EPrDetector = 1
        else:
            raise AttributeError('Instrument type of data file not understood. {} was given.'.format(self.instrument))
        
        self.loadBinning(binning)
        
        EfNormalization = self.instrumentCalibrationEf
        A4Normalization = self.instrumentCalibrationA4#np.array(instrument.get('calib{}/a4offset'.format(str(binning))))
        EdgesNormalization = self.instrumentCalibrationEdges#np.array(instrument.get('calib{}/boundaries'.format(str(binning))))
        Data = self.I#np.array(instrument.get('detector/data'))
        

        detectors = Data.shape[1]
        steps = Data.shape[0]
        
        if self.instrument in ['MULTIFLEXX','FLATCONE']:
            Data.shape = (Data.shape[0],Data.shape[1],-1)

        A4Zero = self.A4Off#file.get('entry/sample/polar_angle_zero')
        
        
        if A4Zero is None:
            A4Zero=0.0
        else:
            A4Zero = np.array(A4Zero)

        A3Zero = self.A3Off
        if A3Zero is None:
            A3Zero=0.0
        else:
            A3Zero = np.deg2rad(np.array(A3Zero))

        A4 = np.deg2rad(A4Normalization)
        A4=A4.reshape(detectors,binning*EPrDetector,order='C')

        PixelEdge = EdgesNormalization.reshape(detectors,EPrDetector,binning,2).astype(int)
        A4File = self.A4
        
        A4File = A4File.reshape((-1,1,1))

        A4Mean = -(A4.reshape((1,detectors,binning*EPrDetector))+np.deg2rad(A4File-A4Zero))
        
        Intensity=np.zeros((Data.shape[0],Data.shape[1],EPrDetector*binning),dtype=int)
        for i in range(detectors): # for each detector
            for j in range(EPrDetector):
                for k in range(binning):
                    Intensity[:,i,j*binning+k] = np.sum(Data[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)

        EfMean = EfNormalization[:,1].reshape(1,A4.shape[0],EPrDetector*binning)
        EfNormalization = EfNormalization[:,0].reshape(1,A4.shape[0],EPrDetector*binning)#(EfNormalization[:,0]*np.sqrt(2*np.pi)*EfNormalization[:,2]).reshape(1,A4.shape[0],EPrDetector*binning)
        A3 = np.deg2rad(np.array(self.A3))+A3Zero #file.get('/entry/sample/rotation_angle/')
        if A3.shape[0]==1:
            A3 = A3*np.ones((steps))
        
        A3.resize((steps,1,1))
        Ei = self.Ei.reshape(-1,1,1)#np.array(instrument.get('monochromator/energy'))
        if False:
            kf = factorsqrtEK*np.sqrt(EfMean)#.reshape(1,detectors,binning*EPrDetector)
            
            ki = factorsqrtEK*np.sqrt(Ei).reshape(-1,1,1)
            # Shape everything into shape (steps,detectors,bins) (if external parameter 
            # is changed, this is assured by A3 reshape)
            Qx = ki-kf*np.cos(A4Mean)
            Qy = -kf*np.sin(A4Mean)
            QX = Qx*np.cos(A3)-Qy*np.sin(A3)
            QY = Qx*np.sin(A3)+Qy*np.cos(A3)
        else:
            UB = self.sample.orientationMatrix
            UBINV = np.linalg.inv(UB)
            HKL,QX,QY = TasUBlib.calcTasQH(UBINV,[np.rad2deg(A3).squeeze(),
            -np.rad2deg(A4Mean)],Ei.squeeze().reshape(-1,1,1),EfMean.squeeze())
            H,K,L = np.swapaxes(np.swapaxes(HKL,1,2),0,3)
            self.sample.B = TasUBlib.calculateBMatrix(self.sample.cell)

        DeltaE = Ei-EfMean
        if DeltaE.shape[0]==1:
            DeltaE = DeltaE*np.ones((steps,1,1))
        Monitor = self.Monitor.reshape((steps,1,1))
        Monitor = Monitor*np.ones((1,detectors,EPrDetector*binning))
        Normalization = EfNormalization*np.ones((steps,1,1))

       
        ###########################
        #Monitor[:,:,:binning] = 0 #
        ###########################


        convFile = DataFile(self) # Copy everything from old file
        updateDict = {'I':Intensity,'Monitor':Monitor,'qx':QX,'qy':QY,'energy':DeltaE,'binning':binning,'Norm':Normalization,
        'h':H,'k':K,'l':L,'type':'nxs','fileLocation':None,'original_file':self,'name':self.name.replace('.hdf','.nxs')}
        convFile.updateProperty(updateDict)
        return convFile


    @_tools.KwargChecker()
    def plotA4(self,binning=None):
        """Method to plot the fitted A4 values of the normalization table

        Kwargs:
    
            - binning (int): Binning for the corresponding normalization table (default self.binning or 8)

        returns:

            - fig (matplotlib figure): Figure into which the A4 values are plotted

        """
        self.loadBinning(binning)
        binning = self.binning
        Norm = (self.instrumentCalibrationEf[:,0]*self.instrumentCalibrationEf[:,2]*np.sqrt(2*np.pi)).reshape((104,8*binning))

        A4 = np.reshape(self.instrumentCalibrationA4,(104,8*binning))
        fig = plt.figure()
        for i in range(104):
            plt.scatter(-A4[i],np.arange(len(A4[i])),c=Norm[i])
        
        plt.title('Instrument calibration')
        plt.ylabel('Pixel')
        plt.xlabel('-A4 [deg]')

        return fig

    @_tools.KwargChecker()
    def plotEf(self,binning=None):
        """Method to plot the fitted Ef values of the normalization table

        Kwargs:
    
            - binning (int): Binning for the corresponding normalization table (default self.binning or 8)

        returns:

            - fig (matplotlib figure): Figure into which the Ef values are plotted

        """
        self.loadBinning(binning)
        
        binning = self.binning
        Ef = self.instrumentCalibrationEf[:,1].reshape(104,8*binning)
        fig = plt.figure()
        for i in range(104):
            plt.scatter(i*np.ones_like(Ef[i]),Ef[i],zorder=10)
        plt.xlabel('Detector number')
        plt.ylabel('Ef [meV]')
        plt.title('Final Energy Individual')
        plt.grid(True)

        return fig

    @_tools.KwargChecker()
    def plotEfOverview(self,binning=None):
        """Method to plot the fitted Ef values of the normalization table

        Kwargs:
    
            - binning (int): Binning for the corresponding normalization table (default self.binning or 8)

        returns:

            - fig (matplotlib figure): Figure into which the Ef values are plotted

        """
        self.loadBinning(binning)
        binning = self.binning
        Ef = self.instrumentCalibrationEf[:,1].reshape(104,8*binning)
        fig = plt.figure()
        plt.imshow(Ef.T,origin='lower')
        plt.xlabel('Detector number')
        plt.ylabel('Pixel number')
        plt.colorbar()
        plt.title('Final Energy Overview')

        return fig

    @_tools.KwargChecker()
    def plotNormalization(self,binning=None):
        """Method to plot the fitted integrated intensities of the normalization table

        Kwargs:
    
            - binning (int): Binning for the corresponding normalization table (default self.binning or 8)

        returns:

            - fig (matplotlib figure): Figure into which the Ef values are plotted

        """
        
        self.loadBinning(binning)
        binning = self.binning
        Norm = (self.instrumentCalibrationEf[:,0]*self.instrumentCalibrationEf[:,2]*np.sqrt(2*np.pi)).reshape((104,8*binning))

        fig = plt.figure()
        plt.imshow(Norm.T,origin='lower')
        plt.xlabel('Detector number')
        plt.ylabel('Pixel number')
        plt.title('Normalization')
        plt.colorbar()
        return fig

    @_tools.KwargChecker()
    def loadBinning(self,binning):
        """Small function to check if current binning is equal to wanted binning and if not reloads to binning wanted"""


        if binning is None or (binning == self.binning and hasattr(self,'instrumentCalibrationEf')):
            binning = self.binning
        else:
            if not binning in self.possibleBinnings:
                raise AttributeError('Wanted binning not in possible binnings!')
            binID = list(self.possibleBinnings).index(binning)
            
            self.instrumentCalibrationEf,self.instrumentCalibrationA4,self.instrumentCalibrationEdges = self.instrumentCalibrations[binID]
            self.binning = binning

    @_tools.KwargChecker()
    def calculateEdgePolygons(self,addEdge=True): # pragma: no cover
        """Method to calculate bounding polygon for all energies. The energies are split using the bin-edges method of DataSet. Hereafter,
        the outer most points are found in polar coordinates and a possible addition is made creating the padded bounding polygon.
        
        Kwargs:
            
            - addEdge (bool/float): If true, padding is found as difference between outer and next outer point. If addEdge is a number, generate padding a padding of this value (default True)
            
        Returns:
            
            - edgePolygon (list): List of shapely polygons of the boundary
            
            - EBins (list): Binning edges in energy
            
        """
        if addEdge:
            if np.isclose(float(addEdge),1.0):
                addEdgeAmount = None
            else:
                addEdgeAmount = addEdge
        
        Energy = self.energy
        
        E = np.sort(np.mean(self.energy,axis=(0,1)))
        dif = np.diff(E)*0.5
        EBins = np.concatenate([[E[0]-dif[0]],E[:-1]+dif,[E[-1]+dif[-1]]])
        steps = Energy.shape[0]
        
        
        edgePolygon = []
        for i in range(len(EBins)-1):
            ELow = EBins[i]
            EHigh= EBins[i+1]
            
            EBool = np.logical_and(Energy>ELow,Energy<EHigh)
            
            x = self.qx[EBool]
            y = self.qy[EBool]
            
            
            r = np.linalg.norm([x,y],axis=0)
            theta = np.arctan2(y,x)
            
            rBins = _tools.binEdges(r,0.00001)
            
            out = -1
            while np.sum(r>rBins[out])<steps:
                out-=1
            rOuter = r>rBins[out]
            inner = 0
            while np.sum(r<rBins[inner])<steps:
                inner+=1
            rInner = r<rBins[inner]
            
            minEdge = []
            maxEdge = []
            _minEdge= []
            _maxEdge= []
            include = 0
            for j in range(len(rBins)-1):
                Bool = np.logical_and(r>rBins[j-include],r<rBins[j+1])
                if np.sum(Bool)<steps-1:
                    include+=1
                    continue
                else:
                    include = 0
                TT = theta[Bool]
                dif = np.diff(TT)
                if np.max(abs(dif))>np.pi*1.9:
                    idx = np.argmax(abs(dif))
                    TT[idx+1:]+=-np.sign(dif[idx])*2*np.pi
                minT,maxT = _tools.minMax(TT)
                _minT,_maxT = _tools.minMax(TT)
                
                if addEdge:
                    mint, maxt = _tools.minMax(TT[np.logical_not(np.logical_or(np.isclose(TT,minT),np.isclose(TT,maxT)))])
                    if addEdgeAmount is None:
                        minT = minT-0.5*(mint-minT)
                        maxT = maxT-0.5*(maxt-maxT)
                    else:
                        minT = minT+np.sign(minT-mint)*addEdgeAmount
                        maxT = maxT+np.sign(maxT-maxt)*addEdgeAmount
                    
                R = np.mean(r[Bool])
                minEdge.append([R*np.cos(minT),R*np.sin(minT)])
                maxEdge.append([R*np.cos(maxT),R*np.sin(maxT)])
                
                _minEdge.append([R*np.cos(_minT),R*np.sin(_minT)])
                _maxEdge.append([R*np.cos(_maxT),R*np.sin(_maxT)])
            
            minEdge = np.array(minEdge).T
            maxEdge = np.array(maxEdge).T
            
            innerPoints = np.array([x[rInner],y[rInner]])
            _innerPoints= np.array([x[rInner],y[rInner]])
            if addEdge:
                if addEdgeAmount is None:
                    RR = rBins[inner]-(rBins[inner+1]-rBins[inner])
                else:
                    RR = rBins[inner]-addEdgeAmount
                Theta = np.arctan2(innerPoints[1],innerPoints[0])
                innerPoints = np.array([np.cos(Theta),np.sin(Theta)])*RR
                    
                
            outerPoints = np.array([x[rOuter],y[rOuter]])
            _outerPoints= np.array([x[rOuter],y[rOuter]])
            if addEdge:
                if addEdgeAmount is None:
                    RR = rBins[out]-(rBins[out-1]-rBins[out])
                else:
                    RR = rBins[out]+addEdgeAmount
                Theta = np.arctan2(outerPoints[1],outerPoints[0])
                outerPoints = np.array([np.cos(Theta),np.sin(Theta)])*RR
            
                
            refvec1,center1 = _tools.calRefVector(innerPoints)
            refvec2,center2 = _tools.calRefVector(outerPoints)
            sInnerEdge = np.array(sorted(innerPoints.T,key=lambda x: _tools.clockwiseangle_and_distance(x,origin=center1,refvec=refvec1))).T
            sOuterEdge = np.array(sorted(outerPoints.T,key=lambda x: _tools.clockwiseangle_and_distance(x,origin=center2,refvec=refvec2))).T
            
            sMinEdge = np.array(sorted(minEdge.T,key=np.linalg.norm)).T
            sMaxEdge = np.array(sorted(maxEdge.T,key=np.linalg.norm)).T
            
            XY = np.concatenate([sInnerEdge,sMinEdge,np.fliplr(sOuterEdge),np.fliplr(sMaxEdge)],axis=1)#np.concatenate([minEdge,rmaxEdge,np.fliplr(maxEdge),rminEdge],axis=1)
            
            edgePolygon.append(shapely.geometry.polygon.Polygon(XY.T))
            
            originalPoints = np.concatenate([_innerPoints.T,_outerPoints.T,_minEdge,_maxEdge],axis=0)
            
            pointsContained = np.sum([edgePolygon[-1].contains(PointS(originalPoints[I][0],originalPoints[I][1])) for I in range(originalPoints.shape[0])])
            if pointsContained!=originalPoints.shape[0]:
                inside = np.array([edgePolygon[-1].contains(PointS(originalPoints[I][0],originalPoints[I][1])) for I in range(originalPoints.shape[0])])
                outside = np.logical_not(inside)
                plt.figure()
                plt.scatter(x,y,c='k')
                plt.scatter(originalPoints[inside][:,0],originalPoints[inside][:,1],c='g')
                plt.scatter(originalPoints[outside][:,0],originalPoints[outside][:,1],c='r',zorder=100)
                plt.plot(np.array(edgePolygon[-1].boundary.coords)[:,0],np.array(edgePolygon[-1].boundary.coords)[:,1],c='r')
                plt.title(i)
                raise AttributeError('Error! {} points are outside the found shape with energy !'.format(np.sum(outside),0.5*(EBins[i]+EBins[i+1])))
        return edgePolygon,EBins


    def saveNXsqom(self,saveFileName):
        """Save converted file into an NXsqom.

        Args:

            - saveFileName (string): File name to be saved into.

        """

        if not self.__hasattr__('original_file'):
            raise AttributeError('Data file does not have link to the original file. This is needed to make a complete copy when creating nxs-files')
        if not self.type =='nxs':
            raise AttributeError('Only nxs typed files can be saved as nxs-files.')

        datafile = self.original_file
        Intensity = self.I # Dont swap axis as they are correct!
        Monitor = self.Monitor
        QX = self.qx
        QY = self.qy
        DeltaE = self.energy 
        binning = self.binning
        Normalization = self.Norm
        H = self.h
        K = self.k
        L = self.l

        if os.path.exists(saveFileName):
            warnings.warn('The file {} exists alread. Old file will be renamed to {}.'.format(saveFileName,saveFileName+'_old'))
            os.rename(saveFileName,saveFileName+'_old')
        fd = hdf.File(saveFileName,'w')
        fs = hdf.File(datafile.fileLocation,'r')
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
        
        rawdata = proc.create_dataset('rawdata',shape=(1,),dtype='S200',data=np.string_(os.path.realpath(datafile.fileLocation)))
        rawdata.attrs['NX_class']=b'NX_CHAR'

        normalizationString = proc.create_dataset('binning',shape=(1,),dtype='int32',data=binning)
        normalizationString.attrs['NX_class']=b'NX_INT'
        
        data = fd.get('entry/data')
        
        fileLength = Intensity.shape
        
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

        en = data.create_dataset('en',shape=(fileLength),dtype='float32',data=DeltaE)
        en.attrs['NX_class']=b'NX_FLOAT'

        h = data.create_dataset('h',shape=(fileLength),dtype='float32',data=H)
        k = data.create_dataset('k',shape=(fileLength),dtype='float32',data=K)
        l = data.create_dataset('l',shape=(fileLength),dtype='float32',data=L)
        for x in [h,k,l]:
            x.attrs['NX_class']=b'NX_FLOAT'

        fd.close()

            
def decodeStr(string):
    try:
        if 'decode' in string.__dir__():
            return string.decode('utf8')
        else:
            return string
    except:
        return string

@_tools.KwargChecker()
def getScanParameter(f):
    """Extract scan parameter from hdf file.

    Args:

        - f (hdf): Open HDF5 file object from which parameters are extracted.

    """
    if f.get('/entry/scanvars') is None:
        return [],[],[]
    scanParameters = [x.decode() for x in f.get('/entry/scanvars')]
    scanValues = []
    scanUnits = []

    if 'a3' in scanParameters:
        a3Id = scanParameters.index("a3")
        scanParameters[a3Id] = 'rotation_angle'
    else:
        a3Id = None

    if 'a4' in scanParameters:
        a4Id = scanParameters.index("a4")
        scanParameters[a4Id] = 'polar_angle'
    else:
        a4Id = None
    if 'ei' in scanParameters:
        eiId = scanParameters.index("ei")
        scanParameters[eiId] = 'incident_energy'
    else:
        eiId = None
    
        

    dataGroup = f.get('/entry/data')
    for d in dataGroup:
        if d in scanParameters:
            SCP = dataGroup[d]
            scanValues.append(np.array(SCP))
            if 'units' in list(SCP.attrs.keys()):
                scanUnits.append(decodeStr(SCP.attrs['units']))
            else:
                scanUnits.append('Unknown')
        
    if not a3Id is None:
        scanParameters[a3Id] = 'a3'
    if not a4Id is None:
        scanParameters[a4Id] = 'a4'
    if not eiId is None:
        scanParameters[eiId] = 'ei'
    
    scanParameters = np.array(scanParameters)
    scanValues = np.array(scanValues)
    scanUnits = np.array(scanUnits)


    return scanParameters,scanValues,scanUnits



@_tools.KwargChecker()
def createEmptyDataFile(A3,A4,Ei,sample,Monitor=50000, A3Off = 0.0, A4Off = 0.0,
                        title='EmptyDataFileTitle', name='EmptyDataFile',
                        temperature = None, electricField = None, magneticField = None,
                        detectors = 104, pixels = 1024, normalizationFiles = None):
    """Create an empty data file with a given sample and measurement parameters.

    Args:

        - A3 (list or float): Value(s) of measurement A3.

        - A4 (list or float): Value(s) of measurement A4.

        - Ei (list or float): Value(s) of measurement Ei.

        - sample (MJOLNIR Sample): Sample measured in data file.

    Kwargs:

        - Monitor (int): Monitor count for datafile (default 50 000)

        - A3Off (float): Offset in A3 used in datafile (default 0.0)

        - A4Off (float): Offset in A4 used in datafile (default 0.0)

        - title (string): Title of datafile (default EmptyDataFileTitle)

        - name (string): name of datafile (default EmptyDataFile)

        - temperature (double): Sample temperature (default None)
        
        - electricField (double): Sample electricField (default None)

        - magneticField (double): Sample magneticField (default None)

        - detectors (int): Number of detectors in spectrometer (default 104)

        - pixels (int): Number of pixels/detector in spectrometer (default 1024)

        - normalizationFiles (list or string): List or string to normalization file (default None)

    .. warning::
        If no normalization file(s) is/are provided, the resulting data file cannot be converted to HKL!

    """
    df = DataFile()
        
    A3 = np.asarray([A3]).flatten()
    A4 = np.asarray([A4]).flatten()
    Ei = np.asarray([Ei]).flatten()
    isChanging = np.array([len(A3)>1,len(A4)>1,len(Ei)>1])
    isChangingData = np.array([A3,A4,Ei])[isChanging]
    
    if np.sum(isChanging)>1:
            # Check if all arrays then have same shape
            if not np.all([x.shape == isChangingData[0].shape for x in isChangingData[1:]]):
                    names = np.array(['A3','A4','Ei'])
                    raise AttributeError('More than one parameter is changing but they do not have the same shape! Chaning: {}'.format(', '.join(str(x) for x in names[isChanging])))
    elif np.sum(isChanging)==0:
        raise AttributeError('No parameter is changing. At least one parameter must be chaning through the scan.')
    steps = len(isChangingData[0])
    
    Monitor = np.asarray([Monitor]*steps)
    
    df.A3Off = Marray(A3Off)
    df.A4Off = Marray(A4Off)
    df.Monitor = Marray(Monitor)
    df.MonitorPreset = Monitor[0]
    df.MonitorMode = 'm'



    df.sample = sample
    df.title = title
    df.name = name
    df.fileLocation = 'Unknown'
    
    df.temperature = temperature
    df.electricField = electricField
    df.magneticField = magneticField
    
    units = np.array(['degree','degree','meV'])
    params= np.array(['a3','a4','ei'])
    df.scanUnits = units[isChanging]
    df.scanParameters = params[isChanging]
    
    df.type = 'hdf'
    df.scanCommand = 'Unknown'
    df.scanValues = isChangingData
    df.instrument = 'CAMEA'
    
    df.Time = df.Monitor*0.0013011099243 # Typical measurement time RITA2 in 2018
    
    df.Ei = Marray(Ei)
    df.A3 = Marray(A3)
    df.A4 = Marray(A4)
    
    df.I = Marray(np.zeros((steps,detectors,pixels)))
    df.binning = None
    
    if not normalizationFiles is None:
            calib = []
            binning = []
            for f in normalizationFiles:
                    data = np.loadtxt(f,skiprows=3,delimiter=',')
                    
                    EfTable = data[:,[3,4,5,6]]
                    A4 = data[:,-1]
                    bound = data[:,[7,8]]
                    calib.append([EfTable,A4,bound])
                    binning.append(len(A4)/(104*8))
            df.instrumentCalibrations = np.array(calib)
            df.possibleBinnings = binning
            df.loadBinning(1)
    
    return df


class Sample(object):
    """Sample object to store all infortion of the sample from the experiment"""
    @_tools.KwargChecker()
    def __init__(self,a=1.0,b=1.0,c=1.0,alpha=90,beta=90,gamma=90,sample=None,name='Unknown',projectionVector1=None, projectionVector2 = None):
        if isinstance(sample,hdf._hl.group.Group):
            self.name = str(np.array(sample.get('name'))[0].decode())
            self.orientationMatrix = np.array(sample.get('orientation_matrix'))*2*np.pi

            self.planeNormal = np.array(sample.get('plane_normal'))
            
            self.polarAngle = np.array(sample.get('polar_angle'))
            self.rotationAngle = np.array(sample.get('rotation_angle'))
            self.unitCell = np.array(sample.get('unit_cell'))
            self.plane_vector1 = np.array(sample.get('plane_vector_1'))
            self.plane_vector2 = np.array(sample.get('plane_vector_2'))
            self.planeNormal = np.cross(self.plane_vector1[:3],self.plane_vector2[:3])
            self.A3Off = np.array([0.0])#
            if not np.isclose(np.linalg.norm(self.plane_vector1[:3].astype(float)),0.0) or not np.isclose(np.linalg.norm(self.plane_vector2[:3].astype(float)),0.0): # If vectors are not zero
                self.projectionVector1,self.projectionVector2 = calcProjectionVectors(self.plane_vector1.astype(float),self.plane_vector2.astype(float))
            else:
                self.projectionVector1,self.projectionVector2 = [np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0])]
            self.initialize()
            self.calculateProjections()

            
            
        elif np.all([a is not None,b is not None, c is not None]):
            self.unitCell = np.array([a,b,c,alpha,beta,gamma])
           
            self.polarAngle = np.array(0)
            self.rotationAngle = np.array(0)
            self.name=name
            if projectionVector1 is None or projectionVector2 is None:
                vector1,vector2 = [np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0])]
            else:
                vector1 = np.array(projectionVector1)
                vector2 = -np.array(projectionVector2)

            self.planeNormal = np.cross(vector1,vector2)
            
            self.planeNormal=self.planeNormal/np.linalg.norm(self.planeNormal)
            
            self.initialize()
            
            # Calcualte angles for plane vectors
            Ei = 5.0 
            k = np.sqrt(Ei)*factorsqrtEK
            H1,K1,L1 = vector1

            STT1 = TasUBlib.calTwoTheta(self.B,[H1,K1,L1,Ei,Ei],1.0)
            
            theta1 = TasUBlib.calcTheta(k,k,STT1)

            self.plane_vector1 = np.array([H1,K1,L1, theta1, STT1, 0.0,0.0, Ei,Ei])

            r2 = TasUBlib.addAuxReflection(self.cell,self.plane_vector1,vector2)
            STT2 = r2[4]
            theta2 = r2[3]
            self.plane_vector2 = r2
            self.orientationMatrix = TasUBlib.calcTasUBFromTwoReflections(self.cell,self.plane_vector1,self.plane_vector2)
            self.projectionVector1,self.projectionVector2 = calcProjectionVectors(self.plane_vector1.astype(float),self.plane_vector2.astype(float))#,self.planeNormal.astype(float))
            self.calculateProjections()
        else:
            print(sample)
            print(a,b,c,alpha,beta,gamma)
            raise AttributeError('Sample not understood')

            
    @property
    def unitCell(self):
        return self._unitCelll

    @unitCell.getter
    def unitCell(self):
        return np.array([self.a,self.b,self.c,self.alpha,self.beta,self.gamma])#self._unitCell

    @unitCell.setter
    def unitCell(self,unitCell):
        self._unitCell = unitCell
        self.a = unitCell[0]
        self.b = unitCell[1]
        self.c = unitCell[2]
        self.alpha = unitCell[3]
        self.beta  = unitCell[4]
        self.gamma = unitCell[5]
        
    @property
    def a(self):
        return self._a

    @a.getter
    def a(self):
        return self._a

    @a.setter
    def a(self,a):
        if a>0:
            self._a = a
        else:
            raise AttributeError('Negative or null given for lattice parameter a')

    @property
    def b(self):
        return self._b

    @b.getter
    def b(self):
        return self._b

    @b.setter
    def b(self,b):
        if b>0:
            self._b = b
        else:
            raise AttributeError('Negative or null given for lattice parameter b')

    @property
    def c(self):
        return self._c

    @c.getter
    def c(self):
        return self._c

    @c.setter
    def c(self,c):
        if c>0:
            self._c = c
        else:
            raise AttributeError('Negative or null given for lattice parameter c')


    @property
    def alpha(self):
        return self._alpha

    @alpha.getter
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self,alpha):
        if alpha>0 and alpha<180:
            self._alpha = alpha
        else:
            raise AttributeError('Negative,null or above 180 degrees given for lattice parameter alpha')

    @property
    def beta(self):
        return self._beta

    @beta.getter
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self,beta):
        if beta>0 and beta<180:
            self._beta = beta
        else:
            raise AttributeError('Negative,null or above 180 degrees given for lattice parameter beta')

    @property
    def gamma(self):
        return self._gamma

    @gamma.getter
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self,gamma):
        if gamma>0 and gamma<180:
            self._gamma = gamma
        else:
            raise AttributeError('Negative,null or above 180 degrees given for lattice parameter gamma')


    def __eq__(self,other):
        if not isinstance(other,type(self)):
            return False
        return np.all([self.name==other.name,np.all(self.unitCell==other.unitCell)])#,\
        #np.all(self.orientationMatrix==other.orientationMatrix)])

    def initialize(self):

        # From http://gisaxs.com/index.php/Unit_cell
        self.realVectorA = np.array([self.a,0,0])
        self.realVectorB = self.b*np.array([cosd(self.gamma),sind(self.gamma),0.0])#np.dot(np.array([self.b,0,0]),rotationMatrix(0,0,self.gamma))
        self.realVectorC = self.c*np.array([cosd(self.beta),(cosd(self.alpha)-cosd(self.beta)*cosd(self.gamma))/sind(self.gamma),
        np.sqrt(1-cosd(self.beta)**2-((cosd(self.alpha)-cosd(self.beta)*cosd(self.gamma))/sind(self.gamma))**2)])#np.dot(np.array([self.c,0,0]),rotationMatrix(0,self.beta,0))
        
        self.volume = np.abs(np.dot(self.realVectorA,np.cross(self.realVectorB,self.realVectorC)))
        self.reciprocalVectorA = 2*np.pi*np.cross(self.realVectorB,self.realVectorC)/self.volume
        self.reciprocalVectorB = 2*np.pi*np.cross(self.realVectorC,self.realVectorA)/self.volume
        self.reciprocalVectorC = 2*np.pi*np.cross(self.realVectorA,self.realVectorB)/self.volume
        ## Ensure that aStar is along the x-axis
        RotMatrix = _tools.rotate2X(self.reciprocalVectorA)
        #self.reciprocalVectorA=np.dot(RotMatrix,self.reciprocalVectorA)
        #self.reciprocalVectorB=np.dot(RotMatrix,self.reciprocalVectorB)
        #self.reciprocalVectorC=np.dot(RotMatrix,self.reciprocalVectorC)

        bv1,bv2,bv3 = self.reciprocalVectorA,self.reciprocalVectorB,self.reciprocalVectorC
        a1,a2,a3,alpha1,alpha2,alpha3= self.unitCell
        
        b1,b2,b3 = [np.linalg.norm(x) for x in [bv1,bv2,bv3]]
        beta1 = np.rad2deg(vectorAngle(bv2,bv3))
        beta2 = np.rad2deg(vectorAngle(bv3,bv1))
        beta3 = np.rad2deg(vectorAngle(bv1,bv2))
        self.cell = [a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]
        self.B = TasUBlib.calculateBMatrix(self.cell)
        #self.reciprocalMatrix = np.array([self.reciprocalVectorA,self.reciprocalVectorB,self.reciprocalVectorC]).T

    def calculateProjections(self):
        """Calculate projections and generate projection angles."""
        checks = np.array(['unitCell','orientationMatrix','projectionVector1','projectionVector2']) # 
        boolcheck = np.logical_not(np.array([hasattr(self,x) for x in checks]))
        if np.any(boolcheck):
            raise AttributeError('Sample object is missing: {}.'.format(', '.join(str(x) for x in checks[boolcheck])))
        
        V1 = self.projectionVector1.copy()
        #V1/=np.linalg.norm(V1)
        V2 = self.projectionVector2.copy()
        #V2/=np.linalg.norm(V2)
        pV1Q = np.dot(self.B,V1)
        pV2Q = np.dot(self.B,V2)
        self.projectionAngle = vectorAngle(pV1Q,pV2Q)

        if np.isclose(0.0,self.projectionAngle):
            raise AttributeError("The provided orientations are equal.")

        
        UB= self.orientationMatrix
        self.UB = UB

        self.orientationMatrixINV = np.linalg.inv(UB)
        

        p23 = np.array([[1,0,0],[0,1,0]]) # To extract Qx, Qy only
        PM = np.array([V1,V2]).T # Projection matrix
        self.PM = PM
        self.convert = np.dot(p23,np.einsum('ij,jk->ik',UB,PM)) # Convert from projX,projY to Qx, Qy
        self.convertHKL = np.dot(p23,UB) # Convert from HKL to Qx, Qy

        # Calculate 'misalignment' of the projection vector 1
        self.theta = -TasUBlib.calcTasMisalignment(UB,self.planeNormal,V1)
        
        self.RotMat = _tools.Rot(self.theta) # Create 2x2 rotation matrix
        
        self.convertinv = np.linalg.inv(self.convert) # Convert from Qx, Qy to projX, projY

        self.convertHKLINV = _tools.invert(self.convertHKL) # Convert from Qx, Qy to HKL

        # Calcualte RotationMatrix for UB as block diagonal matrix. Only Qx,Qy part rotates as'
        # calculated in RotMat
        self.RotMat3D = np.eye(3)
        self.RotMat3D[:2,:2] = self.RotMat
        #self.orientationMatrixINV = np.linalg.inv(np.dot(self.RotMat3D,UB))
        self.orientationMatrixINV = np.linalg.inv(self.UB)
        

    def tr(self,p0,p1):
        """Convert from projX, projY coordinate to Qx,QY coordinate."""
        p0, p1 = np.asarray(p0), np.asarray(p1)
        
        P = np.array([p0,p1])
        
        Pos = np.einsum('ij,j...->i...',self.convertinv,P)
        return Pos[0],Pos[1]
        
    def inv_tr(self, x,y):
        """Convert from Qx,QY  coordinate to projX, projY coordinate."""
        x, y = np.asarray(x), np.asarray(y)
        P = np.array([x,y])
        Pos = np.einsum('ij,j...->i...',self.convert,P)
        return Pos[0],Pos[1]   


    def format_coord(self,x,y): # Convert from Qx,Qy to HKL
        x, y = np.asarray(x), np.asarray(y)
        rlu = self.calculateQxQyToHKL(x,y)#np.dot(self.orientationMatrixINV,np.array([x,y,0]))
        return "h = {0:.3f}, k = {1:.3f}, l = {2:.3f}".format(rlu[0],rlu[1],rlu[2])

    def calculateQxQyToHKL(self,x,y): # convert from Qx,Qy to HKL
        pos = np.array([x,y,np.zeros_like(x)])
        return np.einsum('ij,j...->i...',self.orientationMatrixINV,pos)

    def calculateHKLToQxQy(self,H,K,L): # convert HKL to Qx,Qy
        pos = np.array([H,K,L])
        return np.einsum('ij,j...->i...',self.orientationMatrix,pos)[:2]

    def calculateHKLtoProjection(self,H,K,L):
        HKL = np.array([H,K,L])
        #points = np.einsum('i...,ij...->i...',HKL,self.PM)
        points = np.einsum('ij,j...->i...',_tools.invert(self.PM),HKL)
        return points


    def __str__(self):
        returnStr = 'Sample ' + self.name + '\n'
        #if not self.temperature is None: returnStr+= 'Temperatur: '+str(self.temperature)+'\n'
        #if not self.magneticField is None: returnStr+= 'Magnetic Field: '+str(self.magneticField)+'\n'
        #if not self.electricField is None: returnStr+= 'Electric Field: '+str(self.electricField)+'\n'
        returnStr+= 'Unit cell: \n' + str(self.unitCell) + '\n'
        returnStr+= 'Orientation matrix: \n' + str(self.orientationMatrix) +'\n'

        return returnStr

    @_tools.KwargChecker()
    def CurratAxe(self,Ei,Ef,Bragg,spurionType='Monochromator',HKL=False,Projection=False):
        """Function to calculate Currat-Axe position in QxQy coordinate system.
    
        Args:
            
            - sample (DataFile.Sample): Sample object for which Currat-Axe is to be calculated.
            
            - Ei (float/list): Incoming energy in meV.
            
            - Ef (float/list): Outgoing energy in meV.
            
            - Bragg (list): Bragg peak in HKL or list of.
            
        Kwargs:
            
            - spurionType (str): Either "Monochromator" or "Analyser" for origin of "wrong" energy (default "Monochromator").
            
            - HKL (bool): Whether or not to recalculate to HKL instead of Qx, Qy, Qz (default False).
            
            - Projection (Bool): Whether or not to recalculate to Projection vectors instead of Qx, Qy, Qz (default False).
            
        Returns:
            
            - Position (list): List of size (len(Bragg),len(Ei),len(Ef),3), where last axis is Qx, Qy, Qz
            
        """
        A3Off = self.theta
        
        Ei = np.asarray(Ei).flatten() # Shape (m)
        Ef = np.asarray(Ef).flatten() # Shape (n)
        Bragg = np.asarray(Bragg).reshape(-1,3) # shape (l,3)
        
        QlLocal = []
        if spurionType.lower() == 'monochromator':
            for B in Bragg:
                Ql = []
                Angles = np.array([TasUBlib.calcTasQAngles(self.orientationMatrix,self.planeNormal,1.0,A3Off,np.array([B[0],B[1],B[2],e,e]))[:2] for e in Ef])
                for ei in Ei:
                    Ql.append(np.array([TasUBlib.calcTasQH(self.orientationMatrixINV,angle,ei,e,0) for angle,e in zip(Angles,Ef)])[:,1])
                QlLocal.append(Ql)
        elif spurionType.lower() == 'analyser':
            for B in Bragg:
                Ql = []
                for ei in Ei:
                    Angles2 = np.array(TasUBlib.calcTasQAngles(self.orientationMatrix,self.planeNormal,1.0,A3Off,np.array([B[0],B[1],B[2],ei,ei]))[:2])
                    Ql.append(np.array([TasUBlib.calcTasQH(self.orientationMatrixINV,Angles2,ei,e,A3Off*0.0) for e in Ef])[:,1]) # Extract Qx,Qy
                QlLocal.append(Ql)
        else:
            raise AttributeError('Provided spurionType not understood. Expected "Monochromator" or "Analyser" but recieved "{}".'.format(spurionType))

        returnVal = np.array(QlLocal) # Shape (l,m,n,3)
        
        if HKL == True or Projection == True: # need to calculate HKL for Projection calculation
            returnValShape = np.array(returnVal.shape)
            returnVal = self.calculateQxQyToHKL(returnVal[:,:,:,0].flatten(),returnVal[:,:,:,1].flatten())
            if Projection == True:
                
                toProjection = self.calculateHKLtoProjection(returnVal[0],returnVal[1],returnVal[2])
                returnVal = np.array(self.inv_tr(toProjection[0],toProjection[1]))
                returnValShape[-1]=2 # reshape Qx,Qy,Qz dimension to P1,P2 (3 -> 2)
            
            returnVal.shape = returnValShape # Shape (l,m,n,3) or (l,m,n,2)
        return returnVal
    



@_tools.KwargChecker()
def rotationMatrix(alpha,beta,gamma,format='deg'):
    if format=='deg':
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    Rx = np.array([[1,0,0],[0,np.cos(alpha),-np.sin(alpha)],[0,np.sin(alpha),np.cos(alpha)]])
    Ry = np.array([[np.cos(beta),0,np.sin(beta)],[0,1,0],[-np.sin(beta),0,np.cos(beta)]])
    Rz = np.array([[np.cos(gamma),-np.sin(gamma),0],[np.sin(gamma),np.cos(gamma),0],[0,0,1]])
    return np.dot(Rz,np.dot(Ry,Rx))

def vectorAngle(V1,V2):
    return np.arccos(np.dot(V1,V2.T)/(np.linalg.norm(V1)*np.linalg.norm(V2)))

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

def extractData(files):
    if not isinstance(files,list):
        files = [files]
    I = []
    
    
    Norm = []
    Monitor = []
    a3 = []
    a4 = []
    a3Off = []
    a4Off = []
    Ei = []
    instrumentCalibrationEf = []
    instrumentCalibrationA4 = []
    instrumentCalibrationEdges = []

    if(files[0].type!='hdf'):
        qx = []
        qy = []
        energy = []
        H = []
        K = []
        L = []
    scanParameters = []
    scanParamValue = []
    scanParamUnit = []
    for datafile in files:
        I.append(datafile.I)
        if(files[0].type!='hdf'):
            qx.append(datafile.qx)
            qy.append(datafile.qy)
            energy.append(datafile.energy)
            Norm.append(datafile.Norm)
            H.append(datafile.h)
            K.append(datafile.k)
            L.append(datafile.l)
        scanParameters.append(datafile.scanParameters)
        scanParamValue.append(datafile.scanValues)
        scanParamUnit.append(datafile.scanUnits)
            
        
        Monitor.append(datafile.Monitor)
        if np.array(datafile.A3Off).shape is ():
            datafile.A3Off = 0.0
        a3.append(datafile.A3-datafile.A3Off)
        a3Off.append(datafile.A3Off)
        if np.array(datafile.A4Off).shape is ():
            datafile.A4Off = [0.0]
        a4.append(datafile.A4-datafile.A4Off)
        a4Off.append(datafile.A3Off)
        Ei.append(datafile.Ei)
        instrumentCalibrationEf.append(datafile.instrumentCalibrationEf)
        instrumentCalibrationA4.append(datafile.instrumentCalibrationA4)
        instrumentCalibrationEdges.append(datafile.instrumentCalibrationEdges)
        
    I = Marray(np.concatenate(I,axis=0))
    if(files[0].type!='hdf'):
        qx = Marray(np.concatenate(qx,axis=0))
        qy = Marray(np.concatenate(qy,axis=0))
        H = Marray(np.concatenate(H,axis=0))
        K = Marray(np.concatenate(K,axis=0))
        L = Marray(np.concatenate(L,axis=0))
        energy = Marray(np.concatenate(energy,axis=0))
        Norm = Marray(np.concatenate(Norm,axis=0))
    else:  
        Norm = Marray(Norm)
        
    Monitor = Marray(np.concatenate(Monitor,axis=0))

    a3 = Marray(np.concatenate(a3,axis=0))
    a4 = Marray(np.concatenate(a4,axis=0))

    a3Off = Marray(a3Off)
    a4Off = Marray(a4Off)
    instrumentCalibrationEf = np.array(instrumentCalibrationEf)
    instrumentCalibrationA4 = np.array(instrumentCalibrationA4)
    instrumentCalibrationEdges = np.array(instrumentCalibrationEdges)
    Ei = Marray(Ei)
    if files[0].type!='hdf':
        return I,qx,qy,energy,Norm,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit,H,K,L
    else:
        return I,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit


''' def invert(M):
    """Invert non-square matrices as described on https://en.wikipedia.org/wiki/Generalized_inverse.
    
    Args:
        
        - M (matrix): Matrix in question.
        
    Returns:
        
        - Left or right inverse matrix depending on shape of provided matrix.
    """
    s = M.shape
    if s[0]>s[1]:
        return np.dot(np.linalg.inv(np.dot(M.T,M)),M.T)
    else:
        return np.dot(M.T,np.linalg.inv(np.dot(M,M.T))) '''



    
def calcProjectionVectors(R1,R2,norm=None):
    r1 = R1[:3]
    r2 = R2[:3]
    if not norm is None:
        NV = norm
    else:
        NV = np.cross(r1,r2)
    NV/= np.linalg.norm(NV)
    
    Zeros = np.isclose(NV,0.0)
    if np.sum(Zeros)==3:
        raise AttributeError('The two plane vectors are equivalen, {}, {}!'.format(r1,r2))
    
    
    if np.sum(Zeros) == 2 or np.sum(Zeros)==1: # Easy case where the two vectors are to be along the x, y, or z directions
        if Zeros[0] == True:
            V1 = np.array([1.0,0.0,0.0])
            V2 = np.cross(V1,NV)
        elif Zeros[1]:
            V1 = np.array([0.0,1.0,0.0])
            V2 = np.cross(NV,V1)
        else:
            V1 = np.array([0.0,0.0,1.0])
            V2 = np.cross(NV,V1)
    else: # The tricky case of all vectors having non-zero components.
        V1 = r1
        V2 = r2
            
    V1 = _tools.LengthOrder(V1)
    V2 = _tools.LengthOrder(V2)
    #for V in [V1,V2]: # Flip sign if needed
    #    maxArg = np.argmax(np.abs(V))
    #    V*=np.sign(V[maxArg])
    return V1,V2





# --------------------------- TESTS -------------------------

def test_sample_exceptions():
    try: # No parameters given
        s1 = Sample()
        assert False
    except:
        assert True

    try: # negative parameters given
        s1 = Sample(a=-1,b=1,c=1)
        assert False
    except:
        assert True

    try: # negative parameters given
        s1 = Sample(a=1,b=-1,c=1)
        assert False
    except:
        assert True
    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=-1)
        assert False
    except:
        assert True

    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=1,alpha=200)
        assert False
    except:
        assert True
    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=1,beta=-10)
        assert False
    except:
        assert True
    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=1,gamma=-10)
        assert False
    except:
        assert True

def test_Sample_conversions():
    df = DataFile('Data/camea2018n000178.hdf')
    sample = df.sample

    # Check that tr and inv_tr results in the same
    p0,p1 = 2.25,-1.23
    Qx,Qy = sample.tr(p0,p1)
    P0,P1 = sample.inv_tr(Qx,Qy)

    assert(np.all(np.isclose([p0,p1],[P0,P1])))

    # Check that format_coord prints something
    string = sample.format_coord(Qx,Qy) # format is h = ??, k = ??, l = ??
    hkl = [float(x.split('=')[-1]) for x in string.split(',')]
    print(string)
    print(hkl)

    QxQyFromSample = sample.calculateHKLToQxQy(*hkl)
    print(QxQyFromSample)
    assert(np.all(np.isclose(QxQyFromSample,[Qx,Qy],atol=1e-3)))


def test_Sample_CurratAxe():
    df = DataFile('Data/camea2018n000178.hdf')
    sample = df.sample

    Ei = [5,7,9]
    Ef = np.array([[4,4,4],[5,5,5]])
    Bragg = [[1,0,0],[0,1,1]]


    pos = sample.CurratAxe(Ei,Ef,Bragg)

    assert(pos.shape == (len(Bragg),np.array(Ei).size,Ef.size,3))

    assert(np.all(np.isclose(pos[:,:,:,2],0.0))) # All Qz are to be zero


    POS = np.array([[[[ 1.59337451e-01,  6.96316473e-01,  0.00000000e+00],
            [ 1.59337451e-01,  6.96316473e-01,  0.00000000e+00],
            [ 1.59337451e-01,  6.96316473e-01,  0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00]],

            [[ 4.35859084e-01,  7.63659414e-01,  0.00000000e+00],
            [ 4.35859084e-01,  7.63659414e-01,  0.00000000e+00],
            [ 4.35859084e-01,  7.63659414e-01,  0.00000000e+00],
            [ 2.78156829e-01,  7.17745441e-01,  0.00000000e+00],
            [ 2.78156829e-01,  7.17745441e-01,  0.00000000e+00],
            [ 2.78156829e-01,  7.17745441e-01,  0.00000000e+00]],

            [[ 6.74964311e-01,  8.21890116e-01,  0.00000000e+00],
            [ 6.74964311e-01,  8.21890116e-01,  0.00000000e+00],
            [ 6.74964311e-01,  8.21890116e-01,  0.00000000e+00],
            [ 5.18676001e-01,  7.69828565e-01,  0.00000000e+00],
            [ 5.18676001e-01,  7.69828565e-01,  0.00000000e+00],
            [ 5.18676001e-01,  7.69828565e-01,  0.00000000e+00]]],


        [[[ 1.11373538e+00,  4.08688259e-01,  0.00000000e+00],
            [ 1.11373538e+00,  4.08688259e-01,  0.00000000e+00],
            [ 1.11373538e+00,  4.08688259e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00]],

            [[ 1.33491786e+00,  2.29586062e-01,  0.00000000e+00],
            [ 1.33491786e+00,  2.29586062e-01,  0.00000000e+00],
            [ 1.33491786e+00,  2.29586062e-01,  0.00000000e+00],
            [ 1.19906946e+00,  3.22887295e-01,  0.00000000e+00],
            [ 1.19906946e+00,  3.22887295e-01,  0.00000000e+00],
            [ 1.19906946e+00,  3.22887295e-01,  0.00000000e+00]],

            [[ 1.52617192e+00,  7.47183562e-02,  0.00000000e+00],
            [ 1.52617192e+00,  7.47183562e-02,  0.00000000e+00],
            [ 1.52617192e+00,  7.47183562e-02,  0.00000000e+00],
            [ 1.38306144e+00,  1.59458179e-01,  0.00000000e+00],
            [ 1.38306144e+00,  1.59458179e-01,  0.00000000e+00],
            [ 1.38306144e+00,  1.59458179e-01,  0.00000000e+00]]]])

    assert(np.all(np.isclose(pos,POS)))

    POSAnalyser = np.array([[[[ 1.60279692e-01,  6.22804387e-01,  0.00000000e+00],
            [ 1.60279692e-01,  6.22804387e-01,  0.00000000e+00],
            [ 1.60279692e-01,  6.22804387e-01,  0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00]],

            [[ 4.41363759e-01,  5.77272258e-01,  0.00000000e+00],
            [ 4.41363759e-01,  5.77272258e-01,  0.00000000e+00],
            [ 4.41363759e-01,  5.77272258e-01,  0.00000000e+00],
            [ 2.80013947e-01,  6.06605614e-01,  0.00000000e+00],
            [ 2.80013947e-01,  6.06605614e-01,  0.00000000e+00],
            [ 2.80013947e-01,  6.06605614e-01,  0.00000000e+00]],

            [[ 6.85994178e-01,  5.47926749e-01,  0.00000000e+00],
            [ 6.85994178e-01,  5.47926749e-01,  0.00000000e+00],
            [ 6.85994178e-01,  5.47926749e-01,  0.00000000e+00],
            [ 5.24052917e-01,  5.73796337e-01,  0.00000000e+00],
            [ 5.24052917e-01,  5.73796337e-01,  0.00000000e+00],
            [ 5.24052917e-01,  5.73796337e-01,  0.00000000e+00]]],


        [[[ 1.00477106e+00,  3.48941284e-01,  0.00000000e+00],
            [ 1.00477106e+00,  3.48941284e-01,  0.00000000e+00],
            [ 1.00477106e+00,  3.48941284e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00]],

            [[ 1.06290681e+00,  6.98843240e-02,  0.00000000e+00],
            [ 1.06290681e+00,  6.98843240e-02,  0.00000000e+00],
            [ 1.06290681e+00,  6.98843240e-02,  0.00000000e+00],
            [ 1.03489627e+00,  2.31469030e-01,  0.00000000e+00],
            [ 1.03489627e+00,  2.31469030e-01,  0.00000000e+00],
            [ 1.03489627e+00,  2.31469030e-01,  0.00000000e+00]],

            [[ 1.13033943e+00, -1.67701473e-01,  0.00000000e+00],
            [ 1.13033943e+00, -1.67701473e-01,  0.00000000e+00],
            [ 1.13033943e+00, -1.67701473e-01,  0.00000000e+00],
            [ 1.09633291e+00, -7.27153902e-03,  0.00000000e+00],
            [ 1.09633291e+00, -7.27153902e-03,  0.00000000e+00],
            [ 1.09633291e+00, -7.27153902e-03,  0.00000000e+00]]]])


    posAna = sample.CurratAxe(Ei,Ef,Bragg,spurionType='ANAlYser')

    assert(np.all(np.isclose(posAna[:,:,:,2],0.0))) # All Qz are to be zero
    assert(np.all(np.isclose(posAna,POSAnalyser)))


    posMonoProjection = sample.CurratAxe(Ei,Ef,Bragg,spurionType='monoChrOmaTOR',Projection=True)
    posMonoHKL = sample.CurratAxe(Ei,Ef,Bragg,spurionType='monoChrOmaTOR',HKL=True)
    assert(posMonoProjection.shape == (len(Bragg),np.array(Ei).size,Ef.size,2))
    assert(posMonoHKL.shape == (len(Bragg),np.array(Ei).size,Ef.size,3))

def test_vectorAngle():
    v1 = np.array([1,0,0])
    v2 = np.array([0,1,0])
    v3 = np.array([1,1,1])
    theta1 = vectorAngle(v1,v2)
    theta2 = vectorAngle(v1,v3)
    theta3 = vectorAngle(v2,v3)
    print(theta2)
    assert(np.isclose(theta1,np.pi/2))
    assert(np.isclose(theta2,theta3))
    assert(np.isclose(theta2,0.955316618125))
    
def test_Rotation_matrix():
    M1 = rotationMatrix(0,0,0)

    assert(np.all(np.isclose(M1,np.identity(3))))

    M2 = rotationMatrix(90,0,0)
    M3 = rotationMatrix(np.pi/2,np.pi/3,np.pi/4,format='rad')
    M4 = rotationMatrix(90,60,45)
    
    assert(np.all(np.isclose(M2,np.array([[1,0,0],[0,0,-1],[0,1,0]]))))
    assert(np.all(np.isclose(M3,M4)))
    
    M4_check = np.array([[3.53553391e-01,6.12372436e-01,7.07106781e-01],[3.53553391e-01,6.12372436e-01,-7.07106781e-01],[-8.66025404e-01,5.00000000e-01,3.06161700e-17]])
    assert(np.all(np.isclose(M4,M4_check)))

def test_equality():
    s1 = Sample(1,2,3,90,90,120)
    s2 = Sample(1,2,3,90,90,120)
    assert(s1==s2)

def test_calculateProjections(): # TODO: Update test

    s1 = Sample(np.pi*2,np.pi*2,np.pi*2,90,90,60)

    s1.calculateProjections()

    theta = 2*np.pi/3 # 120 degrees between projection vectors
    AStar = np.linalg.norm(s1.reciprocalVectorA)
    BStar = np.linalg.norm(s1.reciprocalVectorB)
    CStar = np.linalg.norm(s1.reciprocalVectorC)
    ## Projection is given by [[a*,cos(theta)b*],[0,sin(tetha)b*]]
    #assert(np.all(np.isclose(s1.projectionMatrix,np.array([[AStar*1,BStar*np.cos(theta)],[0,BStar*np.sin(theta)]]))))
    #assert(np.all(np.isclose(s1.tr(1,1),np.array([AStar+np.cos(theta)*BStar,np.sin(theta)*BStar]))))

    #point = s1.tr(3.5,7.2)
    #reverse = s1.inv_tr(point[0],point[1])
    #assert(np.all(np.isclose(reverse,np.array([3.5,7.2]))))

    #string = s1.format_coord(point[0],point[1])
    #assert(string=='h = 3.500, k = 7.200, l = 0.000')


def test_DataFile():
    try:
        DF = DataFile('/nope.txt')
        assert False
    except:
        assert True

    try:
        DF= DataFile('Data/CAMEA_Full.xml') # Wrong file
        assert False
    except:
        assert True

    files = ['Data/camea2018n000137.hdf',
             'Data/camea2018n000137.nxs']
    DF1 = DataFile(files[0])
    assertFile(files[1])
    DF2 = DataFile(files[1])
    s = str(DF2)
    sampleS = str(DF2.sample)
    print(str(DF1.sample))
    print(str(DF2.sample))
    assert(DF1.sample == DF2.sample)

def test_DataFile_equility():
    f1 = DataFile('Data/camea2018n000136.hdf')
    print('----------')
    f2 = DataFile('Data/camea2018n000136.hdf')
    assert(f1==f2)
    print('----------')
    f3 = DataFile(f2)
    assert(f1==f3)
    print('----------')
    f3 = DataFile('Data/camea2018n000136.nxs')
    assert(f1==f3)
    print('----------')
    f4 = DataFile('Data/camea2018n000137.hdf')
    assert(f1!=f4)


def test_DataFile_plotA4():
    plt.ioff()
    fileName = 'Data/camea2018n000136.hdf'
    fileName2= 'Data/camea2018n000136.nxs'
    file = DataFile(fileName)
    

    try:
        file.plotA4(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotA4(1)
    fig2 = file.plotA4()
    assertFile(fileName2)
    file2 = DataFile(fileName2)
    try:
        file2.plotA4(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True
    
    file2.plotA4(binning=1)
    plt.close('all')

    
def test_DataFile_plotEf():
    plt.ioff()
    fileName = 'Data/camea2018n000136.hdf'
    fileName2= 'Data/camea2018n000136.nxs'
    assertFile(fileName2)
    file = DataFile(fileName)

    try:
        file.plotEf(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotEf(1)
    fig2 = file.plotEf()
    
    file2 = DataFile(fileName2)
    try:
        file2.plotEf(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    file2.plotEf(binning=1)
    plt.close('all')

def test_DataFile_plotEfOverview():
    plt.ioff()
    fileName = 'Data/camea2018n000136.hdf'
    fileName2= 'Data/camea2018n000136.nxs'
    assertFile(fileName2)

    file = DataFile(fileName)

    try:
        file.plotEfOverview(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotEfOverview(1)
    fig2 = file.plotEfOverview()

    file2 = DataFile(fileName2)
    try:
        file2.plotEfOverview(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    file2.plotEfOverview(binning=1)
    plt.close('all')

def test_DataFile_plotNormalization():
    plt.ioff()
    fileName = 'Data/camea2018n000136.hdf'
    fileName2= 'Data/camea2018n000136.nxs'
    file = DataFile(fileName)
    assertFile(fileName2)

    try:
        file.plotNormalization(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotNormalization(1)
    fig2 = file.plotNormalization()

    file2 = DataFile(fileName2)
    try:
        file2.plotNormalization(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    file2.plotNormalization(binning=1)
    plt.close('all')

def test_DataFile_decodeString():
    a = b'String'
    b = 'String'

    c =1.1 # Float

    assert(decodeStr(a)==decodeStr(b))
    assert(c == decodeStr(c))

def test_DataFile_ScanParameter():

    files = ['Data/camea2018n000136.hdf','Data/camea2018n000136.nxs']
    assertFile(files[1])
    for file in files:
        dfile = DataFile(file)
        assert(dfile.scanParameters[0]=='a3')
        assert(len(dfile.scanParameters)==len(dfile.scanUnits))
        assert(len(dfile.scanParameters)==len(dfile.scanValues))
        assert(len(dfile.scanParameters)==1)
        assert(dfile.scanUnits[0]=='degree')
        ##assert(np.all(dfile.scanValues==np.arange(0,150,1)))


def test_DataFile_Error():
    df = DataFile('Data/camea2018n000136.hdf')

    # Not implimented
    try:
        df+df
        assert False
    except NotImplementedError:
        assert True

    df.instrument = 'WrongInstrument'
    try:
        df.convert(binning=1)
        assert False
    except AttributeError:
        assert True
    
    df2 = DataFile('Data/camea2018n000136.nxs')
    try:
        df.saveNXsqom('Data/saving.nxs')
        assert False
    except AttributeError: # File does not have original_file attribute
        assert True

    df2.type = 'WrongType'
    try:
        df2.saveNXsqom('Data/saving.nxs')
        assert False
    except AttributeError: # Manually change type to wrong
        assert True

def test_DataFile_Sample_UB():
    df = DataFile('Data/camea2018n000136.hdf')
    s = df.sample
    b1,b2,b3 = [np.linalg.norm(x) for x in [s.reciprocalVectorA,s.reciprocalVectorB,s.reciprocalVectorC]]
    a1,a2,a3 = [np.linalg.norm(x) for x in [s.realVectorA,s.realVectorB,s.realVectorC]]
    beta1,beta2,beta3 = vectorAngle(s.reciprocalVectorB,s.reciprocalVectorC),vectorAngle(s.reciprocalVectorC,s.reciprocalVectorA),vectorAngle(s.reciprocalVectorA,s.reciprocalVectorB)
    alpha1,alpha2,alpha3 =  vectorAngle(s.realVectorB,s.realVectorC),vectorAngle(s.realVectorC,s.realVectorA),vectorAngle(s.realVectorA,s.realVectorB)

    cell = s.cell#[a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]
    r1 = s.plane_vector1
    r2 = s.plane_vector2

    UBCalc = TasUBlib.calcTasUBFromTwoReflections(cell,r1,r2)
    comparison = np.isclose(UBCalc,s.orientationMatrix,atol=1e-5) # Assert that they are equal to 1e-5 (Numerical error)
    print(np.round(UBCalc,5))
    print(np.round(s.orientationMatrix,5))
    assert(np.all(comparison))

def test_DataFile_CreateEmpty(): # TODO: Make this test!!!
    nf = np.array(['Data/Normalization_1.calib',
    'Data/Normalization_3.calib','Data/Normalization_8.calib'])

    A3 = np.linspace(0,180,181)
    A4 = -16
    Ei = 5.5
    Monitor = 1e5
    sample = Sample(a=6.0,b=6.0,c=12.2,projectionVector2=[1,0,0],projectionVector1=[0,2,1],gamma=120.,beta=80.,alpha=90.)

    try:
        _ = createEmptyDataFile(A3=10,A4=10,Ei=10,sample=sample) # No change in any parameter
        assert False
    except AttributeError:
        assert True
    
    try:
        _ = createEmptyDataFile(A3=[10,11],A4=[10,11,12],Ei=10,sample=sample) # Two parameters change but not with the same shape
        assert False
    except AttributeError:
        assert True
    

    df = createEmptyDataFile(A3=A3,A4=A4,Ei=Ei,sample=sample,Monitor=Monitor,normalizationFiles = nf)
    
    # Check the contents of df
    assert(df.sample == sample)
    assert(len(df.possibleBinnings)==len(nf))
    #assert(False)

def test_DataFile_Sample_Projection():
    df = DataFile('Data/camea2018n000136.hdf') # In A-B plane
    print(df.sample.projectionVector1,df.sample.projectionVector2)
    assert(np.all(np.isclose(df.sample.projectionVector1,np.array([1.0,0.0,0.0]))))
    assert(np.all(np.isclose(df.sample.projectionVector2,np.array([0.0,1.0,0.0]))))

    df2 = DataFile('Data/camea2018n000178.hdf') # in 1 1 0 and 0 0 1 plane
    assert(np.all(np.isclose(df2.sample.projectionVector1,np.array([0.0,0.0,1.0]))))
    assert(np.all(np.isclose(df2.sample.projectionVector2,np.array([1.0,1.0,0.0]))))


#
#def test_DataFile_BoundaryCalculation(quick):
#    if quick==True:
#        binning = [1,3,8]
#    else:
#        binning = [1]
#    for B in binning:
#        print('Using binning {}'.format(B))
#        df = DataFile('Data/camea2018n000017.hdf')
#        converted = df.convert(binning=B)
#        EP,EBins = converted.calculateEdgePolygons()
#        areas = np.array([e.area for e in EP])
#        assert(np.all(areas>2.0)) # Value found by running algorithm
#        assert(len(EBins)==B*8+1)
        



def assertFile(file):
    """Make sure that file exists for methods to work"""
    if not os.path.isfile(file):
        df = DataFile(file.replace('.nxs','.hdf'))
        con = df.convert(binning=8)
        con.saveNXsqom(file)