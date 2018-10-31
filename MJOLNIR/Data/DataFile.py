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

class DataFile(object):
    """Object to load and keep track of HdF files and their conversions"""
    def __init__(self,fileLocation):
        # Check if file exists
        if isinstance(fileLocation,DataFile): # Copy everything in provided file
            self.updateProperty(fileLocation.__dict__)
        else:
            if not os.path.isfile(fileLocation):
                raise AttributeError('File location does not exist({}).'.format(fileLocation))
            if fileLocation.split('.')[-1]=='nxs':
                self.type='nxs'
                with hdf.File(fileLocation) as f:
                    sample=f.get('/entry/sample')
                    self.sample = Sample(sample=f.get('/entry/sample'))
                    self.I=np.array(f.get('entry/data/intensity')).swapaxes(1,2)
                    self.qx=np.array(f.get('entry/data/qx'))
                    self.qy=np.array(f.get('entry/data/qy'))
                    self.h=np.array(f.get('entry/data/h'))
                    self.k=np.array(f.get('entry/data/k'))
                    self.l=np.array(f.get('entry/data/l'))
                    self.energy=np.array(f.get('entry/data/en'))
                    self.Norm=np.array(f.get('entry/data/normalization'))
                    self.Monitor=np.array(f.get('entry/data/monitor'))
                    instr = getInstrument(f)
                    self.instrument = instr.name.split('/')[-1]
                    self.possibleBinnings = np.array([int(x[-1]) for x in np.array(instr) if x[:5]=='calib'])
                    self.Ei = np.array(instr.get('monochromator/energy'))
                    self.A3 = np.array(f.get('entry/sample/rotation_angle')).reshape(-1)
                    self.A4 = np.array(f.get('entry/sample/polar_angle')).reshape(-1)
                    self.A3Off = np.array(f.get('entry/sample/rotation_angle_zero'))
                    self.A4Off = np.array(f.get('entry/sample/polar_angle_zero'))
                    self.binning = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/binning'))[0]
                    #self.instrumentCalibrationEf = np.array(f.get('entry/calibration/{}_pixels/ef'.format(str(self.binning))))
                    #self.instrumentCalibrationA4 = np.array(f.get('entry/calibration/{}_pixels/a4'.format(str(self.binning))))
                    #self.instrumentCalibrationEdges = np.array(f.get('entry/calibration/{}_pixels/edges'.format(str(self.binning))))
                    Ef = np.array(instr.get('calib{}/final_energy'.format(str(self.binning))))
                    width = np.array(instr.get('calib{}/width'.format(str(self.binning))))
                    bg = np.array(instr.get('calib{}/background'.format(str(self.binning))))
                    amp = np.array(instr.get('calib{}/amplitude'.format(str(self.binning))))
                    self.instrumentCalibrationEf = np.array([amp,Ef,width,bg]).T
                    self.instrumentCalibrationA4 = np.array(instr.get('calib{}/a4offset'.format(str(self.binning))))
                    self.instrumentCalibrationEdges = np.array(instr.get('calib{}/boundaries'.format(str(self.binning))))
                    self.temperature = np.array(sample.get('temperature'))
                    self.magneticField = np.array(sample.get('magnetic_field'))
                    self.electricField = np.array(sample.get('electric_field'))
                    self.scanParameters,self.scanValues,self.scanUnits = getScanParameter(f)
                    self.scanCommand = np.array(f.get('entry/scancommand'))
                    self.original_file = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/rawdata'))[0].decode()

            elif fileLocation.split('.')[-1]=='hdf':
                self.type='h5'
                with hdf.File(fileLocation) as f:
                    sample=f.get('/entry/sample')
                    self.sample = Sample(sample=f.get('/entry/sample'))
                    self.Monitor=np.array(f.get('entry/control/data'))
                    instr = getInstrument(f)
                    self.instrument = instr.name.split('/')[-1]
                    self.possibleBinnings = np.array([int(x[-1]) for x in np.array(instr) if x[:5]=='calib'])
                    self.Ei = np.array(instr.get('monochromator/energy'))
                    self.I = np.array(instr.get('detector/counts')).swapaxes(1,2)
                    self.A3 = np.array(f.get('entry/sample/rotation_angle'))
                    self.A4 = np.array(f.get('entry/sample/polar_angle')).reshape(-1)
                    self.A3Off = np.array(f.get('entry/sample/rotation_angle_zero'))
                    self.A4Off = np.array(f.get('entry/sample/polar_angle_zero'))
                    self.binning=1 # Choose standard binning 1
                    #self.instrumentCalibrationEf = np.array(f.get('entry/calibration/{}_pixels/ef'.format(str(self.binning))))
                    #self.instrumentCalibrationA4 = np.array(f.get('entry/calibration/{}_pixels/a4'.format(str(self.binning))))
                    #self.instrumentCalibrationEdges = np.array(f.get('entry/calibration/{}_pixels/edges'.format(str(self.binning))))
                    Ef = np.array(instr.get('calib{}/final_energy'.format(str(self.binning))))
                    width = np.array(instr.get('calib{}/width'.format(str(self.binning))))
                    bg = np.array(instr.get('calib{}/background'.format(str(self.binning))))
                    amp = np.array(instr.get('calib{}/amplitude'.format(str(self.binning))))
                    self.instrumentCalibrationEf = np.array([amp,Ef,width,bg]).T
                    self.instrumentCalibrationA4 = np.array(instr.get('calib{}/a4offset'.format(str(self.binning))))
                    self.instrumentCalibrationEdges = np.array(instr.get('calib{}/boundaries'.format(str(self.binning))))
                    self.temperature = np.array(sample.get('temperature'))
                    self.magneticField = np.array(sample.get('magnetic_field'))
                    self.electricField = np.array(sample.get('electric_field'))
                    self.scanParameters,self.scanValues,self.scanUnits = getScanParameter(f)
                    self.scanCommand = np.array(f.get('entry/scancommand'))
            else:
                raise AttributeError('File is not of type nxs or h5.')
            self.name = fileLocation.split('/')[-1]
            self.fileLocation = os.path.abspath(fileLocation)
            self.sample.calculateProjections()
            for key in ['magneticField','temperature','electricField']:
                if self.__dict__[key].dtype ==object: # Is np nan object
                    self.__dict__[key] = None
        

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

    def updateProperty(self,dictionary):
        if isinstance(dictionary,dict):
            for key in dictionary.keys():
                self.__setattr__(key,dictionary[key])

    def __eq__(self,other):
        return len(self.difference(other))==0
    
    def difference(self,other,keys = set(['sample','instrument','Ei','I','A3','A4','binning','scanParameters'])):
        """Return the difference between two data files by keys"""
        dif = []
        if not set(self.__dict__.keys()) == set(other.__dict__.keys()): # Check if same generation and type (h5 or nxs)
            return list(set(self.__dict__.keys())-set(other.__dict__.keys()))

        comparisonKeys = set(['sample','instrument','Ei','I','A3','A4','binning','scanParameters'])
        for key in comparisonKeys:
            skey = self.__dict__[key]
            okey = other.__dict__[key]
            if isinstance(skey,np.ndarray):
                try:
                    if not np.all(np.isclose(skey,okey)):
                        if not np.all(np.isnan(skey),np.isnan(okey)):
                            dif.append(key)
                except (TypeError, AttributeError):
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
            A4Zero = np.deg2rad(np.array(A4Zero))

        
        A3Zero = self.A3Off#file.get('entry/sample/rotation_angle_zero')
        if A3Zero is None:
            A3Zero=0.0
        else:
            A3Zero = np.deg2rad(np.array(A3Zero))

        A4 = np.deg2rad(A4Normalization)+A4Zero
        A4=A4.reshape(detectors,binning*EPrDetector,order='C')

        PixelEdge = EdgesNormalization.reshape(detectors,EPrDetector,binning,2).astype(int)
        factorsqrtEK = 0.694692
        
        A4File = self.A4#np.array(instrument.get('detector/polar_angle'))
        
        A4File = A4File.reshape((-1,1,1))

        A4Mean = A4.reshape((1,detectors,binning*EPrDetector))-np.deg2rad(A4File)
        
        Intensity=np.zeros((Data.shape[0],Data.shape[1],EPrDetector*binning),dtype=int)
        for i in range(detectors): # for each detector
            for j in range(EPrDetector):
                for k in range(binning):
                    Intensity[:,i,j*binning+k] = np.sum(Data[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)

        EfMean = EfNormalization[:,1].reshape(1,A4.shape[0],EPrDetector*binning)
        EfNormalization = (EfNormalization[:,0]*np.sqrt(2*np.pi)*EfNormalization[:,2]).reshape(1,A4.shape[0],EPrDetector*binning)

        

        kf = factorsqrtEK*np.sqrt(EfMean)#.reshape(1,detectors,binning*EPrDetector)
        Ei = self.Ei#np.array(instrument.get('monochromator/energy'))
        
        
        ki = factorsqrtEK*np.sqrt(Ei).reshape(-1,1,1)

        
        A3 = np.deg2rad(np.array(self.A3))+A3Zero #file.get('/entry/sample/rotation_angle/')
        
        if A3.shape[0]==1:
            A3 = A3*np.ones((steps))
        
        A3.resize((steps,1,1))
        
        # Shape everything into shape (steps,detectors,bins) (if external parameter is changed, this is assured by A3 reshape)
        Qx = ki-kf*np.cos(A4Mean)
        Qy = -kf*np.sin(A4Mean)
        QX = Qx*np.cos(A3)-Qy*np.sin(A3)
        QY = Qx*np.sin(A3)+Qy*np.cos(A3)
        DeltaE = Ei-EfMean
        if DeltaE.shape[0]==1:
            DeltaE = DeltaE*np.ones((steps,1,1))
        Monitor = self.Monitor.reshape((steps,1,1))#np.array(file.get('/entry/control/data'),dtype=int).reshape((steps,1,1))
        Monitor = Monitor*np.ones((1,detectors,EPrDetector*binning))
        Normalization = EfNormalization*np.ones((steps,1,1))

        shapes = QX.shape
        sample = self.sample
        pos = sample.inv_tr(QX.flatten(),QY.flatten())
        H,K,L = (sample.orientationMatrix[0].reshape(-1,1)*pos[0].reshape(1,-1)+sample.orientationMatrix[1].reshape(-1,1)*pos[1].reshape(1,-1)).reshape(3,*shapes)


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


        if binning is None or binning == self.binning:
            binning = self.binning
        else:
            with hdf.File(self.fileLocation) as f:
                # Check if binning is in file
                instr = getInstrument(f)
                
                binningsPossible = np.array([int(x[-1]) for x in np.array(instr) if x[:5]=='calib'])#np.array([int(x.split('_')[0]) for x in np.array(f.get('entry/calibration'))])
                if not binning in binningsPossible:
                    raise AttributeError('The provided binning ({}) is not present in the data file.'.format(binning))
                
                #self.instrumentCalibrationEf = np.array(f.get('entry/calibration/{}_pixels/ef'.format(str(self.binning))))
                #self.instrumentCalibrationA4 = np.array(f.get('entry/calibration/{}_pixels/a4'.format(str(self.binning))))
                #self.instrumentCalibrationEdges = np.array(f.get('entry/calibration/{}_pixels/edges'.format(str(self.binning))))
                Ef = np.array(instr.get('calib{}/final_energy'.format(str(binning))))
                width = np.array(instr.get('calib{}/width'.format(str(binning))))
                bg = np.array(instr.get('calib{}/background'.format(str(binning))))
                amp = np.array(instr.get('calib{}/amplitude'.format(str(binning))))
                self.instrumentCalibrationEf = np.array([amp,Ef,width,bg]).T
                self.instrumentCalibrationA4 = np.array(instr.get('calib{}/a4offset'.format(str(binning))))
                if self.instrumentCalibrationA4 is None:
                    print('self.instrumentCalibrationA4 is NONE with binning {}!!'.format(binning))
                self.instrumentCalibrationEdges = np.array(instr.get('calib{}/boundaries'.format(str(binning))))
                self.binning = binning 
        #return self

    @_tools.KwargChecker()
    def calculateEdgePolygons(self,addEdge=True):
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
                minT,maxT = minMax(TT)
                _minT,_maxT = minMax(TT)
                
                if addEdge:
                    mint, maxt = minMax(TT[np.logical_not(np.logical_or(np.isclose(TT,minT),np.isclose(TT,maxT)))])
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
            
                
            refvec1,center1 = calRefVector(innerPoints)
            refvec2,center2 = calRefVector(outerPoints)
            sInnerEdge = np.array(sorted(innerPoints.T,key=lambda x: clockwiseangle_and_distance(x,origin=center1,refvec=refvec1))).T
            sOuterEdge = np.array(sorted(outerPoints.T,key=lambda x: clockwiseangle_and_distance(x,origin=center2,refvec=refvec2))).T
            
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


    def saveNXsqom(self,savefilename):
        if not self.__hasattr__('original_file'):
            raise AttributeError('Data file does not have link to the original file. This is needed to make a complete copy when creating nxs-files')
        if not self.type =='nxs':
            raise AttributeError('Only nxs typed files can be saved as nxs-files.')

        datafile = self.original_file
        Intensity = self.I.swapaxes(1,2)
        Monitor = self.Monitor
        QX = self.qx
        QY = self.qy
        DeltaE = self.energy 
        binning = self.binning
        Normalization = self.Norm
        H = self.h
        K = self.k
        L = self.l

        if os.path.exists(savefilename):
            warnings.warn('The file {} exists alread. Old file will be renamed to {}.'.format(savefilename,savefilename+'_old'))
            os.rename(savefilename,savefilename+'_old')
        fd = hdf.File(savefilename,'w')
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
    scanParameters = str(np.array(f.get('/entry/scanvars'))).split()
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


class Sample(object):
    """Sample object to store all infortion of the sample from the experiment"""
    @_tools.KwargChecker()
    def __init__(self,a=None,b=None,c=None,alpha=90,beta=90,gamma=90,sample=None,name='Unknown'):
        if isinstance(sample,hdf._hl.group.Group):
            self.name = str(np.array(sample.get('name'))[0])
            self.orientationMatrix = np.array(sample.get('orientation_matrix'))
            self.planeNormal = np.array(sample.get('plane_normal'))
            self.polarAngle = np.array(sample.get('polar_angle'))
            self.rotationAngle = np.array(sample.get('rotation_angle'))
            self.unitCell = np.array(sample.get('unit_cell'))
            
        elif np.all([a is not None,b is not None, c is not None]):
            self.unitCell = np.array([a,b,c,alpha,beta,gamma])
            self.orientationMatrix = np.array([[1,0,0],[0,1,0]])
            self.planeNormal = np.array([0,0,1])
            self.polarAngle = np.array(0)
            self.rotationAngle = np.array(0)
            self.name=name
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

    def calculateProjections(self):
        """Calculate projections and generate projection angles."""
        try:
            self.unitCell
            self.orientationMatrix
        except:
            raise AttributeError('No unit cell is set for sample object.')
        self.realVectorA = np.array([self.a,0,0])
        self.realVectorB = np.dot(np.array([self.b,0,0]),rotationMatrix(0,0,self.gamma))
        self.realVectorC = np.dot(np.array([self.c,0,0]),rotationMatrix(0,self.beta,0))
        
        self.volume = np.abs(np.dot(self.realVectorA,np.cross(self.realVectorB,self.realVectorC)))
        self.reciprocalVectorA = 2*np.pi*np.cross(self.realVectorB,self.realVectorC)/self.volume
        self.reciprocalVectorB = 2*np.pi*np.cross(self.realVectorC,self.realVectorA)/self.volume
        self.reciprocalVectorC = 2*np.pi*np.cross(self.realVectorA,self.realVectorB)/self.volume
        ## Ensure that aStar is along the x-axis
        RotMatrix = rotate2X(self.reciprocalVectorA)
        #angle = vectorAngle(self.reciprocalVectorA,np.array([1,0,0])) # TODO: make general!!!
        self.reciprocalVectorA=np.dot(self.reciprocalVectorA,RotMatrix)
        self.reciprocalVectorB=np.dot(self.reciprocalVectorB,RotMatrix)
        self.reciprocalVectorC=np.dot(self.reciprocalVectorC,RotMatrix)

        reciprocalMatrix = np.array([self.reciprocalVectorA,self.reciprocalVectorB,self.reciprocalVectorC])
        self.projectionVector1 = np.array(np.dot(self.orientationMatrix[0],reciprocalMatrix))
        self.projectionVector2 = np.array(np.dot(self.orientationMatrix[1],reciprocalMatrix))
        self.projectionAngle = vectorAngle(self.projectionVector1,self.projectionVector2)

        if np.isclose(0,self.projectionAngle):
            raise AttributeError("The provided orientations are equal.")

        self.projectionMatrix = np.array([\
        [np.linalg.norm(self.projectionVector1),np.cos(self.projectionAngle)*np.linalg.norm(self.projectionVector2)]\
        ,[0,np.sin(self.projectionAngle)*np.linalg.norm(self.projectionVector2)]])

    def tr(self,h,k):
        """Convert from curved coordinate to rectlinear coordinate."""
        h, k = np.asarray(h), np.asarray(k)
        Pos = np.dot(self.projectionMatrix,[h,k])
        return Pos[0],Pos[1]
        
    def inv_tr(self, x,y):
        """Convert from rectlinear coordinate to curved coordinate."""
        x, y = np.asarray(x), np.asarray(y)
        Pos = np.dot(np.linalg.inv(self.projectionMatrix),[x,y])
        return Pos[0],Pos[1]   


    def format_coord(self,x,y):
        pos = self.inv_tr(x,y)
        rlu = self.orientationMatrix[0]*pos[0]+self.orientationMatrix[1]*pos[1]
        return "h = {0:.3f}, k = {1:.3f}, l = {2:.3f}".format(rlu[0],rlu[1],rlu[2])

    def __str__(self):
        returnStr = 'Sample ' + self.name + '\n'
        #if not self.temperature is None: returnStr+= 'Temperatur: '+str(self.temperature)+'\n'
        #if not self.magneticField is None: returnStr+= 'Magnetic Field: '+str(self.magneticField)+'\n'
        #if not self.electricField is None: returnStr+= 'Electric Field: '+str(self.electricField)+'\n'
        returnStr+= 'Unit cell: \n' + str(self.unitCell) + '\n'
        returnStr+= 'Orientation matrix: \n' + str(self.orientationMatrix) +'\n'

        return returnStr

    



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

    if(files[0].type!='h5'):
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
        if(files[0].type!='h5'):
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
        
    I = np.array(I)
    if(files[0].type!='h5'):
        qx = np.array(qx)
        qy = np.array(qy)
        H = np.array(H)
        K = np.array(K)
        L = np.array(L)
        energy = np.array(energy)

    scanParameters = np.array(scanParameters)
    scanParamValue = np.array(scanParamValue)
    scanParamUnit = np.array(scanParamUnit)

    Norm = np.array(Norm)
    Monitor = np.array(Monitor)

    a3 = np.array(a3)
    a4 = np.array(a4)

    a3Off = np.array(a3Off)
    a4Off = np.array(a4Off)
    instrumentCalibrationEf = np.array(instrumentCalibrationEf)
    instrumentCalibrationA4 = np.array(instrumentCalibrationA4)
    instrumentCalibrationEdges = np.array(instrumentCalibrationEdges)
    Ei = np.array(Ei)
    if files[0].type!='h5':
        return I,qx,qy,energy,Norm,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit,H,K,L
    else:
        return I,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit
def rotMatrix(v,theta): # https://en.wikipedia.org/wiki/Rotation_matrix
    v/=np.linalg.norm(v)
    m11 = np.cos(theta)+v[0]**2*(1-np.cos(theta))
    m12 = v[0]*v[1]*(1-np.cos(theta))-v[2]*np.sin(theta)
    m13 = v[0]*v[2]*(1-np.cos(theta))+v[1]*np.sin(theta)
    m21 = v[0]*v[1]*(1-np.cos(theta))+v[2]*np.sin(theta)
    m22 = np.cos(theta)+v[1]**2*(1-np.cos(theta))
    m23 = v[1]*v[2]*(1-np.cos(theta))-v[0]*np.sin(theta)
    m31 = v[0]*v[2]*(1-np.cos(theta))-v[1]*np.sin(theta)
    m32 = v[1]*v[2]*(1-np.cos(theta))+v[0]*np.sin(theta)
    m33 = np.cos(theta)+v[2]**2*(1-np.cos(theta))
    return np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])


def rotate2X(v):
    if np.isclose(v[2]/np.linalg.norm(v),1): # v is along z
        return rotMatrix([0,1,0],np.pi/2)
    # Find axis perp to v and proj v into x-y plane -> rotate 2 plane and then to x
    vRotInPlane = np.array([-v[1],v[0],0])
    vPlan = np.array([v[0],v[1],0])
    ## TODO: Check this!
    #if np.isclose(np.dot(v,vPlan)/(np.linalg.norm(v)*np.linalg.norm(vPlan)),1.0):
    #    return rotMatrix([1,0,0],0.0)
    theta = np.arccos(np.dot(v,vPlan)/(np.linalg.norm(v)*np.linalg.norm(vPlan)))
    R = rotMatrix(vRotInPlane,theta)
    v2 = np.dot(R,v)
    theta2 = np.arccos(np.dot(v2,np.array([1,0,0]))/np.linalg.norm(v2))
    R2 = rotMatrix(np.array([0,0,1.0]),-theta2)
    
    Rotation = np.dot(R2,R)
    return Rotation



def clockwiseangle_and_distance(point,origin=[0,0],refvec = [0,1]):
    """Sort points clockwise. Taken from https://stackoverflow.com/questions/41855695/sorting-list-of-two-dimensional-coordinates-by-clockwise-angle-using-python
    
    Args:
        
        - point (list): List of points in 2D of size 2xN
        
    Kwargs:
        
        - origin (list): Location of origin from which the points are to be sorted (default [0,0])
        
        - refvec (list): Vector direction for definition of zero point (default [0,1])
        
    """
    # Vector between point and the origin: v = p - o
    vector = [point[0]-origin[0], point[1]-origin[1]]
    # Length of vector: ||v||
    lenvector = math.hypot(vector[0], vector[1])
    # If length is zero there is no angle
    if lenvector == 0:
        return -math.pi, 0
    # Normalize vector: v/||v||
    normalized = [vector[0]/lenvector, vector[1]/lenvector]
    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    # Negative angles represent counter-clockwise angles so we need to subtract them 
    # from 2*pi (360 degrees)
    if angle < 0:
        return 2*math.pi+angle, lenvector
    # I return first the angle because that's the primary sorting criterium
    # but if two vectors have the same angle then the shorter distance should come first.
    return angle, lenvector

def calRefVector(points):
    """ Calcualte reference vector as vector pointing from mean point to geometric center. For half moon shape this is anti-radially.
    
    Args:
        
        - points (list): list of points for which reference vector is calcualted, shape is 2xN
        
    Returns:
        
        vector: Reference vector
    
    """
    center = np.mean(points,axis=1).reshape(2,1)
    argMinDist = np.argmin(np.linalg.norm(center-points,axis=0))
    return center-points[:,argMinDist].reshape(2,1),center


def minMax(x,axis=0):
    """Return minimal and maximal of list.
    
    Args:
        
        - x (list): Object from which min and max is to be found.
        
    Kwargs:
        
        - axis (int): Axis or axes along which to operate (default 0)
      
    Returns:
        
        - min: Minimal value
        
        - max: Maximal value
        
    """
    return np.min(x),np.max(x)


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

def test_calculateProjections():

    s1 = Sample(np.pi*2,np.pi*2,np.pi*2,90,90,60)

    s1.orientationMatrix = np.array([[1,0,0],[0,1,0]])
    s1.calculateProjections()

    theta = 2*np.pi/3 # 120 degrees between projection vectors
    AStar = np.linalg.norm(s1.reciprocalVectorA)
    BStar = np.linalg.norm(s1.reciprocalVectorB)
    CStar = np.linalg.norm(s1.reciprocalVectorC)
    ## Projection is given by [[a*,cos(theta)b*],[0,sin(tetha)b*]]
    assert(np.all(np.isclose(s1.projectionMatrix,np.array([[AStar*1,BStar*np.cos(theta)],[0,BStar*np.sin(theta)]]))))
    assert(np.all(np.isclose(s1.tr(1,1),np.array([AStar+np.cos(theta)*BStar,np.sin(theta)*BStar]))))

    point = s1.tr(3.5,7.2)
    reverse = s1.inv_tr(point[0],point[1])
    assert(np.all(np.isclose(reverse,np.array([3.5,7.2]))))

    string = s1.format_coord(point[0],point[1])
    assert(string=='h = 3.500, k = 7.200, l = 0.000')


def test_DataFile():
    try:
        DF = DataFile('/nope.txt')
        assert False
    except:
        assert True
    files = ['TestData/1024/Magnon_ComponentA3Scan.hdf',
             'TestData/1024/Magnon_ComponentA3Scan.nxs']
    DF1 = DataFile(files[0])
    assertFile(files[1])
    DF2 = DataFile(files[1])
    s = str(DF2)
    sampleS = str(DF2.sample)
    print(str(DF1.sample))
    print(str(DF2.sample))
    assert(DF1.sample == DF2.sample)

def test_DataFile_equility():
    f1 = DataFile('TestData/1024/Magnon_ComponentA3Scan.hdf')
    f2 = DataFile('TestData/1024/Magnon_ComponentA3Scan.hdf')
    assert(f1==f2)

    f3 = DataFile(f2)
    assert(f1==f3)

def test_DataFile_rotations():
    vectors = [np.array([0,0,3.0]),np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0]),np.random.rand(3)]
    rotations = [rotate2X(v) for v in vectors]
    rotVector = [np.dot(rotations[i],vectors[i]) for i in range(len(vectors))]
    for i in range(len(rotVector)):
        assert(np.isclose(np.linalg.norm(vectors[i]),np.linalg.norm(rotVector[i])))
        print(rotVector[i][0],np.linalg.norm(rotVector[i]))
        assert(np.isclose(rotVector[i][0],np.linalg.norm(rotVector[i])))

def test_DataFile_plotA4():
    plt.ioff()
    fileName = 'TestData/1024/Magnon_ComponentA3Scan.hdf'
    fileName2= 'TestData/1024/Magnon_ComponentA3Scan.nxs'
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
    fileName = 'TestData/1024/Magnon_ComponentA3Scan.hdf'
    fileName2= 'TestData/1024/Magnon_ComponentA3Scan.nxs'
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
    fileName = 'TestData/1024/Magnon_ComponentA3Scan.hdf'
    fileName2= 'TestData/1024/Magnon_ComponentA3Scan.nxs'
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
    fileName = 'TestData/1024/Magnon_ComponentA3Scan.hdf'
    fileName2= 'TestData/1024/Magnon_ComponentA3Scan.nxs'
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

    assert(decodeStr(a)==decodeStr(b))


def test_DataFile_ScanParameter():

    files = ['TestData/1024/Magnon_ComponentA3Scan.hdf','TestData/1024/Magnon_ComponentA3Scan.nxs']
    assertFile(files[1])
    for file in files:
        dfile = DataFile(file)
        assert(dfile.scanParameters[0]=='a3')
        assert(len(dfile.scanParameters)==len(dfile.scanUnits))
        assert(len(dfile.scanParameters)==len(dfile.scanValues))
        assert(len(dfile.scanParameters)==1)
        assert(dfile.scanUnits[0]=='degrees')
        assert(np.all(dfile.scanValues==np.arange(0,150,1)))

def test_DataFile_BoundaryCalculation(quick):
    if quick==True:
        binning = [1,3,8]
    else:
        binning = [1]
    for B in binning:
        print('Using binning {}'.format(B))
        df = DataFile('TestData/1024/Magnon_ComponentA3Scan.hdf')
        converted = df.convert(binning=B)
        EP,EBins = converted.calculateEdgePolygons()
        areas = np.array([e.area for e in EP])
        assert(np.all(areas>2.0)) # Value found by running algorithm
        assert(len(EBins)==B*8+1)
        



def assertFile(file):
    """Make sure that file exists for methods to work"""
    if not os.path.isfile(file):
        df = DataFile(file.replace('.nxs','.hdf'))
        con = df.convert(binning=8)
        con.saveNXsqom(file)