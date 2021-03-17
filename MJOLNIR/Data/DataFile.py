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
import MJOLNIR
import datetime
import math
#import shapely
# from shapely.geometry import Polygon as PolygonS, Point as PointS
from MJOLNIR import TasUBlibDEG as TasUBlib
from MJOLNIR._tools import Marray
from MJOLNIR.Data import Mask

import MJOLNIR.Data.Sample
import re
import copy
import platform

multiFLEXXDetectors = 31*5
reFloat = r'-?\d*\.\d*'
factorsqrtEK = 0.694692
supportedRawFormats = ['hdf','MultiFLEXX','FlatCone']
supportedConvertedFormats = ['nxs']

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
            elif os.path.splitext(fileLocation)[1]=='': # No extension
                self.type = 'MultiFLEXX'
            else:
                raise AttributeError('File is not of type nxs or hdf.')
            self.name = os.path.basename(fileLocation)
            self.fileLocation = os.path.abspath(fileLocation)		
            self._binning = 1
            self._mask = False
            self.absolutNormalized = False

            if not self.type == 'MultiFLEXX':
                with hdf.File(fileLocation,mode='r') as f:
                    sample=f.get('/entry/sample')
                    self.sample = MJOLNIR.Data.Sample.Sample(sample=f.get('/entry/sample'))
                    instr = getInstrument(f)
                    if self.type == 'hdf':
                        if np.shape(np.array(instr.get('detector/counts'))) == ():
                            raise AttributeError('Data File {} has no data in {}/detector/counts. The file might be empty.'.format(self.name,instr.name))
                        self.I = np.array(instr.get('detector/counts')).swapaxes(1,2)
                    else:
                        self.I=np.array(f.get('entry/data/intensity'))
                        self.counts = np.array(instr.get('detector/counts')).swapaxes(1,2)
                        self.qx=np.array(f.get('entry/data/qx'))
                        self.qy=np.array(f.get('entry/data/qy'))
                        self.h=np.array(f.get('entry/data/h'))
                        self.k=np.array(f.get('entry/data/k'))
                        self.l=np.array(f.get('entry/data/l'))
                        self.energy=np.array(f.get('entry/data/en'))
                        self.Norm=np.array(f.get('entry/data/normalization'))
                    self.MonitorMode = np.array(f.get('entry/control/mode'))[0].decode()
                    self.MonitorPreset=np.array(f.get('entry/control/preset'))                
                    if self.type == 'hdf':
                        self.Monitor = np.array(f.get('entry/control/data'))
                        if not self.MonitorMode == 't' and len(self.Monitor)>1: # If not counting on time and more than one point saved
                            if self.Monitor.flatten()[0]!=self.MonitorPreset: # For all data in 2018 with wrong monitor saved
                                self.Monitor = np.ones_like(self.Monitor)*self.MonitorPreset ### TODO: Make Mark save the correct monitor!!
                    else:
                        self.Monitor=np.array(f.get('entry/data/monitor'))
                    self.Time = np.array(f.get('entry/control/time'))
                    self.startTime = np.array(f.get('entry/start_time'))[0]
                    self.endTime = np.array(f.get('entry/end_time'))[0]
                    self.experimentIdentifier = np.array(f.get('entry/experiment_identifier'))[0]
                    self.comment = np.array(f.get('entry/comment'))[0]
                    self.proposalId = np.array(f.get('entry/proposal_id'))[0]
                    self.proposalTitle = np.array(f.get('entry/proposal_title'))[0]

                    self.localContactName = np.array(f.get('entry/local_contact/name'))[0]
                    
                    self.proposalUserName = np.array(f.get('entry/proposal_user/name'))[0]
                    self.proposalUserEmail = np.array(f.get('entry/proposal_user/email'))[0]

                    self.userName = np.array(f.get('entry/user/name'))[0]
                    self.userEmail = np.array(f.get('entry/user/email'))[0]
                    self.userAddress = np.array(f.get('entry/user/address'))[0]
                    self.userAffiliation = np.array(f.get('entry/user/affiliation'))[0]

                    # Monochromator

                    attributes = ['type','d_spacing','horizontal_curvature','vertical_curvature',
                        'horizontal_curvature_zero','vertical_curvature_zero',
                        'gm','gm_zero','tlm','tlm_zero','tum','tum_zero','rotation_angle','rotation_angle_zero']
                    
                    values = ['monochromator'+x for x in ['Type','DSpacing','HorizontalCurvature',
                            'VerticalCurvature','HorizontalCurvatureZero','VerticalCurvatureZero',
                            'GM','GMZero','TLM','TLMZero','TUM','TUMZero','RotationAngle','RotationAngleZero']]
                    for att,value in zip(attributes,values):
                        loadedValue = f.get('entry/CAMEA/monochromator/{}'.format(att))
                        if not loadedValue is None: # if it does exist
                            setattr(self,value,np.array(loadedValue))

                    # MonochromatorSlit

                    attributes = [x+zero for x in ['bottom','left','right','top'] for zero in ['','_zero']]+\
                        ['x_gap','y_gap']
                    
                    values = ['monochromatorSlit'+x+zero for x in ['Bottom','Left','Right','Top'] for zero in ['','Zero']]+\
                        ['monochromatorSlitXGap','monochromatorSlitYGap']
                    for att,value in zip(attributes,values):
                        setattr(self,value,np.array(f.get('entry/CAMEA/monochromator_slit/{}'.format(att))))

                    # Analyzer
                    # analyzer_selection
                    attributes = ['d_spacing','nominal_energy','polar_angle','polar_angle_offset']
                    values = ['analyzer'+x.replace('_',' ').title().replace(' ','') for x in attributes]
                    for att,value in zip(attributes,values):
                        setattr(self,value,np.array(f.get('entry/CAMEA/analyzer/{}'.format(att))))
                    self.analyzerType = np.array(f.get('entry/CAMEA/analyzer/type'))[0]
                    self.analyzerSelection = int(np.array(f.get('entry/CAMEA/analyzer/analyzer_selection'))[0])
                    self.detectorSelection = int(np.array(f.get('entry/CAMEA/detector/detector_selection'))[0])

                    instr = getInstrument(f)
                    self.instrument = instr.name.split('/')[-1]
                    self.possibleBinnings = np.array([int(x[-1]) for x in np.array(instr) if x[:5]=='calib'])
                    self.Ei = np.array(instr.get('monochromator/energy'))
                    self.A3 = np.array(f.get('entry/sample/rotation_angle'))
                    self.A4 = np.array(instr.get('analyzer/polar_angle')).reshape(-1)
                    try:
                        self.A3Off = self.sample.A3Off#np.array(f.get('entry/sample/rotation_angle_zero'))  
                    except:
                        self.A3Off = [0.0]
                    self.A4Off = np.array(instr.get('analyzer/polar_angle_offset'))
                    if self.type == 'hdf':
                        self.binning=np.max(self.possibleBinnings).astype(int) # Choose standard binning max
                    else:
                        self.binning = int(np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/binning'))[0])
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
                    self.instrumentCalibrations = np.array(calibrations,dtype=object)
                    self.loadBinning(self.binning)
                        
                    self.temperature = np.array(sample.get('temperature'))
                    self.magneticField = np.array(sample.get('magnetic_field'))
                    self.electricField = np.array(sample.get('electric_field'))
                    self.scanParameters,self.scanValues,self.scanUnits,self.scanDataPosition = getScanParameter(f)
                    self.scanCommand = np.array(f.get('entry/scancommand'))[0]
                    if self.type == 'nxs':
                        self.original_file = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/rawdata'))[0].decode()
                    self.title = np.array(f.get('entry/title'))[0]

                    self.absoluteTime = np.array(f.get('entry/control/absolute_time'))
                    self.protonBeam = np.array(f.get('entry/proton_beam/data'))

                    try:
                        self.scanParameters,self.scanValues,self.scanUnits,self.scanDataPosition = getScanParameter(f)
                    except:
                        pass

                    if self.type == 'hdf':
                        ###################
                        #self.I[:,:,:150]=0#
                        ###################
                        pass
                    self.mask = np.zeros_like(self.I,dtype=bool)
                    if self.binning == 8:
                        self.mask[:,:,:2] = True
            else: # type is multiFLEXX
                self.loadMultiFLEXXData(fileLocation)
            
            self.scanSteps = self.scanValues.shape[1]
            if True:
                for attr in dir(self):
                    # If attribute is a function or property, skip it
                    if hasattr(getattr(self,attr),'__call__') or isinstance(getattr(self,attr),property):
                        continue
                    value = getattr(self,attr) 
                    if hasattr(value,'shape'): # Does it have a shape method -> Numpy array
                        if value.shape == (): # It is 0 D but either bytes or simple list
                            if isinstance(value,np.bytes_):
                                setattr(self,attr,value.decode()) # decode into text
                            else:
                                value = np.array([value]) # recast into 1D-array
                                setattr(self,attr,value) 

            for key in ['magneticField','temperature','electricField']:

                if hasattr(getattr(self,key),'dtype'):
                    if self.__dict__[key].dtype ==object: # Is np nan object
                        self.__dict__[key] = None
                else:
                    setattr(self,key,None)
        else: # Create an empty data set
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
            self._A3Off = np.array([0.0])
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
            self._A4Off = np.array([0.0])
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
        if not hasattr(A3,'shape'):
            self._A3 = np.array(A3)
        else:
            if A3.shape == ():
                self._A3 = np.array([0.0])
            else:
                if A3[0] is None:
                    self._A3 = np.array([0.0])
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
        if not hasattr(A4,'shape'):
            self._A4 = np.array(A4)
        else:
            if A4.shape == ():
                self._A4 = [0.0]
            else:
                self._A4 = A4

    @property
    def binning(self):
        return self._binning

    @binning.getter
    def binning(self):
        try:
            len(self._binning)
        except TypeError:
            return self._binning
        else:
            return self._binning[0]
        

    @binning.setter
    def binning(self,value):
        try:
            len(value)
        except TypeError:
            pass
        else:
            value = value[0]
        if hasattr(self,'possibleBinnings'):
            if not value in self.possibleBinnings:
                raise AttributeError('Wanted binning ({}) not allowed in {}. Possible binnings are: {}'.format(value,self.name,self.possibleBinnings))
            
            self._binning = value
        
            self.loadBinning(value)
        else:
            
            self._binning = value

    @property
    def mask(self):
        return self._mask

    @mask.getter
    def mask(self):
        return self._mask

    @mask.setter
    def mask(self,mask):
        if isinstance(mask,Mask.MaskingObject):
            coordsInside = np.array([hasattr(self,coord) for coord in mask.coordinates])
            if not np.all(coordsInside):
                names = np.array(mask.coordinates)
                raise AttributeError('Provided mask has coord(s) not in DataFile:','\n'.join(names[np.logical_not(coordsInside)]))
            self._maskingObject = mask
        elif hasattr(self,'_maskingObject'): # if _maskingObject exists from ealier, delete it
            del self._maskingObject
        if hasattr(self,'I'): # if identity has    
            if hasattr(self,'_maskingObject'):
                self._mask = self._maskingObject(self) # Generate boolean mask
            else:
                if not np.all(self.I.shape == mask.shape):
                    raise AttributeError('Shape of provided mask {} does not match shape of data {}.'.format(mask.shape,self.I.shape))
                self._mask = mask
        else:
            self._mask = mask

    def updateProperty(self,dictionary):
        if isinstance(dictionary,dict):
            for key in dictionary.keys():
                self.__setattr__(key,copy.deepcopy(dictionary[key]))

    def __eq__(self,other):
        return len(self.difference(other))==0
    
    def difference(self,other,keys = set(['sample','instrument','Ei','I','_A3','_A4','_binning','scanParameters','Monitor'])):
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
    def loadMultiFLEXXData(self,fileLocation,calibrationFile=None):
        """"Dedicated loader for MultiFLEXX data.

        Args:

            - fileLocation (str): Location of file

        Kwargs:

            - calibrationFile (str): Location of calibration file (default None: uses shipped calibration)

        """
        
        
        
        #if calibrationFile is None: # Use shipped calibration
        #    this_dir, _ = os.path.split(__file__)
        #    calibrationFile = os.path.realpath(os.path.join(this_dir,"..", "CalibrationMultiFLEXX.csv"))

        self.fileLocation = fileLocation
        self.possibleBinnings = [1] # Standard value (1 energy/detector)
        self.binning = 1

        with open(fileLocation) as f:
            dataString = f.readlines()
        dataString = ''.join(dataString)
        ## Load header data:
        headerStart = dataString.find('VVVVVVVVVVVVVVVVVVVV')
        headerEnd = dataString.find('DATA_:')
        headerString = dataString[headerStart:headerEnd].split('\n')[1:-1]
        dataPart = dataString[headerEnd:].split('\n')[1:-1]
        # Markers of multiple lines of parameters
        multiLines = ['param','varia','zeros','posqe']
        stringParameters = ['file_']
        parameters = {}


        for line in headerString:
            splitLine = line.split(': ')
            if isinstance(splitLine,list):
                description,value = splitLine
            else:
                continue
            description = description.lower()
            
            if not description in multiLines:
                ## Single value
                if not description in stringParameters:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                description=description.replace('_','')
                parameters[description] = value
                #print(description,value)
            elif description in ['varia','param','posqe']:
                
                keyValuePattern = re.compile(r'\w*\s*=\s*'+reFloat)
                KVPairs = keyValuePattern.findall(value)
                
                for pair in KVPairs:
                    key,val = pair.split('=')

                    key = key.strip()
                    
                    try:
                        val = np.array(val.strip(),dtype=float)
                    except ValueError:
                        continue
                    setattr(self,key,val)
            elif description == 'zeros':
                keyValuePattern = re.compile(r'\w+\s+=\s*'+reFloat)
                KVPairs = keyValuePattern.findall(value)
                
                for pair in KVPairs:
                    key,val = pair.split('=')
                    
                    key = key.strip()
                    try:
                        val = np.array(val.strip(),dtype=float)
                    except ValueError:
                        continue
                    key+='Off'
                    setattr(self,key,val)

        self.__dict__.update(parameters)

        ## Format scanCommand
        scanParameters = [x.split('=')[0].strip() for x in self.steps.split(',')]
        self.scanParameters = [p.capitalize() for p in scanParameters]

        self.scanUnits = []
        for param in self.scanParameters:
            if param.lower() in ['a3','a4','sgu','sgl']:
                unit= 'degree'
            elif param.lower() == 'ei':
                unit = 'meV'
            else:
                unit = 'unknown'
            self.scanUnits.append(unit)

        scanSteps = [float(x.split('=')[1]) for x in self.steps.split(',')]

        scanCommandNPMN = re.compile(r'\w+\s+\d+')
        NPMN = scanCommandNPMN.findall(self.comnd)[-2:]
        for paramVal in NPMN:
            param,val = paramVal.split()
            if param.lower() == 'np':
                self.steps = int(val)
            elif param.lower() == 'mn':
                self.MonitorPreset = int(val)
                self.MonitorMode = 'm'
            elif param.lower() == 'ti':
                self.MonitorPreset = int(val)
                self.MonitorMode = 't'

        def extractSample(obj,sample=None):
            if sample is None:
                sample = dict()
            
            nameConversion = [['AA','alpha'],
                            ['BB','beta'],
                            ['CC','gamma'],
                            ['AS','a'],
                            ['BS','b'],
                            ['CS','c'],
                            ]
            
            for oldName,newName in nameConversion:
                sample[newName] = getattr(obj,oldName)
                delattr(obj,oldName)
                
            Ei = 10.0 # TODO: What is the good solution here? Dummy incoming energy needed to calcualte UB
            k = np.sqrt(Ei)*factorsqrtEK
            q1 = [getattr(obj,x) for x in ['AX','AY','AZ']]
            q2 = [getattr(obj,x) for x in ['BX','BY','BZ']]
            cell = TasUBlib.calcCell([sample[x] for x in ['a','b','c','alpha','beta','gamma']])
            B = TasUBlib.calculateBMatrix(cell)
            A41 = TasUBlib.calTwoTheta(B,[*q1,Ei,Ei],-1)
            A31 = TasUBlib.calcTheta(k,k,A41)
            A42 = TasUBlib.calTwoTheta(B,[*q2,Ei,Ei],-1)
            A32 = TasUBlib.calcTheta(k,k,A42)

            planeVector1 = q1
            planeVector1.append(A31) # A3 
            planeVector1.append(A41) # A4
            [planeVector1.append(0.0) for _ in range(2)]# append values for gonios set to zero
            planeVector1.append(Ei)
            planeVector1.append(Ei)

            planeVector2 = q2
            planeVector2.append(A32) # A3 
            planeVector2.append(A42) # A4 
            [planeVector2.append(0.0) for _ in range(2)]# append values for gonios set to zero
            planeVector2.append(Ei)
            planeVector2.append(Ei)

            # add correct angle in theta between the two reflections
            between = TasUBlib.tasAngleBetweenReflections(B,np.array(planeVector1),np.array(planeVector2))

            planeVector2[3]+=between

            sample['projectionVector1']=np.array(planeVector1)
            sample['projectionVector2']=np.array(planeVector2)

            return sample


        def updateKeyName(obj,key,newName):
            if not hasattr(obj,key):
                return
            setattr(obj,newName,getattr(obj,key))
            delattr(obj,key)


        ## Find correct Ei

        KIPattern = re.compile(r'PARAM:\s*(FX=\s*\d*.?\d*),\s*(KFIX=\s*'+reFloat+')')


        EI = np.power(self.KFIX/0.694692,2)

        if np.isclose(self.FX,1.0): # If FX == 1 subtract, otherwise add
            self.EI = EI-self.EN
        else:
            self.EI = EI+self.EN



        sampleDict = extractSample(self)


        ## Find labels for data
        labels = dataPart[0].split()
        stepData = np.zeros((self.steps,len(labels)))
        singleDetectorData = np.zeros(self.steps,dtype = int)
        ## Extract the intensity data
        if dataPart[3].strip() == 'endflat': # Data is mix of steps and data in each line
            self.pntData = np.array([np.array(x.split(),dtype=float) for x in dataPart[1::3] if x[:8]!='Finished'])
            self.multiData = np.array([np.array(x.split()[1:],dtype=int) for x in dataPart[2::3]])
            self.multiData = self.multiData[:,:155] # Remove additional 5 zeros in data file
            self.type = 'MultiFLEXX'
            if 'Finished' in dataPart[-1]:
                self.endTime = dataPart[-1].replace('Finished at ','')
            else:
                self.endTime = 'N/A'
        else: # Assume flatcone, data is split in steps and then flatcone data
            # Split is marked with 'MULTI:'
            splitMarker = np.array([x=='MULTI:' for x in dataPart])
            split = splitMarker.argmax()
            if split !=0: # We have flatCone data!
                self.pntData = np.array([np.array(x.split(),dtype=float) for x in dataPart[1:split]])
                self.multiData = np.array([np.array(x.split(),dtype=int) for x in dataPart[split+1:]])
                self.type = 'FlatCone'
            else:
                self.pntData = np.array([np.array(x.split(),dtype=float) for x in dataPart[1:] if x[:8]!='Finished'])
                self.type='1D'
                

        ## Extract the data from above
        for parameter,value in zip(labels[1:],self.pntData[:,1:].T):
            setattr(self,parameter,np.array(value))

        dataLen = self.pntData.shape[0]
        if not dataLen == self.steps: # Scan stopped prematurely
            self.steps = dataLen

        if hasattr(self,'multiData'):
            updateKeyName(self,'multiData','I')
        else:
            updateKeyName(self,'CNTS','I')


        nameSwaps = [['file','name'],
                    ['MAG','magneticField'],
                    ['CNTS','ISingleDetector'],
                    ['M1','Monitor'],
                    ['TIME','Time'],
                    ['TT','temperature'],
                    ['comnd','scanCommand'],
                    ['instr','instrument'],
                    ['EI','Ei'],
                    ['local','localContact'],
                    ['user','userName'],
                    ['DA','analyzerDSpacing'],
                    ['expno','experimentIdentifier'],
                    ['localContact','localContactName'],
                    ['date','startTime']
                    ]

        for pair in nameSwaps:
            updateKeyName(self,pair[0],pair[1])


        self.scanValues = np.array([getattr(self,param) for param in self.scanParameters])

        ## Create sample from sample dictionary
        self.sample = MJOLNIR.Data.Sample.Sample(**sampleDict)
        self.sample.name = self.title
        self.comment = 'N/A'

        ## Create calibration data from file
        if calibrationFile is None: # Use shipped calibration
            if self.type in ['MultiFLEXX','FlatCone']:
                if self.type =='MultiFLEXX':
                    calibrationFile = MJOLNIR.__multiFLEXXNormalization__
                    detectors = 155
                    mask = np.zeros_like(self.I,dtype=bool)
                    self._mask = mask
                else:
                    calibrationFile = MJOLNIR.__flatConeNormalization__
                    detectors = 31
                    self.mask = False
                calibrationData = np.genfromtxt(calibrationFile,skip_header=1,delimiter=',')
                amplitude = calibrationData[:,3]
                background = calibrationData[:,6] 
                bound = calibrationData[:,7:9]
                final_energy = calibrationData[:,4]
                width = calibrationData[:,5]
                A4 = -calibrationData[:,9]
                
                EfTable = np.array([amplitude,final_energy,width,background]).T
                calibrations=[[EfTable,A4,bound]]
                bound =  np.array(detectors*[0,1],dtype=int).reshape(-1,2)
                
                self.instrumentCalibrationEdges = bound
                self.instrumentCalibrationEf = EfTable
                self.instrumentCalibrationA4 = A4
                
                self.instrumentCalibrations = calibrations
                if self.type == 'MultiFLEXX':
                    mask = np.zeros_like(self.I.data,dtype=bool)
                    mask[:,np.isnan(self.instrumentCalibrationEf[:,0])] = True
                    self.mask=mask
            elif self.type == '1D':
                pass
            else:
                import warnings
                warnings.warn('Instrument of type "{}" is not supported yet....'.format(self.type))
                #raise NotImplementedError('Instrument of type "{}" is not supported yet....'.format(self.type))

        # Create Dasel

        self.analyzerSelection = 0
        self.detectorSelection = 0
        ### Weird additional changes:
        if self.type == 'MultiFLEXX':
            self.A3Off = np.array([90.0])
            self.A4Off = np.array([0.0])
            


        if not hasattr(self,'electricField'):
            self.electricField = None
        if not hasattr(self,'magneticField'):
            self.magneticField = None
            

        
    @_tools.KwargChecker()
    def convert(self,binning=None):
        if self.instrument == 'CAMEA':
            EPrDetector = 8 
            if binning is None:
                binning = self.binning
        elif self.type in ['MultiFLEXX','FlatCone']:
            EPrDetector = 1
            if binning is None:
                binning = self.binning
        else:
            raise AttributeError('Instrument type of data file not understood. {} was given.'.format(self.instrument))
        self.loadBinning(binning)
        
        EfNormalization = self.instrumentCalibrationEf.copy()
        A4Normalization = self.instrumentCalibrationA4.copy()#np.array(instrument.get('calib{}/a4offset'.format(str(binning))))
        EdgesNormalization = self.instrumentCalibrationEdges.copy()#np.array(instrument.get('calib{}/boundaries'.format(str(binning))))
        Data = self.I.copy()#np.array(instrument.get('detector/data'))


        detectors = Data.shape[1]
        steps = Data.shape[0]
        
        if self.type in ['MultiFLEXX','FlatCone']:
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
        A4File = self.A4.copy()
        
        A4File = A4File.reshape((-1,1,1))

        A4Mean = (A4.reshape((1,detectors,binning*EPrDetector))+np.deg2rad(A4File-A4Zero))
        
        Intensity=np.zeros((Data.shape[0],Data.shape[1],EPrDetector*binning),dtype=int)
        for i in range(detectors): # for each detector
            for j in range(EPrDetector):
                for k in range(binning):
                    Intensity[:,i,j*binning+k] = np.sum(Data[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)

        
        EfMean = EfNormalization[:,1].reshape(1,A4.shape[0],EPrDetector*binning)
        EfNormalization = EfNormalization[:,0]#.reshape(1,A4.shape[0],EPrDetector*binning)#
        #EfNormalization = EfNormalization[:,0]*(np.sqrt(2*np.pi)*EfNormalization[:,2])
        
        EfNormalization.shape = (1,A4.shape[0],EPrDetector*binning)
        A3 = np.deg2rad(np.array(self.A3).copy())+A3Zero #file.get('/entry/sample/rotation_angle/')
        if A3.shape[0]==1:
            A3 = A3*np.ones((steps))
        
        A3.resize((steps,1,1))
        Ei = self.Ei.copy().reshape(-1,1,1)#np.array(instrument.get('monochromator/energy'))
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
            HKL,QX,QY = TasUBlib.calcTasQH(UBINV,[np.rad2deg(A3),
                np.rad2deg(A4Mean)],Ei,EfMean)
            H,K,L = np.swapaxes(np.swapaxes(HKL,1,2),0,3)
            self.sample.B = TasUBlib.calculateBMatrix(self.sample.cell)

        DeltaE = Ei-EfMean
        if DeltaE.shape[0]==1:
            DeltaE = DeltaE*np.ones((steps,1,1))
        Monitor = self.Monitor.copy().reshape((steps,1,1))
        Monitor = Monitor*np.ones((1,detectors,EPrDetector*binning))
        Normalization = EfNormalization*np.ones((steps,1,1))

       
        ###########################
        #Monitor[:,:,:binning] = 0 #
        ###########################


        convFile = DataFile(self) # Copy everything from old file
        updateDict = {'I':Intensity,'Monitor':Monitor,'qx':QX,'qy':QY,'energy':DeltaE,'binning':binning,'Norm':Normalization,
        'h':H,'k':K,'l':L,'type':'nxs','fileLocation':None,'original_file':self,'name':self.name.replace('.hdf','.nxs')}
        convFile.updateProperty(updateDict)

        if convFile.type == 'nxs' and convFile.binning == 8:
            convFile.mask = np.zeros_like(convFile.I,dtype=bool)
            convFile.mask[:,:,:2] = True


        if convFile.instrument == 'CAMEA':
            defectTubes = np.arange(104)[np.any(np.isnan(convFile.instrumentCalibrationA4.reshape(104,-1)),axis=1)] 

            if len(defectTubes)>0: # if any tubes are defect
                    if len(defectTubes)>1:
                        warnings.warn('Detector tubes {} masked'.format(', '.join([str(x) for x in defectTubes])))
                    else:
                        warnings.warn('Detector tube {} masked'.format(defectTubes[0]))

                    newMask = np.repeat(np.isnan(convFile.instrumentCalibrationA4.reshape(104,-1))[np.newaxis],len(convFile.I),axis=0)
                    if np.all(convFile.mask.shape==newMask.shape):
                        newMask = np.logical_or(convFile.mask,newMask)
                    convFile.mask = newMask

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
        for a4,N in zip(A4,Norm):
            plt.scatter(-a4,np.arange(len(a4)),c=N)
        
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


        if binning is None or not hasattr(self,'instrumentCalibrations'):
            return
        
        if not binning in self.possibleBinnings:
            raise AttributeError('Wanted binning not in possible binnings!')
        binID = list(self.possibleBinnings).index(binning)
        self.instrumentCalibrationEf,self.instrumentCalibrationA4,self.instrumentCalibrationEdges = self.instrumentCalibrations[binID]
        try:
            len(binning)
        except TypeError:
            pass
        else:
            binning = binning[0]
        self._binning = binning

        self.instrumentCalibrationEf.shape = (-1,4)
        self.instrumentCalibrationA4.shape = (-1)
        self.instrumentCalibrationEdges.shape = (-1,2)
        

    # @_tools.KwargChecker()
    # def calculateEdgePolygons(self,addEdge=True): # pragma: no cover
    #     """Method to calculate bounding polygon for all energies. The energies are split using the bin-edges method of DataSet. Hereafter,
    #     the outer most points are found in polar coordinates and a possible addition is made creating the padded bounding polygon.
        
    #     Kwargs:
            
    #         - addEdge (bool/float): If true, padding is found as difference between outer and next outer point. If addEdge is a number, generate padding a padding of this value (default True)
            
    #     Returns:
            
    #         - edgePolygon (list): List of shapely polygons of the boundary
            
    #         - EBins (list): Binning edges in energy
            
    #     """
    #     if addEdge:
    #         if np.isclose(float(addEdge),1.0):
    #             addEdgeAmount = None
    #         else:
    #             addEdgeAmount = addEdge
        
    #     Energy = self.energy
        
    #     E = np.sort(np.mean(self.energy,axis=(0,1)))
    #     dif = np.diff(E)*0.5
    #     EBins = np.concatenate([[E[0]-dif[0]],E[:-1]+dif,[E[-1]+dif[-1]]])
    #     steps = Energy.shape[0]
        
        
    #     edgePolygon = []
    #     for ELow,EHigh in zip(EBins,EBins[1:]):#range(len(EBins)-1):
            
    #         EBool = np.logical_and(Energy>ELow,Energy<EHigh)
            
    #         x = self.qx[EBool]
    #         y = self.qy[EBool]
            
            
    #         r = np.linalg.norm([x,y],axis=0)
    #         theta = np.arctan2(y,x)
            
    #         rBins = _tools.binEdges(r,0.00001)
            
    #         out = -1
    #         while np.sum(r>rBins[out])<steps:
    #             out-=1
    #         rOuter = r>rBins[out]
    #         inner = 0
    #         while np.sum(r<rBins[inner])<steps:
    #             inner+=1
    #         rInner = r<rBins[inner]
            
    #         minEdge = []
    #         maxEdge = []
    #         _minEdge= []
    #         _maxEdge= []
    #         include = 0
    #         for j in range(len(rBins)-1):
    #             Bool = np.logical_and(r>rBins[j-include],r<rBins[j+1])
    #             if np.sum(Bool)<steps-1:
    #                 include+=1
    #                 continue
    #             else:
    #                 include = 0
    #             TT = theta[Bool]
    #             dif = np.diff(TT)
    #             if np.max(abs(dif))>np.pi*1.9:
    #                 idx = np.argmax(abs(dif))
    #                 TT[idx+1:]+=-np.sign(dif[idx])*2*np.pi
    #             minT,maxT = _tools.minMax(TT)
    #             _minT,_maxT = _tools.minMax(TT)
                
    #             if addEdge:
    #                 mint, maxt = _tools.minMax(TT[np.logical_not(np.logical_or(np.isclose(TT,minT),np.isclose(TT,maxT)))])
    #                 if addEdgeAmount is None:
    #                     minT = minT-0.5*(mint-minT)
    #                     maxT = maxT-0.5*(maxt-maxT)
    #                 else:
    #                     minT = minT+np.sign(minT-mint)*addEdgeAmount
    #                     maxT = maxT+np.sign(maxT-maxt)*addEdgeAmount
                    
    #             R = np.mean(r[Bool])
    #             minEdge.append([R*np.cos(minT),R*np.sin(minT)])
    #             maxEdge.append([R*np.cos(maxT),R*np.sin(maxT)])
                
    #             _minEdge.append([R*np.cos(_minT),R*np.sin(_minT)])
    #             _maxEdge.append([R*np.cos(_maxT),R*np.sin(_maxT)])
            
    #         minEdge = np.array(minEdge).T
    #         maxEdge = np.array(maxEdge).T
            
    #         innerPoints = np.array([x[rInner],y[rInner]])
    #         _innerPoints= np.array([x[rInner],y[rInner]])
    #         if addEdge:
    #             if addEdgeAmount is None:
    #                 RR = rBins[inner]-(rBins[inner+1]-rBins[inner])
    #             else:
    #                 RR = rBins[inner]-addEdgeAmount
    #             Theta = np.arctan2(innerPoints[1],innerPoints[0])
    #             innerPoints = np.array([np.cos(Theta),np.sin(Theta)])*RR
                    
                
    #         outerPoints = np.array([x[rOuter],y[rOuter]])
    #         _outerPoints= np.array([x[rOuter],y[rOuter]])
    #         if addEdge:
    #             if addEdgeAmount is None:
    #                 RR = rBins[out]-(rBins[out-1]-rBins[out])
    #             else:
    #                 RR = rBins[out]+addEdgeAmount
    #             Theta = np.arctan2(outerPoints[1],outerPoints[0])
    #             outerPoints = np.array([np.cos(Theta),np.sin(Theta)])*RR
            
                
    #         refvec1,center1 = _tools.calRefVector(innerPoints)
    #         refvec2,center2 = _tools.calRefVector(outerPoints)
    #         sInnerEdge = np.array(sorted(innerPoints.T,key=lambda x: _tools.clockwiseangle_and_distance(x,origin=center1,refvec=refvec1))).T
    #         sOuterEdge = np.array(sorted(outerPoints.T,key=lambda x: _tools.clockwiseangle_and_distance(x,origin=center2,refvec=refvec2))).T
            
    #         sMinEdge = np.array(sorted(minEdge.T,key=np.linalg.norm)).T
    #         sMaxEdge = np.array(sorted(maxEdge.T,key=np.linalg.norm)).T
            
    #         XY = np.concatenate([sInnerEdge,sMinEdge,np.fliplr(sOuterEdge),np.fliplr(sMaxEdge)],axis=1)#np.concatenate([minEdge,rmaxEdge,np.fliplr(maxEdge),rminEdge],axis=1)
            
    #         edgePolygon.append(shapely.geometry.polygon.Polygon(XY.T))
            
    #         originalPoints = np.concatenate([_innerPoints.T,_outerPoints.T,_minEdge,_maxEdge],axis=0)
            
    #         pointsContained = np.sum([edgePolygon[-1].contains(PointS(originalPoints[I][0],originalPoints[I][1])) for I in range(originalPoints.shape[0])])
    #         if pointsContained!=originalPoints.shape[0]:
    #             inside = np.array([edgePolygon[-1].contains(PointS(originalPoints[I][0],originalPoints[I][1])) for I in range(originalPoints.shape[0])])
    #             outside = np.logical_not(inside)
    #             plt.figure()
    #             plt.scatter(x,y,c='k')
    #             plt.scatter(originalPoints[inside][:,0],originalPoints[inside][:,1],c='g')
    #             plt.scatter(originalPoints[outside][:,0],originalPoints[outside][:,1],c='r',zorder=100)
    #             plt.plot(np.array(edgePolygon[-1].boundary.coords)[:,0],np.array(edgePolygon[-1].boundary.coords)[:,1],c='r')
    #             plt.title(i)
    #             raise AttributeError('Error! {} points are outside the found shape with energy !'.format(np.sum(outside),0.5*(EBins[i]+EBins[i+1])))
    #     return edgePolygon,EBins


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
            if os.path.exists(saveFileName+'_old'):
                os.remove(saveFileName+'_old')
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
        qx.attrs['units']=b'1/angstrom'
        
        qy = data.create_dataset('qy',shape=(fileLength),dtype='float32',data=QY)
        qy.attrs['NX_class']=b'NX_FLOAT'
        qy.attrs['units']=b'1/angstrom'

        en = data.create_dataset('en',shape=(fileLength),dtype='float32',data=DeltaE)
        en.attrs['NX_class']=b'NX_FLOAT'
        en.attrs['units']=b'mev'

        h = data.create_dataset('h',shape=(fileLength),dtype='float32',data=H)
        k = data.create_dataset('k',shape=(fileLength),dtype='float32',data=K)
        l = data.create_dataset('l',shape=(fileLength),dtype='float32',data=L)
        for x in [h,k,l]:
            x.attrs['NX_class']=b'NX_FLOAT'
            x.attrs['units']=b'rlu'

        fd.close()

    def updateCalibration(self,calibrationFile,overwrite=False):
        """Update calibrations for the data file. Does not save the changes.
        
        Args:
            
            - calibrationFile (string or list): calibration file, as generated from MJOLNIR.Geometry.Instrument or list of these.
            
        Kwargs:
            
            - overwrite (bool): If true, previous binnings will be overwritten if new files contain same binning (default False)
        
    .. note::
        Changes performed by this method is not saved to disk. If this is wanted, use the saveNXsqom method.      
            
        """
        
        try:
            len(calibrationFile)
        except TypeError:
            calibrationFile = [calibrationFile]
            
        currentBinnings = self.possibleBinnings
        calibrations = {}
        newBinnings = []
        
        for binning in self.possibleBinnings:
            self.loadBinning(binning)
            calibrations[binning] = [self.instrumentCalibrationEf,self.instrumentCalibrationA4,self.instrumentCalibrationEdges]
        
        
        express = re.compile(r'\w*\spixel')
        
        for f in calibrationFile:
            with  open(f) as file:
                line = file.readline()
                
            pixel = int(express.findall(line)[0].split(' ')[0])
            if overwrite == False and pixel in currentBinnings: # Do not overwrote current values
                warnings.warn('Binning {} from file "{}" is skipped as overwrite is set to False.'.format(pixel,f))
                continue
                
            newBinnings.append(pixel)
            data = np.loadtxt(f,skiprows=3,delimiter=',')
            EfTable = data[:,[3,4,5,6]]
            A4 = data[:,-1]
            bound = data[:,[7,8]]
            calibrations[pixel] = [EfTable,A4,bound]

        self.instrumentCalibrations = np.array([c for c in calibrations.values()])
        self.possibleBinnings = np.array(list(calibrations.keys()))

    def updateSampleParameters(self,unitCell):
        """Update unit cell parameters and corresponding UB matrix

        Args:

            - unitCell (list): List of cell parameters (a,b,c,alpha,beta,gamma)

        """
        self.sample.updateSampleParameters(unitCell=unitCell) 
        

    def saveHDF(self,saveFileName):
        """Save current HDF file object into an HDF file.

        Args:

            - saveFileName (string): File name to be saved into.

        """
        
        
        def addMetaData(self,entry):
            dset = entry.create_dataset('start_time',(1,),dtype='<S70')
            dset[0] = np.string_(self.startTime)

            dset = entry.create_dataset('end_time',(1,),dtype='<S70')
            dset[0] = np.string_(self.endTime)
            
            dset = entry.create_dataset('experiment_identifier',(1,),dtype='<S70')
            dset[0] = self.experimentIdentifier.encode('utf8')

            dset = entry.create_dataset('instrument',(1,),dtype='<S70')
            dset[0] = self.instrument.title().upper().encode('utf8')

            dset = entry.create_dataset('comment',(1,),data=np.string_(self.comment))

            dset = entry.create_dataset('title',(1,),data=np.string_(self.title))

            dset = entry.create_dataset('proposal_id',(1,),data=np.string_(self.proposalId))

            dset = entry.create_dataset('proposal_title',(1,),data=np.string_(self.proposalTitle))

            cont = entry.create_group('local_contact')
            cont.attrs['NX_class'] = np.string_('NXuser')
            dset = cont.create_dataset('name',(1,),data=np.string_(self.localContactName))

            us = entry.create_group('proposal_user')
            us.attrs['NX_class'] = np.string_('NXuser')
            dset = us.create_dataset('name',(1,),data=np.string_(self.proposalUserName))
            dset = us.create_dataset('email',(1,),data=np.string_(self.proposalUserEmail))

            pus = entry.create_group('user')
            pus.attrs['NX_class'] = np.string_('NXuser')
            dset = pus.create_dataset('name',(1,),data=np.string_(self.userName))
            dset = pus.create_dataset('email',(1,),data=np.string_(self.userEmail))
            dset = pus.create_dataset('address',(1,),data=np.string_(self.userAddress))
            dset = pus.create_dataset('affiliation',(1,),data=np.string_(self.userAffiliation))

            

        def addMono(self,inst):
            mono = inst.create_group('monochromator')
            mono.attrs['NX_class'] = np.string_('NXmonochromator')
            
                
            dset = mono.create_dataset('type',(1,),dtype='S70')
            dset[0] = getattr(self,'monochromatorType')
            
            attributes = ['d_spacing','horizontal_curvature','vertical_curvature',
                'horizontal_curvature_zero','vertical_curvature_zero',
                'gm','gm_zero','tlm','tlm_zero','tum','tum_zero']
            units = ['angstrom']+['meter']*4+['degree']*6
            
            
            values = ['monochromator'+x for x in ['DSpacing','HorizontalCurvature',
                    'VerticalCurvature','HorizontalCurvatureZero','VerticalCurvatureZero',
                    'GM','GMZero','TLM','TLMZero','TUM','TUMZero']]
            
            for att,val,unit in zip(attributes,values,units):
                if val in self.__dict__:
                    dset = mono.create_dataset(att,(1,),'float32')
                    dset[0] = getattr(self,val)
                    dset.attrs['units'] = unit



            monoSlit = inst.create_group('monochromator_slit')
            monoSlit.attrs['NX_class'] = np.string_('NXmonochromatorslit')


            attributes = [x+zero for x in ['bottom','left','right','top'] for zero in ['','_zero']]
            values = ['monochromatorSlit'+x+zero for x in ['Bottom','Left','Right','Top'] for zero in ['','Zero']]
            
            for att,value in zip(attributes,values):
                val =  getattr(self,value)
                if not val.dtype == 'O':
                    dset = monoSlit.create_dataset(att,(1,),'float32')
                    dset[0] = val
                    dset.attrs['units'] = np.string_('mm')

            for value,att in zip(['monochromatorSlitXGap','monochromatorSlitYGap'],['x_gap','y_gap']):
                dset = monoSlit.create_dataset(att,(3,),'float32')
                dset[:] = getattr(self,value) # Might not work
                dset.attrs['units'] = 'mm'
        
        def addAna(self,inst):
            ana = inst.create_group('analyzer')
            ana.attrs['NX_class'] = np.string_('NXcrystal')
            
            attributes = ['d_spacing','nominal_energy','polar_angle','polar_angle_offset']
            values = ['analyzer'+x.replace('_',' ').title().replace(' ','') for x in attributes]
            units = ['anstrom','mev','degree','degree']
            for att,value,unit in zip(attributes,values,units):
                dset = ana.create_dataset(att,(1,),'float32')
                dset[0] = getattr(self,value)
                if not unit is None:
                    dset.attrs['units'] = np.string_(unit)
                
            dset = ana.create_dataset('type',data = np.array([np.string_(self.analyzerType)]))
            dset = ana.create_dataset('analyzer_selection',(1,),'int32',data=self.analyzerSelection)
            


        def addDetector(inst):
            det = inst.create_group('detector')
            det.attrs['NX_class'] = np.string_('NXdetector')

            
        def addSample(self,entry):
            sam = entry.create_group('sample')
            sam.attrs['NX_class'] = np.string_('NXsample')
            dset = sam.create_dataset('name',(1,),data=np.string_(self.sample.name))

            ub = self.sample.orientationMatrix/(2*np.pi) # 2pi is for change in convention
            
            dset = sam.create_dataset('orientation_matrix',data=ub)
            dset = sam.create_dataset('plane_vector_1',data=self.sample.plane_vector1)
            dset = sam.create_dataset('plane_vector_2',data=self.sample.plane_vector2)

            normal = self.sample.planeNormal
            dset = sam.create_dataset('plane_normal',data=normal)

            cell = np.array(self.sample.unitCell,dtype='float32')
            dset = sam.create_dataset('unit_cell',data=cell)

            dset = sam.create_dataset('azimuthal_angle',data=self.sample.azimuthalAngle)
            dset.attrs['units']=np.string_('degree')
            dset = sam.create_dataset('x',data=self.sample.x)
            dset.attrs['units']=np.string_('degree')
            dset = sam.create_dataset('y',data=self.sample.y)
            dset.attrs['units']=np.string_('degree')

            if hasattr(self,'temperature'):
                if not self.temperature is None:
                    dset = sam.create_dataset('temperature',data=self.temperature,dtype='float32')
                    dset.attrs['units'] = np.string_('K')

            if hasattr(self,'magneticField'):
                if not self.magneticField is None:
                    dset = sam.create_dataset('magnetic_field',data=self.magneticField,dtype='float32')
                    dset.attrs['units'] = np.string_('T')

            if hasattr(self,'electricField'):
                if not self.electricField is None:
                    dset = sam.create_dataset('electric_field',data=self.electricField,dtype='float32')
                    dset.attrs['units'] = np.string_('V') # TODO: Check if this unit is correct.

            for attr,value in zip(['sgu','sgl'],['sgu','sgl']):
                dset = sam.create_dataset(attr,(1,),data=getattr(self.sample,value))
                dset.attrs['units']=np.string_('degree')
                dset = sam.create_dataset(attr+'_zero',(1,),data=getattr(self.sample,value+'Zero'))
                dset.attrs['units']=np.string_('degree')
            
        def makeTheta(self):
            
            k = np.sqrt(self.Ei/2.072)
            fd = np.pi/(k*self.monochromatorDSpacing[0])
            theta = np.degrees(np.arcsin(fd))
            
            return theta,2*theta
        
            
        def storeScanData(self,entry):
            nxdata = entry.create_group('data')
            nxdata.attrs['NX_class'] = np.string_('NXdata')
            
            det = entry['CAMEA/detector']
            dset = det.create_dataset('counts',data=self.I.swapaxes(1,2), compression="gzip", compression_opts=6)
            dset.attrs['target'] = np.string_('/entry/CAMEA/detector/counts')
            nxdata['counts'] = dset
            
            dset = det.create_dataset('detector_selection',(1,),'int32',data=self.detectorSelection)
            
            dset = det.create_dataset('summed_counts',data=np.sum(self.I,axis=(1,2)))
            dset.attrs['target'] = np.string_('/entry/CAMEA/detector/summed_counts')
            nxdata['summed_counts'] = dset
            
            sam = entry['sample']

            dset = sam.create_dataset('rotation_angle',data=self.A3,dtype='float32')
            dset_zero = sam.create_dataset('rotation_angle_zero',data=self.A3Off,dtype='float32')

            dset.attrs['units'] = np.string_('degree')
            dset_zero.attrs['units'] = np.string_('degree')
            
            dset = sam.create_dataset('polar_angle',data=self.A4,dtype='float32')
            dset_zero = sam.create_dataset('polar_angle_zero',data=self.A4Off,dtype='float32')

            dset.attrs['units'] = np.string_('degree')
            dset_zero.attrs['units'] = np.string_('degree')
            dset.attrs['units'] = np.string_('degree')
            dset_zero.attrs['units'] = np.string_('degree')
            

            mono = entry['CAMEA/monochromator']
            
            dset = mono.create_dataset('energy',data=self.Ei,dtype='float32')
            dset.attrs['units'] = np.string_('mev')

            dset = mono.create_dataset('rotation_angle',data=self.monochromatorRotationAngle,dtype='float32')
            dset.attrs['units'] = np.string_('degree')
            dset = mono.create_dataset('rotation_angle_zero',data=self.monochromatorRotationAngleZero,dtype='float32')
            dset.attrs['units'] = np.string_('degree')


            entry.create_dataset('scancommand',(1,),data=np.string_(self.scanCommand))
            entry.create_dataset('scanvars',data=np.string_([x.encode('utf8') for x in self.scanParameters]))
            
            # save the correct scan variables 

            for variable,pos in zip(self.scanParameters,self.scanDataPosition):
                positionRelativeEntry = '/'.join([x for x in pos.split('/')[2:]])
                original = entry.get(positionRelativeEntry)
                nxdata[variable] = original
                nxdata[variable].attrs['target'] = np.string_('/entry/'+positionRelativeEntry)


            control = entry.create_group('control')
            control.attrs['NX_class'] = np.string_('NXmonitor')
            mons = self.Monitor
            control.create_dataset('data',data=mons,dtype='int32')
            dset = control.create_dataset('preset',(1,),dtype='int32')
            dset[0] = self.MonitorPreset
            dset = control.create_dataset('mode',(1,),data=np.string_(self.MonitorMode))
            time = self.Time
            dset = control.create_dataset('time',data=time,dtype='float32')
            dset.attrs['units'] = np.string_('seconds')

            time =  self.absoluteTime
            dset = control.create_dataset('absolute_time',data=time,dtype='float32')
            dset.attrs['units'] = np.string_('seconds')
            
            pb = entry.create_group('proton_beam')
            pb.attrs['NX_class'] = np.string_('NXmonitor')
            vals = self.protonBeam
            dset = pb.create_dataset('data',data=vals,dtype='int32')

        with hdf.File(saveFileName,'w') as f:
            
            f.attrs['file_name'] = np.string_(saveFileName)
            
            
            import datetime,time
            cT = datetime.datetime.now()
            
            f.attrs['file_time'] = np.string_('{}-{}-{}T{}:{}:{}{:+02.0f}:00'.format(cT.year,cT.month,cT.day,cT.hour,cT.minute,cT.second,-time.timezone/(60*60)))
            
            entry = f.create_group('entry')
            entry.attrs['NX_class'] = np.string_('NXentry')
        
            
            #------------ Instrument
            inst = entry.create_group(b'CAMEA')
            inst.attrs['NX_class'] = np.string_('NXinstrument')
            
            
        
            attribute = ['a4offset','amplitude','background','boundaries','final_energy','width']
            for calibration,binning in zip(self.instrumentCalibrations,self.possibleBinnings):
                pixelCalib = inst.create_group('calib{}'.format(binning))
                Etable, A4, bound = calibration
                amp,Ef,width,bg = Etable.T
                
                values = [A4,amp,bg,bound,Ef,width]
                dtypes = ['float32','float32','float32','int','float32','float32']
                units = ['degree',None,None,None,'mev','mev']
                for att,value,dtype,unit in zip(attribute,values,dtypes,units):
                    dset = pixelCalib.create_dataset(att,data=value,dtype=dtype)
                    if not unit is None:
                        dset.attrs['units']=np.string_(unit)
                
            
            addMetaData(self,entry)
            addMono(self,inst)
            addAna(self,inst)
            addDetector(inst)
            addSample(self,entry)
            storeScanData(self,entry)


            
def decodeStr(string):
    #try:
    if hasattr(string,'decode'):
        return string.decode('utf8')
    else:
        return string
    #except:
    #    return string

@_tools.KwargChecker()
def getScanParameter(f):
    """Extract scan parameter from hdf file.

    Args:

        - f (hdf): Open HDF5 file object from which parameters are extracted.

    """
    if f.get('/entry/data') is None:
        return [],[],[]
    
    scanParameters = []
    scanValues = []
    scanUnits = []
    scanDataPosition = []

    
    for item in f.get('/entry/data/'):
        if item in ['scanvar']: # only present in some files, but has to be ignored
            continue
        if not item in ['counts','summed_counts','en','h','intensity','k','l','monitor',
        'normalization','qx','qy'] and item[-4:]!='zero':
            scanParameters.append(item)
            fItem = f.get('/entry/data/{}'.format(item))
            scanUnits.append(decodeStr(fItem.attrs['units']))
            scanValues.append(np.array(fItem))
            try:
                scanDataPosition.append(decodeStr(fItem.attrs['target']))
            except:
                pass


    return scanParameters,np.array(scanValues),scanUnits,scanDataPosition



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
    isChangingData = np.array([A3,A4,Ei],dtype=object)[isChanging]
    
    if np.sum(isChanging)>1:
            # Check if all arrays then have same shape
            if not np.all([x.shape == isChangingData[0].shape for x in isChangingData[1:]]):
                    names = np.array(['A3','A4','Ei'])
                    raise AttributeError('More than one parameter is changing but they do not have the same shape! Changing: {}'.format(', '.join(str(x) for x in names[isChanging])))
    elif np.sum(isChanging)==0:
        raise AttributeError('No parameter is changing. At least one parameter must be changing through the scan.')
    steps = len(isChangingData[0])
    
    Monitor = np.asarray([Monitor]*steps)
    
    df.A3Off = np.array([A3Off])
    df.A4Off = np.array([A4Off])
    df.Monitor = np.array(Monitor)
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
    
    df.Ei = np.array(Ei)
    df.A3 = np.array(A3)
    df.A4 = np.array(A4)
    
    df.I = np.array(np.zeros((steps,detectors,pixels)))
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
            df.instrumentCalibrations = np.array(calib,dtype=object)
            df.possibleBinnings = binning
            df.loadBinning(1)
    
    
    df.mask = np.zeros_like(df.I,dtype=bool)
    return df



    





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
    if(files[0].type=='nxs'):
        qx = []
        qy = []
        energy = []
        H = []
        K = []
        L = []
    scanParameters = []
    scanParamValue = []
    scanParamUnit = []
    mask = []
    for datafile in files:
        I.append(datafile.I)
        if(datafile.type=='nxs'):
            qx.append(datafile.qx)
            qy.append(datafile.qy)
            energy.append(datafile.energy)
            Norm.append(datafile.Norm)
            H.append(datafile.h)
            K.append(datafile.k)
            L.append(datafile.l)
        mask.append(datafile.mask)
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
        a4Off.append(datafile.A4Off)
        Ei.append(datafile.Ei)
        instrumentCalibrationEf.append(datafile.instrumentCalibrationEf)
        instrumentCalibrationA4.append(datafile.instrumentCalibrationA4)
        instrumentCalibrationEdges.append(datafile.instrumentCalibrationEdges)
        
    I = Marray(I)#np.concatenate(I,axis=0))
    if(files[0].type=='nxs'):
        qx = Marray(qx)#np.concatenate(qx,axis=0))
        qy = Marray(qy)#np.concatenate(qy,axis=0))
        H = Marray(H)#np.concatenate(H,axis=0))
        K = Marray(K)#np.concatenate(K,axis=0))
        L = Marray(L)#np.concatenate(L,axis=0))
        energy = Marray(energy)#np.concatenate(energy,axis=0))
        Norm = Marray(Norm)#np.concatenate(Norm,axis=0))
    else:  
        #print(Norm)
        Norm = Marray(Norm)
        
    Monitor = Marray(Monitor)#np.concatenate(Monitor,axis=0))

    a3 = Marray(a3)#np.concatenate(a3,axis=0))
    a4 = Marray(a4)#np.concatenate(a4,axis=0))

    a3Off = Marray(a3Off)
    a4Off = Marray(a4Off)
    instrumentCalibrationEf = np.array(instrumentCalibrationEf)
    instrumentCalibrationA4 = np.array(instrumentCalibrationA4)
    instrumentCalibrationEdges = np.array(instrumentCalibrationEdges)
    Ei = Marray(Ei)
    if files[0].type=='nxs':
        return I,qx,qy,energy,Norm,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit,H,K,L,mask
    else:
        return I,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit,mask

def assertFile(file):
    """Make sure that file exists for methods to work"""
    if not os.path.isfile(file):
        df = DataFile(file.replace('.nxs','.hdf'))
        con = df.convert(binning=8)
        con.saveNXsqom(file)