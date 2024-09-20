import sys, os
from typing import DefaultDict
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
from MJOLNIR._tools import PointerArray
from MJOLNIR.Data import Mask

import MJOLNIR.Data.Sample
import re
import copy
import platform
from collections import defaultdict

multiFLEXXDetectors = 31*5
reFloat = r'-?\d*\.\d*'
reInt   = r'-?\d'
factorsqrtEK = 0.694692
supportedRawFormats = ['hdf','','dat']
supportedInstruments = ['CAMEA','MultiFLEXX','FlatCone','Bambus']
supportedConvertedFormats = ['nxs']

def cosd(x):
    return np.cos(np.deg2rad(x))

def sind(x):
    return np.sin(np.deg2rad(x))


## Dictionary for holding hdf position of attributes. HDFTranslation['a3'] gives hdf position of 'a3'
HDFTranslation = {'sample':'/entry/sample',
                  'sampleName':'/entry/sample/name',
                  'unitCell':'/entry/sample/unit_cell',
                  'intensity':'entry/data/intensity',
                  'qx':'entry/data/qx',
                  'qy':'entry/data/qy',
                  'QH':'entry/data/h',
                  'QK':'entry/data/k',
                  'QL':'entry/data/l',
                  'energy':'entry/data/en',
                  'normalization':'entry/data/normalization',
                  'mode':'entry/control/mode',
                  'preset':'entry/control/preset',
                  'startTime':'entry/start_time',
                  'hdfMonitor':'entry/control/data',
                  'monitor':'entry/monitor_2/data',
                  'time':'entry/control/time',
                  'endTime':'entry/end_time',
                  'experimentalIdentifier':'entry/experiment_identifier',
                  'comment':'entry/comment',
                  'proposal':'entry/proposal_id',
                  'proposalTitle':'entry/proposal_title',
                  'localContact':'entry/local_contact/name',
                  'proposalUser':'entry/proposal_user/name',
                  'proposalEmail':'entry/proposal_user/email',
                  'user':'entry/user/name',
                  'email':'entry/user/email',
                  'address':'entry/user/address',
                  'affiliation':'entry/user/affiliation',
                  'A3':'entry/sample/rotation_angle',
                  'temperature':'entry/sample/temperature',
                  'magneticField':'entry/sample/magnetic_field',
                  'electricField':'entry/sample/electric_field',
                  'scanCommand':'entry/scancommand',
                  'title':'entry/title',
                  'absoluteTime':'entry/control/absolute_time',
                  'protonBeam':'entry/proton_beam/data',
                  'singleDetector1':'entry/CAMEA/segment_1/data',
                  'singleDetector8':'entry/CAMEA/segment_8/data'
}

HDFTranslationNICOSAlternative = {
                   'temperature':['entry/sample/Ts','entry/sample/Ts/value','entry/sample/se_ts','entry/sample/T','entry/sample/temperature'],
                   'temperature_log':['entry/sample/temperature_log/value','entry/sample/T_log/value'],
                   'temperature_time_log':['entry/sample/temperature_log/time','entry/sample/temperature_log/time','entry/sample/T_log/time'],
                   'magneticField':['entry/sample/B','entry/sample/B/value','entry/sample/magnetic_field'],
                   'ei':'entry/CAMEA/monochromator/energy',
                   'hdfMonitor':['entry/monitor_2/data','entry/control/data']
}

## Default dictionary to perform on loaded data, i.e. take the zeroth element, swap axes, etc

HDFTranslationFunctions = defaultdict(lambda : [])
HDFTranslationFunctions['mode'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['sampleName'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['startTime'] = [['__getitem__',[0]]]
HDFTranslationFunctions['endTime'] = [['__getitem__',[0]]]
HDFTranslationFunctions['experimentalIdentifier'] = [['__getitem__',[0]]]
HDFTranslationFunctions['comment'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['proposal'] = [['__getitem__',[0]]]
HDFTranslationFunctions['proposalTitle'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['localContact'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['proposalUser'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['proposalEmail'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['user'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['email'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['address'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['affiliation'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['scanCommand'] = [['__getitem__',[0]],['decode',['utf8']]]
HDFTranslationFunctions['title'] = [['__getitem__',[0]],['decode',['utf8']]]



HDFInstrumentTranslation = {
                   'A4':'analyzer/polar_angle',
                   'ei':'monochromator/energy',
                   'A4Offset':'analyzer/polar_angle_offset',
                   'counts':'detector/counts'
}

HDFInstrumentTranslationNICOSAlternative = {
                    'counts':'detector/data',
                    'A4':'analyzer/polar_angle'#_raw'

}

HDFInstrumentTranslationFunctions = DefaultDict(lambda : [])
HDFInstrumentTranslationFunctions['counts'] = [['swapaxes',[1,2]]]
HDFInstrumentTranslationFunctions['A4'] = [['reshape',[-1]]]

extraAttributes = ['name','fileLocation','twoTheta']

possibleAttributes = list(HDFTranslation.keys())+list(HDFInstrumentTranslation.keys())+extraAttributes
possibleAttributes.sort(key=lambda v: v.lower())


def getHDFEntry(f,prop,fromNICOS=False):
    if fromNICOS:
        if prop in HDFTranslationNICOSAlternative:
            proplist = HDFTranslationNICOSAlternative[prop]
            if isinstance(proplist,list):
                for subprop in proplist:
                    value = f.get(subprop)
                    if not value is None:
                        break
            else:
                value = f.get(HDFTranslationNICOSAlternative[prop])
            return value

    return f.get(HDFTranslation[prop])

def getHDFInstrumentEntry(instr,prop,fromNICOS=False):
    if fromNICOS:
        if prop in HDFInstrumentTranslationNICOSAlternative:
            value = instr.get(HDFInstrumentTranslationNICOSAlternative[prop]) 
            return value
    
    return instr.get(HDFInstrumentTranslation[prop])

analyzerLimits = {'CAMEA':7,
                  'Bambus': 4,
                  'MultiFLEXX': 4,
                  'FlatCone':1}

detectorLimits = {'CAMEA':103,
                  'Bambus':19,
                  'MultiFLEXX':31,
                  'FlatCone':32}

class DataFile(object):
    """Object to load and keep track of HdF files and their conversions"""
    def __init__(self,fileLocation=None,McStas=False):
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
            elif fileLocation.split('.')[-1]=='dat': # Could be both Bambus and MultiFLEXX
                
                if isMultiFLEXX(fileLocation):
                    self.type = 'MultiFLEXX'
                else:
                    self.type = 'Bambus'
            else:
                raise AttributeError('File is not of type nxs or hdf.')
            self.name = os.path.basename(fileLocation)
            self.fileLocation = os.path.abspath(fileLocation)
            self._binning = 1
            self._mask = False
            self.absolutNormalized = False

            if self.type in ['hdf','nxs']:
                with hdf.File(fileLocation,mode='r') as f:
                    # Find out if file is a NICOS file
                    # NICOS has data saved in /entry/data/data while six has /entry/data/counts
                    if McStas:
                        self.fromNICOS = 0
                    else:
                        self.fromNICOS = checkNICOS(f)
                    instr = getInstrument(f)
                    if self.fromNICOS is None:
                        raise AttributeError('Data File {} has no data in {}/detector/counts. The file might be empty.'.format(self.name,instr.name))
                    self.sample = MJOLNIR.Data.Sample.Sample(sample=getHDFEntry(f,'sample'),recalculateUB=self.fromNICOS)
                    
                    if self.type == 'hdf':
                        value = np.array(getHDFInstrumentEntry(instr,'counts',fromNICOS=self.fromNICOS))
                        fileShape = np.shape(value)
                        if fileShape == ():
                            raise AttributeError('Data File {} has no data in {}/detector/counts. The file might be empty.'.format(self.name,instr.name))
                        
                        if len(fileShape) == 2:
                            value = value.reshape(1,*fileShape)
                        
                        self.I = value.swapaxes(1,2)
                    else:
                        self.I=np.array(getHDFEntry(f,'intensity'))
                        self.counts = np.array(getHDFInstrumentEntry(instr,'counts')).swapaxes(1,2)
                        self.qx=np.array(getHDFEntry(f,'qx'))
                        self.qy=np.array(getHDFEntry(f,'qy'))
                        self.h=np.array(getHDFEntry(f,'QH'))
                        self.k=np.array(getHDFEntry(f,'QK'))
                        self.l=np.array(getHDFEntry(f,'QL'))
                        self.energy=np.array(getHDFEntry(f,'energy'))
                        self.Norm=np.array(getHDFEntry(f,'normalization'))
                    self.MonitorMode = np.array(getHDFEntry(f,'mode',fromNICOS=self.fromNICOS))[0].decode()
                    self.MonitorPreset=np.array(getHDFEntry(f,'preset',fromNICOS=self.fromNICOS))
                    if len(self.MonitorPreset)>1:
                        self.MonitorPreset = self.MonitorPreset[0]             
                    self.startTime = np.array(getHDFEntry(f,'startTime',fromNICOS=self.fromNICOS))[0].decode()
                    self.totalCounts = self.I.sum(axis=(1,2))
                    if self.type == 'hdf':
                        self.Monitor = np.array(getHDFEntry(f,'hdfMonitor',fromNICOS=self.fromNICOS))
                        if not self.MonitorMode == 't' and len(self.Monitor)>1: # If not counting on time and more than one point saved
                            if self.Monitor.flatten()[0]!=self.MonitorPreset and self.startTime[:4]=='2018': # For all data in 2018 with wrong monitor saved
                                self.Monitor = np.ones_like(self.Monitor)*self.MonitorPreset ### TODO: Make Mark save the correct monitor!!
                        if len(self.Monitor)==len(self.I)+1: # Monitor is exactly one longer than I when opening a running file under NICOS
                            self.Monitor = self.Monitor[:-1]
                    else:
                        self.Monitor=np.array(getHDFEntry(f,'monitor',fromNICOS=self.fromNICOS))
                    self.Time = np.array(getHDFEntry(f,'time',fromNICOS=self.fromNICOS))
                    self.endTime = np.array(getHDFEntry(f,'endTime',fromNICOS=self.fromNICOS))[0]
                    expIdx =getHDFEntry(f,'experimentalIdentifier',fromNICOS=self.fromNICOS)
                    if not expIdx is None:
                        self.experimentIdentifier = np.array(expIdx)[0]
                    else:
                        self.experimentIdentifier = 'UNKNOWN'
                    self.comment = np.array(getHDFEntry(f,'comment',fromNICOS=self.fromNICOS))[0]
                    self.proposalId = np.array(getHDFEntry(f,'proposal',fromNICOS=self.fromNICOS))[0]
                    self.proposalTitle = np.array(getHDFEntry(f,'proposalTitle',fromNICOS=self.fromNICOS))[0]

                    localContact = getHDFEntry(f,'localContact',fromNICOS=self.fromNICOS)
                    if not localContact is None:
                        self.localContactName = np.array(localContact)[0]
                    else:
                        self.localContactName = 'UNKNOWN'
                    
                    proposalUserName = getHDFEntry(f,'proposalUser',fromNICOS=self.fromNICOS)
                    if not proposalUserName is None:
                        self.proposalUserName = np.array(proposalUserName)[0]
                        self.proposalUserEmail = np.array(getHDFEntry(f,'proposalEmail'))[0]
                    else:
                        self.proposalUserName = 'UNKNOWN'
                        self.proposalUserEmail = 'UNKNOWN'

                    self.userName = np.array(getHDFEntry(f,'user',fromNICOS=self.fromNICOS))[0]
                    self.userEmail = np.array(getHDFEntry(f,'email',fromNICOS=self.fromNICOS))[0]
                    userAddress = getHDFEntry(f,'address',fromNICOS=self.fromNICOS)
                    if not userAddress is None:
                        self.userAddress = np.array(userAddress)[0]
                    else:
                        self.userAddress = 'UNKNOWN'

                    userAffiliation = getHDFEntry(f,'affiliation',fromNICOS=self.fromNICOS)
                    if not userAffiliation is None:
                        self.userAffiliation = np.array(userAffiliation)[0]
                    else:
                        self.userAffiliation = 'UNKNOWN'
                    
                    
                    
                    try:
                        self.singleDetector1 = np.array(getHDFEntry(f,'singleDetector1',fromNICOS=self.fromNICOS))[0]
                    except:
                        pass

                    try:
                        self.singleDetector8 = np.array(getHDFEntry(f,'singleDetector8',fromNICOS=self.fromNICOS))[0]
                    except:
                        pass
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
                    attributes = ['d_spacing','nominal_energy','polar_angle','polar_angle_offset','polar_angle_raw']
                    values = ['analyzer'+x.replace('_',' ').title().replace(' ','') for x in attributes]
                    for att,value in zip(attributes,values):
                        setattr(self,value,np.array(f.get('entry/CAMEA/analyzer/{}'.format(att))))
                    self.analyzerType = np.array(f.get('entry/CAMEA/analyzer/type'))[0]
                    self.analyzerSelection = int(np.array(f.get('entry/CAMEA/analyzer/analyzer_selection'))[0])
                    self.detectorSelection = int(np.array(f.get('entry/CAMEA/detector/detector_selection'))[0])

                    instr = getInstrument(f)
                    self.instrument = instr.name.split('/')[-1]
                    self.possibleBinnings = np.array([int(x[-1]) for x in np.array(instr) if x[:5]=='calib'])
                    self.Ei = np.array(getHDFInstrumentEntry(instr,'ei',fromNICOS=self.fromNICOS))
                    self.A3 = np.array(getHDFEntry(f,'A3',fromNICOS=self.fromNICOS))
                    self.A4 = np.array(getHDFInstrumentEntry(instr,'A4',fromNICOS=self.fromNICOS)).reshape(-1)
                    self.A4Off = np.array(getHDFInstrumentEntry(instr,'A4Offset',fromNICOS=self.fromNICOS))
                    if self.fromNICOS:
                        self.twotheta = copy.deepcopy(self.A4) - self.A4Off
                    else:
                        self.twotheta = self.A4-self.A4Off

                    # As of 2022 and implementation of NICOS temperature and magnetic field is saved
                    # in another subfolder in the data file (e.g. entry/sample/B/value)
                    self.temperature = np.array(getHDFEntry(f,'temperature',fromNICOS=self.fromNICOS))
                    self.magneticField = np.array(getHDFEntry(f,'magneticField',fromNICOS=self.fromNICOS))
                    self.electricField = np.array(getHDFEntry(f,'electricField',fromNICOS=self.fromNICOS))
                    if self.fromNICOS and self.temperature is None: # Only use interpolated temperature if no temperature was found
                        self.temperature_log = np.array(getHDFEntry(f,'temperature_log',fromNICOS=self.fromNICOS))
                        self.temperature_time_log = np.array(getHDFEntry(f,'temperature_time_log',fromNICOS=self.fromNICOS))
                        if not self.temperature_log is None:
                            timeSteps = self.absoluteTime-self.absoluteTime[0]
                            try:
                                self.temperature_log = self.temperature_log[:len(self.temperature_time_log)]

                                self.temperature =[np.mean(self.temperature_log[np.logical_and(self.temperature_time_log>tStart,self.temperature_time_log<tStop)]) for tStart,tStop in zip(timeSteps,timeSteps[1:])]
                                ## Above adds all but last temperature interval 
                                self.temperature.append(np.mean(self.temperature_log[self.temperature_time_log>timeSteps[-1]]))
                                self.temperature = np.array(self.temperature)
                            except TypeError: # no length of temperature_time_log, i.e. the log is not present
                                self.temperature = np.asarray(len(self.absoluteTime)*[None])

                    self.scanCommand = np.array(getHDFEntry(f,'scanCommand',fromNICOS=self.fromNICOS))[0]
                    try:
                        self.scanParameters,self.scanValues,self.scanUnits,self.scanDataPosition = getScanParameter(self,f)
                    except KeyError:
                        warnings.warn("Couldn't load scan parameters")
                        self.scanParameters = []
                        self.scanValues= np.array([[None]])
                        self.scanUnits =''
                        self.scanDataPosition = ''
                    if len(self.scanParameters) == 1 and self.scanParameters[0] == 'rotation_angle' and len(self.A4)>1 and np.all(np.isclose(self.A4,self.A4[0],atol=0.01)): 
                        # If all A4 values are the same, out 2t has been written on instrument computer
                        # and because of this, six saves 2t and not a4. Solution: set A4Offset to 0 and
                        # use only the first value for A4 when an A3 (=rotation_angle) scan is made and 
                        self.A4 = np.array([self.A4[0]])
                        self.analyzerPolarAngle = self.A4
                        self.A4Off = np.array([0.0])
                        self.analyzerPolarAngleOffset = self.A4Off
                    else:
                        self.A4Off = np.array(getHDFInstrumentEntry(instr,'A4Offset',fromNICOS=self.fromNICOS))
                        

                    try:
                        self.A3Off = self.sample.A3Off#np.array(f.get('entry/sample/rotation_angle_zero'))  
                    except:
                        self.A3Off = [0.0]
                    if self.type == 'hdf':
                        
                        if len(self.possibleBinnings) == 0:
                            self.possibleBinnings = [None]
                            self.binning = None
                        else:
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

                    
                    
                    
                    if self.type == 'nxs':
                        self.original_fileLocation = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/rawdata'))[0].decode()
                    self.title = np.array(getHDFEntry(f,'title'))[0]

                    self.absoluteTime = np.array(getHDFEntry(f,'absoluteTime',fromNICOS=self.fromNICOS))
                    self.protonBeam = np.array(getHDFEntry(f,'protonBeam',fromNICOS=self.fromNICOS))
                    


                    if self.type == 'hdf':
                        ###################
                        #self.I[:,:,:150]=0#
                        ###################
                        pass
                    self.mask = np.zeros_like(self.I,dtype=bool)
                    if self.binning == 8:
                        self.mask[:,:,:2] = True
            elif self.type == 'MultiFLEXX': # type is multiFLEXX
                self.loadMultiFLEXXData(fileLocation)
            elif self.type == 'Bambus':
                self.loadBambusData(fileLocation)
            try:
                self.scanSteps = self.scanValues.shape[1]
            except:
                pass
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
            self.instrument = 'CAMEA'
            self.type = 'hdf'
            self.I = np.array([])
            self._binning = 1
            self._mask = False
            self.absolutNormalized = False
            self.scanParameters = ''
            self.scanUnits = ''
            self.scanValues = np.array([])
            self.Monitor = np.array([])
            self._A3Off = 0.0
            self._A3 = np.array([])
            self._A4Off = 0.0
            self._A4 = np.array([])
            self.Ei = np.array([])
            self.possibleBinnings = [1,3,5,8]
            calibrations = []
            for binning in self.possibleBinnings:
                fileName = getattr(MJOLNIR,'__CAMEANormalizationBinning{}__'.format(binning))
                calib = np.loadtxt(fileName,delimiter=',',skiprows=3)
                calibrations.append([calib[:,3:7],calib[:,-1],calib[:,7:9]])
            self.instrumentCalibrations = np.array(calibrations,dtype=object)
            self.loadBinning(self.binning)


        if self.instrument == 'CAMEA':
            self.EPrDetector = 8 
        elif self.type in ['MultiFLEXX','FlatCone','Bambus']:
            self.EPrDetector = 1
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
    def dasel(self):
        if hasattr(self,'detectorSelection') and hasattr(self,'analyzerSelection'):
            return self.detectorSelection,self.analyzerSelection
        else:
            return None,None

    @dasel.getter
    def dasel(self):
        if hasattr(self,'detectorSelection') and hasattr(self,'analyzerSelection'):
            return self.detectorSelection,self.analyzerSelection
        else:
            return None,None

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
        if hasattr(self,'I') and len(self.I.shape) > 1:
            if not np.all(self.I.shape == mask.shape):
                raise AttributeError('Shape of provided mask {} does not match shape of data {}.'.format(mask.shape,self.I.shape))
        self._mask = mask

    @property
    def shape(self):
        return self.I.shape

    @shape.getter
    def shape(self):
        return self.I.shape

    @shape.setter
    def shape(self,shape):
        raise AttributeError('Shape of a DataFile cannot be set.')

    @property
    def size(self):
        return self.I.size

    @size.getter
    def size(self):
        return self.I.size

    @size.setter
    def size(self,size):
        raise AttributeError('size of a DataFile cannot be set.')

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
        
        self.fileLocation = fileLocation
        self.possibleBinnings = [1] # Standard value (1 energy/detector)
        self.binning = 1
        self.instrument = 'MultiFLEXX'
        self.type = self.instrument


        ## No dasel
        self.analyzerSelection = 0
        self.detectorSelection = 0

        with open(fileLocation) as f:
            dataString = f.readlines()
        if len(dataString) == 0:
            raise AttributeError('DataFile not found! Looking for "{:}"'.format(f))


        # Format for starting time is 2021-11-26 10:45:05
        searchString = r'([0-9]+)-([0-9]+)-([0-9]+)\s([0-9]+):([0-9]+):([0-9]+)'

        matches = re.findall(searchString,dataString[0])
        if len(matches) == 0:
            self.startTime = 'UNKNONW'
        else:
            self.startTime = '-'.join(matches[0][:3])+' '+':'.join(matches[0][3:])
        self.endTime = 'UNKNOWN'

        if calibrationFile is None:
            calibrationFile = MJOLNIR.__multiFLEXXNormalization__

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

        ## Load header data:

        ## Format is "#    param : value"

        parameters = {}

        for I,line in enumerate(dataString):
            if line[:2] == '# ':
                param,*value = [x.strip() for x in line.split(':')]
                param = param.replace('#','').strip()
                if len(value)>1:
                    value = ':'.join(value)
                parameters[param] = value[0]
            else:
                if line.find('DATA_:') >-1:
                    dataline = I
                    break

        # Convert from ch number in data file to detector,energy (A1,A2,...,20D,20E) LEGACY from Bambus
        detectorMap = np.arange(155)#np.array([ 0,  1,  2, 27,  3, 50, 51, 52, 77, 53,  4,  5,  6, 28,  7, 54, 55,
            #56, 78, 57,  8,  9, 10, 29, 11, 58, 59, 60, 79, 61, 12, 13, 14, 30,
            #15, 62, 63, 64, 80, 65, 16, 17, 18, 31, 19, 66, 67, 68, 81, 69, 20,
            #21, 22, 32, 23, 70, 71, 72, 82, 73, 24, 25, 26, 33, 37, 74, 75, 76,
            #83, 87, 38, 39, 40, 34, 41, 88, 89, 90, 84, 91, 42, 43, 44, 35, 45,
            #92, 93, 94, 85, 95, 46, 47, 48, 36, 49, 96, 97, 98, 86, 99])


        # Load data
        titles = dataString[dataline+1].strip().split()
        units =  dataString[dataline+2].strip().split()


        # remove the '# ' part of the first title and unit
        titles = titles[1:]#[0] = titles[0].replace('#','').strip()
        units = units[1:]#units[0] = units[0].replace('#','').strip()


        # Find splitter in data set destinguishing between title and chanels
        
        splitter = titles.index('ch0')-1


        # Shape becomes [scan points, timer + mon1 + mon2 + chsum + 100 detectors + events]
        data = np.asarray([np.asarray(line.strip().split(),dtype=float) for line in dataString[dataline+3:-1]])
        scanData = data[:,:splitter]
        data = data[:,splitter+1:]
        
        #data = [[float(x) for x in line.strip().split()[splitter+1:]] for line in dataString[dataline+3:-1]]
        #scanData = [[float(x) for x in line.strip().split()[:splitter]] for line in dataString[dataline+3:-1]]

        if len(data[-1]) == 0: # empty last line
            data = np.array(data[:-1])
            scanData = np.array(scanData[:-1])
        else:
            data = np.array(data)
            scanData = np.array(scanData)
            
        
        #reshape into [scan points, 31 wedges, 5 energies]
        self.I = np.array([d[detectorMap] for d in data]).reshape(-1,31*5,1)


        self.timer = data[:,0].astype(float)
        self.mon1,self.mon2,self.ctot,self.events = data[:,[1,2,3,-1]].astype(int).T

        self.scanParameters = np.asarray(titles[:splitter])
        self.scanUnits = np.asarray(units[:splitter])
        self.scanValues = np.asarray(scanData[:,:splitter]).T.astype(float)

        self.scanParameters[0] = self.scanParameters[0].replace('#','').strip()
        self.scanUnits[0] = self.scanUnits[0].replace('#','').strip()

        ### Move last scan parameter to first position in lists

        self.scanParameters = self.scanParameters[[-1,*range(len(self.scanParameters)-1)]]
        self.scanUnits = self.scanUnits[[-1,*range(len(self.scanParameters)-1)]]
        self.scanValues = self.scanValues[[-1,*range(len(self.scanParameters)-1)]]


        self.__dict__.update(parameters)
        
        #### Import sample correctly

        sample = {}
        def extractFromParenthesis(string):
            return list(map(float,re.findall(r'(?<=\().*?(?=\))',string)[0].split(',')))


        cell = np.concatenate([extractFromParenthesis(dat) for dat in [self.Sample_lattice,self.Sample_angles]])

        for param,value in zip(['a','b','c','alpha','beta','gamma'],cell):
            sample[param] = value


        q1 = np.array(extractFromParenthesis(self.Sample_orient1))
        q2 = np.array(extractFromParenthesis(self.Sample_orient2))



        Ei = 10.0 # TODO: What is the good solution here? Dummy incoming energy needed to calcualte UB
        k = np.sqrt(Ei)*factorsqrtEK

        Cell = TasUBlib.calcCell(cell)
        B = TasUBlib.calculateBMatrix(Cell)

        A3offset = float(self.Sample_psi0.split(' ')[0])

        A41 = TasUBlib.calTwoTheta(B,[*q1,Ei,Ei],-1)
        A31 = TasUBlib.calcTheta(k,k,A41)+A3offset
        A42 = TasUBlib.calTwoTheta(B,[*q2,Ei,Ei],-1)
        A32 = TasUBlib.calcTheta(k,k,A42)

        planeVector1 = list(q1)
        planeVector1.append(A31) # A3 
        planeVector1.append(A41) # A4
        [planeVector1.append(0.0) for _ in range(2)]# append values for gonios set to zero
        planeVector1.append(Ei)
        planeVector1.append(Ei)

        planeVector2 = list(q2)
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

        sample['name'] = self.Sample_samplename


        self.sample = MJOLNIR.Data.Sample.Sample(**sample)

        for sP,sV in zip(self.scanParameters,self.scanValues):
            setattr(self,sP,sV)


        monoQ = np.array(float(self.mono_value.split(' ')[0]))
        self.Ei = np.array([float(self.Ei_value.split(' ')[0])])
        # self.A4 = np.array([float(self.stt_value.split(' ')[0])])

        # if isinstance(self.sth_value,str):
        #     self._A3 = np.array(float(self.sth_value.split(' ')[0]))
        # else:
        #     self._A3 = self.sth_value

        self._A4Off = 0.0
        self._A3Off = 0.0

        # self.countingTime = np.array([self.etime[0],*np.diff(self.etime)])


        nameSwaps = [['filename','name'],
                    ['info','scanCommand'],
                    ['Exp_remark','comment'],
                    
                    #['MAG','magneticField'],
                    ['mon1','Monitor'],
                    ['timer','Time'],
                    #['TT','temperature'],
                    #['comnd','scanCommand'],
                    #['instr','instrument'],
                    #['EI','Ei'],
                    
                    ['Exp_localcontact','localContact'],
                    ['Exp_proposal','proposalID'],
                    ['Exp_title','title'],
                    ['Exp_users','userName'],

                    
                    ['timer','countingTime'],
                    
                    #['DA','analyzerDSpacing'],
                    #['expno','experimentIdentifier'],
                    ['localContact','localContactName'],
                    #['date','startTime']
                    ]


        if not hasattr(self,'electricField'):
            self.electricField = None
        if not hasattr(self,'magneticField'):
            self.magneticField = None

        if not hasattr(self,'temperature'):
            self.temperature = None
            


        ## Format scanCommand

        def updateKeyName(obj,key,newName):
            if not hasattr(obj,key):
                return
            setattr(obj,newName,getattr(obj,key))
            delattr(obj,key)


        for pair in nameSwaps:
            updateKeyName(self,pair[0],pair[1])
            
        self.twotheta = self.stt
        self.A4 = self.stt
        self.A3 = self.sth

        self._mask = np.zeros_like(self.I)
        self.EPrDetector = 1

        self.Monitor = np.ones_like(self.Monitor)


    def loadBambusData(self,fileLocation,calibrationFile=None):
        self.fileLocation = fileLocation
        self.possibleBinnings = [1] # Standard value (1 energy/detector)
        self.binning = 1
        self.instrument = 'Bambus'
        self.type = self.instrument
        

        ## No dasel
        self.analyzerSelection = 0
        self.detectorSelection = 0

        with open(fileLocation) as f:
            dataString = f.readlines()



        # Format for starting time is 2021-11-26 10:45:05
        searchString = r'([0-9]+)-([0-9]+)-([0-9]+)\s([0-9]+):([0-9]+):([0-9]+)'
        
        matches = re.findall(searchString,dataString[0])
        if len(matches) == 0:
            self.startTime = 'UNKNONW'
        else:
            self.startTime = '-'.join(matches[0][:3])+' '+':'.join(matches[0][3:])
        self.endTime = 'UNKNOWN'

        if calibrationFile is None:
            calibrationFile = MJOLNIR.__bambusNormalization__

        detectors = 20
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

        ## Load header data:

        ## Format is "#    param : value"

        parameters = {}

        for I,line in enumerate(dataString):
            if line[:2] == '# ':
                param,*value = [x.strip() for x in line.split(':')]
                param = param.replace('#','').strip()
                if len(value)>1:
                    value = ':'.join(value)
                parameters[param] = value[0]
            else:
                if line.find('### Scan data') >-1:
                    dataline = I
                    break

        # Convert from ch number in data file to detector,energy (A1,A2,...,20D,20E)
        detectorMap = np.array([ 0,  1,  2, 27,  3, 50, 51, 52, 77, 53,  4,  5,  6, 28,  7, 54, 55,
            56, 78, 57,  8,  9, 10, 29, 11, 58, 59, 60, 79, 61, 12, 13, 14, 30,
            15, 62, 63, 64, 80, 65, 16, 17, 18, 31, 19, 66, 67, 68, 81, 69, 20,
            21, 22, 32, 23, 70, 71, 72, 82, 73, 24, 25, 26, 33, 37, 74, 75, 76,
            83, 87, 38, 39, 40, 34, 41, 88, 89, 90, 84, 91, 42, 43, 44, 35, 45,
            92, 93, 94, 85, 95, 46, 47, 48, 36, 49, 96, 97, 98, 86, 99])


        # Load data
        titles = dataString[dataline+1].strip().split('\t')
        units =  dataString[dataline+2].strip().split('\t')


        # remove the '# ' part of the first title and unit
        titles[0] = titles[0].replace('#','').strip()
        units[0] = units[0].replace('#','').strip()


        # Find splitter in data set ';'
        
        splitter = titles.index(';')
        

        # Shape becomes [scan points, timer + mon1 + mon2 + chsum + 100 detectors + events]
        data = [[float(x) for x in line.split('\t')[splitter+1:]] for line in dataString[dataline+3:-1]]
        scanData = [[float(x) for x in line.split('\t')[:splitter]] for line in dataString[dataline+3:-1]]
        
        if len(data[-1]) == 0: # empty last line
            data = np.array(data[:-1])
            scanData = np.array(scanData[:-1])
        else:
            data = np.array(data)
            scanData = np.array(scanData)
            

        #reshape into [scan points, 20 wedges, 5 energies]
        self.I = np.array([d[detectorMap] for d in data[:,4:-1]]).reshape(-1,20*5,1)


        self.timer = data[:,0].astype(float)
        self.mon1,self.mon2,self.ctot,self.events = data[:,[1,2,3,-1]].astype(int).T

        self.scanParameters = np.asarray(titles[:splitter])
        self.scanUnits = np.asarray(units[:splitter])
        self.scanValues = np.asarray(scanData[:,:splitter]).T.astype(float)

        self.scanParameters[0] = self.scanParameters[0].replace('#','').strip()
        self.scanUnits[0] = self.scanUnits[0].replace('#','').strip()

        ### Move last scan parameter to first position in lists

        self.scanParameters = self.scanParameters[[-1,*range(len(self.scanParameters)-1)]]
        self.scanUnits = self.scanUnits[[-1,*range(len(self.scanParameters)-1)]]
        self.scanValues = self.scanValues[[-1,*range(len(self.scanParameters)-1)]]


        self.__dict__.update(parameters)




        #### Import sample correctly

        sample = {}
        def extractFromParenthesis(string):
            return list(map(float,re.findall(r'(?<=\().*?(?=\))',string)[0].split(',')))


        cell = np.concatenate([extractFromParenthesis(dat) for dat in [self.Sample_lattice,self.Sample_angles]])

        for param,value in zip(['a','b','c','alpha','beta','gamma'],cell):
            sample[param] = value


        q1 = np.array(extractFromParenthesis(self.Sample_orient1))
        q2 = np.array(extractFromParenthesis(self.Sample_orient2))



        Ei = 10.0 # TODO: What is the good solution here? Dummy incoming energy needed to calcualte UB
        k = np.sqrt(Ei)*factorsqrtEK

        Cell = TasUBlib.calcCell(cell)
        B = TasUBlib.calculateBMatrix(Cell)

        A3offset = float(self.Sample_psi0.split(' ')[0])

        A41 = TasUBlib.calTwoTheta(B,[*q1,Ei,Ei],-1)
        A31 = TasUBlib.calcTheta(k,k,A41)+A3offset
        A42 = TasUBlib.calTwoTheta(B,[*q2,Ei,Ei],-1)
        A32 = TasUBlib.calcTheta(k,k,A42)

        planeVector1 = list(q1)
        planeVector1.append(A31) # A3 
        planeVector1.append(A41) # A4
        [planeVector1.append(0.0) for _ in range(2)]# append values for gonios set to zero
        planeVector1.append(Ei)
        planeVector1.append(Ei)

        planeVector2 = list(q2)
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

        sample['name'] = self.Sample_samplename


        self.sample = MJOLNIR.Data.Sample.Sample(**sample)

        for sP,sV in zip(self.scanParameters,self.scanValues):
            setattr(self,sP,sV)
        

        monoQ = np.array(float(self.mono_value.split(' ')[0]))
        self.Ei = np.array(float(self.Ei_value.split(' ')[0]))
        self.A4 = np.array([float(self.stt_value.split(' ')[0])])

        if isinstance(self.sth_value,str):
            self._A3 = np.array(float(self.sth_value.split(' ')[0]))
        else:
            self._A3 = self.sth_value

        self.A4Off = 0.0
        self.A3Off = None

        self.countingTime = np.array([self.etime[0],*np.diff(self.etime)])
        

        nameSwaps = [['filename','name'],
                    ['info','scanCommand'],
                    ['Exp_remark','comment'],
                    
                    #['MAG','magneticField'],
                    ['mon1','Monitor'],
                    ['timer','Time'],
                    #['TT','temperature'],
                    #['comnd','scanCommand'],
                    #['instr','instrument'],
                    #['EI','Ei'],
                    
                    ['Exp_localcontact','localContact'],
                    ['Exp_proposal','proposalID'],
                    ['Exp_title','title'],
                    ['Exp_users','userName'],

                    
                    
                    
                    #['DA','analyzerDSpacing'],
                    #['expno','experimentIdentifier'],
                    ['localContact','localContactName'],
                    #['date','startTime']
                    ]


        if not hasattr(self,'electricField'):
            self.electricField = None
        if not hasattr(self,'magneticField'):
            self.magneticField = None

        if not hasattr(self,'temperature'):
            self.temperature = None
            
        self.twotheta = self.A4


        ## Format scanCommand

        def updateKeyName(obj,key,newName):
            if not hasattr(obj,key):
                return
            setattr(obj,newName,getattr(obj,key))
            delattr(obj,key)


        for pair in nameSwaps:
            updateKeyName(self,pair[0],pair[1])

        self._mask = np.zeros_like(self.I)

    @_tools.KwargChecker()
    def calcualteDataIndexFromDasel(self,detectorSelection=None,analyzerSelection=None):
        if detectorSelection is None:
            detectorSelection = self.detectorSelection
        if analyzerSelection is None:
            analyzerSelection = self.analyzerSelection

        if self.instrument == 'CAMEA':
            if detectorSelection > detectorLimits['CAMEA']:
                raise AttributeError('Provided detectorSelection is out of range. Recieved {} and instrument has {} detectors enumerated [0-{}]'.format(detectorSelection,detectorLimits['CAMEA']+1,detectorLimits['CAMEA']))
            if analyzerSelection > analyzerLimits['CAMEA']:
                raise AttributeError('Provided analyzerSelection is out of range. Recieved {} and instrument has {} analyzers enumerated [0-{}]'.format(analyzerSelection,analyzerLimits['CAMEA']+1,analyzerLimits['CAMEA']))
            analyzerPixels = self.instrumentCalibrationEdges[analyzerSelection+detectorSelection*8] 
            return detectorSelection,range(int(analyzerPixels[0]),int(analyzerPixels[1]))
        
        else:
            if detectorSelection > detectorLimits[self.instrument]:
                raise AttributeError('Provided detectorSelection is out of range. Recieved {} and instrument has {} detectors enumerated [0-{}]'.format(detectorSelection,detectorLimits[self.instrument]+1,detectorLimits[self.instrument]))
            if analyzerSelection > analyzerLimits['Bambus']:
                raise AttributeError('Provided analyzerSelection is out of range. Recieved {} and instrument has {} analyzers enumerated [0-{}]'.format(analyzerSelection,analyzerLimits[self.instrument]+1,analyzerLimits[self.instrument]))
            # Calcualte detector number from detector and analyzer. As there are 5 energies(analyzers) per wedge
            return detectorSelection*(analyzerLimits[self.instrument]+1)+analyzerSelection,[0] # Last index is to be able to sum over
    
    
    @_tools.KwargChecker()
    def convert(self,binning=None,printFunction=None):
        if self.instrument == 'CAMEA' or self.type in ['MultiFLEXX','FlatCone','Bambus']:
            if binning is None:
                binning = self.binning
        else:
            raise AttributeError('Instrument type of data file not understood. {} was given.'.format(self.instrument))
        self.loadBinning(binning)

        if printFunction is None:
            printFunction = warnings.warn
        
        EfNormalization = self.instrumentCalibrationEf.copy()
        A4Normalization = self.instrumentCalibrationA4.copy()#np.array(instrument.get('calib{}/a4offset'.format(str(binning))))
        EdgesNormalization = self.instrumentCalibrationEdges.copy()#np.array(instrument.get('calib{}/boundaries'.format(str(binning))))
        self.counts = self.I.copy()
        Data = self.I.copy()#np.array(instrument.get('detector/data'))


        detectors = Data.shape[1]
        steps = Data.shape[0]
        
        if self.type in ['MultiFLEXX','FlatCone','Bambus']:
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
        A4=A4.reshape(detectors,binning*self.EPrDetector,order='C')

        PixelEdge = EdgesNormalization.reshape(detectors,self.EPrDetector,binning,2).astype(int)
        A4File = self.A4.copy()
        
        A4File = A4File.reshape((-1,1,1))

        A4Mean = (A4.reshape((1,detectors,binning*self.EPrDetector))+np.deg2rad(A4File-A4Zero))
        
        Intensity=np.zeros((Data.shape[0],Data.shape[1],self.EPrDetector*binning),dtype=int)
        for i in range(detectors): # for each detector
            for j in range(self.EPrDetector):
                for k in range(binning):
                    Intensity[:,i,j*binning+k] = np.sum(Data[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)

        
        EfMean = EfNormalization[:,1].reshape(1,A4.shape[0],self.EPrDetector*binning)
        EfNormalization = EfNormalization[:,0]#.reshape(1,A4.shape[0],EPrDetector*binning)#
        #EfNormalization = EfNormalization[:,0]*(np.sqrt(2*np.pi)*EfNormalization[:,2])
        
        EfNormalization.shape = (1,A4.shape[0],self.EPrDetector*binning)
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
        Monitor = Monitor*np.ones((1,detectors,self.EPrDetector*binning))
        Normalization = EfNormalization*np.ones((steps,1,1))
        mask = np.zeros_like(Intensity) # TODO: Redo???
       
        ###########################
        #Monitor[:,:,:binning] = 0 #
        ###########################


        convFile = DataFile(self) # Copy everything from old file
        del convFile.counts
        updateDict = {'I':Intensity,'Monitor':Monitor,'qx':QX,'qy':QY,'energy':DeltaE,'binning':binning,'Norm':Normalization,
        'h':H,'k':K,'l':L,'type':'nxs','fileLocation':None,'original_fileLocation':self.fileLocation,'name':self.name.replace('.hdf','.nxs'),'mask':mask}
        convFile.updateProperty(updateDict)

        if convFile.type == 'nxs' and convFile.binning == 8:
            convFile.mask = np.zeros_like(convFile.I,dtype=bool)
            convFile.mask[:,:,:2] = True


        if convFile.instrument == 'CAMEA':
            defectTubes = np.arange(104)[np.any(np.isnan(convFile.instrumentCalibrationA4.reshape(104,-1)),axis=1)] 

            if len(defectTubes)>0: # if any tubes are defect
                    if len(defectTubes)>1:
                        printFunction('Detector tubes {} masked'.format(', '.join([str(x) for x in defectTubes])))
                    else:
                        printFunction('Detector tube {} masked'.format(defectTubes[0]))

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
        



    def saveNXsqom(self,saveFileName):
        """Save converted file into an NXsqom.

        Args:

            - saveFileName (string): File name to be saved into.

        """

        if not self.__hasattr__('original_fileLocation'):
            raise AttributeError('Data file does not have link to the original file. This is needed to make a complete copy when creating nxs-files')
        if not self.type =='nxs':
            raise AttributeError('Only nxs typed files can be saved as nxs-files.')

        datafile = self.original_fileLocation
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
        with hdf.File(saveFileName,'w') as fd:
            with hdf.File(datafile,'r') as fs:
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
                
                rawdata = proc.create_dataset('rawdata',shape=(1,),dtype='S200',data=np.string_(os.path.realpath(datafile)))
                rawdata.attrs['NX_class']=b'NX_CHAR'

                normalizationString = proc.create_dataset('binning',shape=(1,),dtype='int32',data=binning)
                normalizationString.attrs['NX_class']=b'NX_INT'
                
                data = fd.get('entry/data')
                
                fileLength = Intensity.shape
                
                Int = data.create_dataset('intensity',shape=(fileLength),dtype='int32',data=Intensity)
                Int.attrs['NX_class']='NX_INT'

                if self.fromNICOS:
                    counts = np.array(fd.get('entry/data/data'))
                    Int = data.create_dataset('counts',dtype='int32',data=counts)
                    Int.attrs['NX_class']='NX_INT'
                    instr = getInstrument(fd)
                    Int = instr.create_dataset('detector/counts',dtype='int32',data=counts)
                    Int.attrs['NX_class']='NX_INT'
                
                monitor = data.create_dataset('monitor',shape=(fileLength),dtype='int32',data=Monitor)
                monitor.attrs['NX_class']=b'NX_INT'
                
                if fd.get('entry/monitor_2') is None:
                    mon = fd.create_group('entry/monitor_2')
                    monitor = mon.create_dataset('data',shape=(fileLength),dtype='int32',data=Monitor)
                    monitor.attrs['NX_class']=b'NX_INT'
                else:
                    pass

                

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

                #fd.close()

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
            if binning is None:
                continue
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

        self.instrumentCalibrations = np.array([c for c in calibrations.values()],dtype=object)
        self.possibleBinnings = np.array(list(calibrations.keys()))

    def updateSampleParameters(self,unitCell):
        """Update unit cell parameters and corresponding UB matrix

        Args:

            - unitCell (list): List of cell parameters (a,b,c,alpha,beta,gamma)

        """
        self.sample.updateSampleParameters(unitCell=unitCell) 

    def calculateCurratAxeMask(self,BraggPeaks,dqx=None,dqy=None,dH=None,dK=None,dL=None,spurionType='both',maskInside=True):
        """Generate an elliptical mask centered on the Currat-Axe spurion.
        
        Args:
        
            - BraggPeaks (list): List of Bragg peaks to be used for mask. Shape is given by [[H1,K1,L1], [H2,K2,L2], ... [HN,KN,LN]]
            
        Kwargs:

            - dqx (float): Radius used for masking along qx (default None)

            - dqy (float): Radius used for masking along qy (default None)

            - dH (float): Radius used for masking along H (default None)

            - dK (float): Radius used for masking along K (default None)

            - dL (float): Radius used for masking along L (default None)

            - spurionType (str): Either monochromator, analyser or both (default 'both')

            - maskInside (bool): If true, points inside is masked otherwise outside (default True)

        Returns:

            - mask (list): Boolean numpy array with shape equal to self.I.shape

        Note:

            If either dqx or dqy is None, utilizes the dH, dK, dL instead.

        """

        if np.any([dqx is None,dqy is None]):
            if np.any([x is None for x in [dH,dK,dL]]):
                raise AttributeError('Provided masking radius not understood. Either dqx and dqy are to be given or dH, dK, and dL.\nRecieved: \n{}'.format('\n'.join(['    {}\t= {}'.format(x,v) for x,v in zip(['dqx','dqy','dH','dK','dL'],[dqx,dqy,dH,dK,dL])])))

            rlu = True
            factor = np.array([1.0,dH/dK,dH/dL]).reshape(3,1,1,1)
            
        else:
            rlu = False
            factor = dqx/dqy

        if not spurionType.lower() in ['monochromator','analyser','both']:
            raise AttributeError('Provided spurion type not understood. Received '+spurionType+' but expected '+', '.join(['monochromator','analyser','both']))

        s = self.sample
    
        Ei = self.Ei[0]
        Ef = self.instrumentCalibrationEf[:,1].reshape(-1,self.EPrDetector*self.binning).mean(axis=0)
        if spurionType.lower() in ['monochromator','both']:
            monoQx,monoQy = s.CurratAxe(Ei=Ei,Ef=Ef,Bragg=BraggPeaks,HKL=False)[:,0,:,:2].transpose(2,0,1)
        if spurionType.lower() in ['analyser','both']:
            anaQx,anaQy = s.CurratAxe(Ei=Ei,Ef=Ef,Bragg=BraggPeaks,HKL=False,spurionType='Analyser')[:,0,:,:2].transpose(2,0,1) # Into shape (2,len(Bragg),len(Ef))
        
        if not rlu:
            qx = self.qx[:,:,:]
            qy = self.qy[:,:,:]
            monoInside = None
            anaInside = None
            # Reshape to fit qx,qy from data file
            if spurionType.lower() in ['monochromator','both']:
                monoQx = monoQx.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                monoQy = monoQy.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                for localMonoQx,localMonoQy in zip(monoQx,monoQy):
                    if monoInside is None:
                        monoInside = np.linalg.norm([qx-localMonoQx,(qy-localMonoQy)*factor],axis=0)<dqx
                        monoInside.dtype = bool
                    else:
                        monoInside += np.linalg.norm([qx-localMonoQx,(qy-localMonoQy)*factor],axis=0)<dqx
            
            if spurionType.lower() in ['analyser','both']:

                anaQx = anaQx.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                anaQy = anaQy.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                for localAnaQx,localAnaQy in zip(anaQx,anaQy):
                    if anaInside is None:
                        anaInside = np.linalg.norm([qx-localAnaQx,(qy-localAnaQy)*factor],axis=0)<dqx
                        anaInside.dtype = bool
                    else:
                        anaInside += np.linalg.norm([qx-localAnaQx,(qy-localAnaQy)*factor],axis=0)<dqx
            
        else:
            H = self.h[:,:,:]
            K = self.k[:,:,:]
            L = self.l[:,:,:]

            monoInside = None
            anaInside = None

            if spurionType.lower() in ['monochromator','both']:
                monoH,monoK,monoL = s.calculateQxQyToHKL(monoQx,monoQy)
                
                monoH = monoH.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                monoK = monoK.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                monoL = monoL.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                
                for localMonoH,localMonoK,localMonoL in zip(monoH,monoK,monoL):
                    if monoInside is None:
                        monoInside = np.linalg.norm(np.array([H-localMonoH,K-localMonoK,L-localMonoL])*factor,axis=0)<dH
                        monoInside.dtype = bool
                    else:
                        monoInside += np.linalg.norm(np.array([H-localMonoH,K-localMonoK,L-localMonoL])*factor,axis=0)<dH

            if spurionType.lower() in ['analyser','both']:
                anaH,anaK,anaL = s.calculateQxQyToHKL(anaQx,anaQy)
                
                anaH = anaH.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                anaK = anaK.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)
                anaL = anaL.transpose(1,0).reshape(1,1,len(Ef),-1).transpose(3,0,1,2)

                for localAnaH,localAnaK,localAnaL in zip(anaH,anaK,anaL):
                    if anaInside is None:
                        anaInside = np.linalg.norm(np.asarray([H-localAnaH,K-localAnaK,L-localAnaL])*factor,axis=0)<dH
                        anaInside.dtype = bool
                    else:
                        anaInside += np.linalg.norm(np.asarray([H-localAnaH,K-localAnaK,L-localAnaL])*factor,axis=0)<dH

        if spurionType.lower() == 'both':
            mask = np.logical_or(monoInside,anaInside)
        elif spurionType.lower() == 'monochromator':
            mask = monoInside
        else:
            mask = anaInside
        if not maskInside:
            mask = np.logical_not(mask)
        return mask

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
            if self.fromNICOS: 
                attributes += ['x_gap','y_gap']
                values += ['monochromatorSlit'+x+'Gap' for x in ['X','Y']]
            
            for att,value in zip(attributes,values):
                val =  getattr(self,value)
                if not val.dtype == 'O':
                    dset = monoSlit.create_dataset(att,(1,),'float32')
                    dset[0] = val
                    dset.attrs['units'] = np.string_('mm')

        
        def addAna(self,inst):
            ana = inst.create_group('analyzer')
            ana.attrs['NX_class'] = np.string_('NXcrystal')
            
            attributes = ['d_spacing','nominal_energy','polar_angle','polar_angle_offset']+self.fromNICOS*['polar_angle_raw']
            values = ['analyzer'+x.replace('_',' ').title().replace(' ','') for x in attributes]
            units = ['anstrom','mev','degree','degree']+self.fromNICOS*['degree']


            for att,value,unit in zip(attributes,values,units):
                data = getattr(self,value)
                dset = ana.create_dataset(att,(len(data),),'float32')
                dset[:len(data)] = data
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
            if hasattr(self,'monochromatorRotationAngleZero'):
                v = self.monochromatorRotationAngleZero
            else:
                v = 0.0
            dset = mono.create_dataset('rotation_angle_zero',data=v,dtype='float32')
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
            if time[0] == np.array(None):
                time = [0.0]
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

            if hasattr(self,'singleDetector1'): # If the single detectors have been loaded
                for idx in ['1','8']:
                    segment = inst.create_group('segment_'+idx)
                    dset = segment.create_dataset('data',data=getattr(self,'singleDetector'+idx),dtype='int32')
                    dset.attrs['units']=np.string_('counts')
            
            
        
            attribute = ['a4offset','amplitude','background','boundaries','final_energy','width']
            for calibration,binning in zip(self.instrumentCalibrations,self.possibleBinnings):
                if binning is None: continue
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
def getScanParameter(self,f):

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

    scanItems = f.get('/entry/data/')
    for item in scanItems:
        if item in ['scanvar']: # only present in some files, but has to be ignored
            continue
        if not item in ['counts','summed_counts','en','h','intensity','k','l','monitor',
        'normalization','qx','qy','data'] and item[-4:]!='zero':
            scanParameters.append(item)
            fItem = f.get('/entry/data/{}'.format(item))
            scanUnits.append(decodeStr(fItem.attrs['units']))
            scanValues.append(np.array(fItem))
            try:
                scanDataPosition.append(decodeStr(fItem.attrs['target']))
            except:
                pass
    
    if len(scanParameters) == 0: # no data was stored... NICOS again, but due to manual scan
        
        scanItems = [x for x in np.array(f.get('/entry/scanvars'))[0].decode().split(',') if len(x)>0 ]
        ## In a manual scan, no scanvars are saved.... guess that it is either A3, A4 or Ei
        # Check first temperature and magnetic field


        if len(scanItems) == 0:
            try:
                if np.abs(np.diff(self.temperature)).mean()>0.01:
                    scanParameters = ['T']
                    scanValues = np.array([self.temperature])
                    scanUnits = ['K']
                    scanDataPosition = ['/entry/sample/temperature']
                elif np.abs(np.diff(self.magneticField)).mean()>0.01:
                    scanParameters = ['B']
                    scanValues = np.array([self.magneticField])
                    scanUnits = ['T']
                    scanDataPosition = ['/entry/sample/B']
                elif np.abs(np.diff(self.Ei)).mean()>0.001: # it is an energy scan!
                    scanParameters = ['Ei']
                    scanValues = np.array([self.Ei])
                    scanUnits = ['meV']
                    scanDataPosition = ['/entry/CAMEA/monochromator/energy']
                elif np.abs(np.diff(self.A3)).mean()>0.001:
                    scanParameters = ['A3']
                    scanValues = np.array([self.A3])
                    scanUnits = ['deg']
                    scanDataPosition = ['/entry/sample/rotation_angle']
                elif np.abs(np.diff(self.twotheta)).mean()>0.001:
                    scanParameters = ['A4']
                    scanValues = np.array([self.twotheta])
                    scanUnits = ['deg']
                    scanDataPosition = ['/entry/analyzer/polar_angle_raw']
                else:
                    pass
                    #raise AttributeError('Scan values from Datafile ("{}") cannot be determined'.format(self.name))
            except ValueError:
                pass
        else:
            for item in scanItems:
                if item in ['a3']:
                    item = item.upper()
                
                if item == 's2t':
                    item = 'A4'
                

                if item == 'A4':
                    fItem = getHDFInstrumentEntry(getInstrument(f),item,self.fromNICOS)
                elif item != 'CAMEA':
                    fItem = getHDFEntry(f,item,self.fromNICOS)
                    #fItem = f.get(position)    
                
                if item == 'CAMEA':
                    # We are doing a QScan
                    scanCommand = self.scanCommand.decode()

                    t = scanCommand.split('(')[0]
                    if not t.lower() == 'qcscan': continue
                    start = scanCommand.split('(')[2].split(')')[0]
                    step = scanCommand.split('(')[3].split(')')[0]
                    _,steps,monitor = scanCommand.split(')')[2].split(',')

                    centre = np.asarray([float(x) for x in start.split(',')])
                    step = np.asarray([float(x) for x in step.split(',')])
                    dist =np.linalg.norm(step)
                    step = _tools.LengthOrder(step)
                    dist*=1.0/np.linalg.norm(step)
                    steps = int(steps)
                    monitor = int(monitor.split('=')[1])
                    #constant = np.isclose(step,0)
                    textCentre = '['+', '.join([str(x) for x in centre])+'] '
                    textStep = _tools.generateLabelDirection(step,labels=['H','K','L','E'])
                    
                    scanParameters.append(textCentre+' + '+textStep)
                    scanUnits.append('RLU, meV')
                    scanValues.append(np.linspace(-dist*steps,dist*steps,steps*2+1))
                    
                else:
                    scanParameters.append(item)
                
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
    df.binning = 1
    
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



def shallowRead(files,parameters,fromNICOS=None):
    """Read a list of paramters from hdf file with minimal overhead
    
    Args:
        
        - files (list): List of files
        
        - parameters (list): List of parameters
        
    Returns:
        
        - list: Parameters in list after file, with asked properties in a sublist for each file
        
    """
    parameters = np.array(parameters)
    values = []
    possibleAttributes.sort(key=lambda v: v.lower())
    possible = []
    for p in parameters:
        possible.append(p in possibleAttributes)
    
    if not np.all(possible):
        if np.sum(np.logical_not(possible))>1:
            raise AttributeError('Parameters {} not found'.format(parameters[np.logical_not(possible)]))
        else:
            raise AttributeError('Parameter {} not found'.format(parameters[np.logical_not(possible)]))
    
    for file in files:
        vals = []
        if os.path.splitext(file)[-1] == '.hdf': # if an hdf file
            with hdf.File(file,mode='r') as f:
                if fromNICOS is None:
                    NICOS = checkNICOS(f)
                else:
                    NICOS = fromNICOS
                instr = getInstrument(f)
                for p in parameters:
                    if p == 'name':
                        v = os.path.basename(file)
                        vals.append(v)
                        continue
                    elif p == 'fileLocation':
                        v = os.path.dirname(file)
                        vals.append(v)
                        continue
                    elif p == 'twoTheta':
                        A4 = np.array(getHDFInstrumentEntry(instr,'A4',fromNICOS=NICOS))
                        for func,args in HDFInstrumentTranslationFunctions['A4']:
                            A4 = getattr(A4,func)(*args)
                        A4Offset = np.array(getHDFInstrumentEntry(instr,'A4Offset',fromNICOS=NICOS))
                        for func,args in HDFInstrumentTranslationFunctions['A4Offset']:
                            A4Offset = getattr(A4Offset,func)(*args)
                        vals.append(A4-A4Offset)
                        continue
                    elif p in HDFTranslation:
                        v = np.array(getHDFEntry(f,p,fromNICOS=NICOS))
                        TrF= HDFTranslationFunctions
                    elif p in HDFInstrumentTranslation:
                        v = np.array(getHDFInstrumentEntry(instr,p,fromNICOS=NICOS))
                        TrF= HDFInstrumentTranslationFunctions
                    else:
                        raise AttributeError('Parameter "{}" not found'.format(p))
                    for func,args in TrF[p]:
                        try:
                            v = getattr(v,func)(*args)
                        except IndexError as e:
                            v = 'Not In File'+func
                            break
                    
                    vals.append(v)
                values.append(vals)
        else:
            df = DataFile(file)
            for p in parameters:
                
                try:
                    vals.append(getattr(df,p))
                except:
                    if 'sample' in p:
                        try:
                            vals.append(getattr(df.sample,p.replace('sample','').lower()))
                        except:
                            vals.append('Not Found :S')
                    else:
                        try:
                            vals.append(getattr(df,p.lower()))
                        except:
                            vals.append('Not Found :S')
            values.append(vals)
            


    return values    





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
    I = PointerArray('I',files)
    Norm = PointerArray('Norm',files)
    Monitor = PointerArray('Monitor',files)
    energy = PointerArray('energy',files)
    
    
    a3 = []
    a4 = []
    a3Off = PointerArray('a3Off',files)
    a4Off = PointerArray('a4Off',files)
    Ei = PointerArray('Ei',files)
    instrumentCalibrationEf = PointerArray('instrumentCalibrationEf',files)
    instrumentCalibrationA4 = PointerArray('instrumentCalibrationA4',files)
    instrumentCalibrationEdges = PointerArray('instrumentCalibrationEdges',files)

    qx = PointerArray('qx',files)
    qy = PointerArray('qy',files)
    H = PointerArray('H',files)
    K = PointerArray('K',files)
    L = PointerArray('L',files)
    temperature = PointerArray('temperature',files)

    
    #mask = []
    
    scanParameters = []
    scanParamValue = []
    scanParamUnit = []
    #mask = []
    for datafile in files:
        #mask.append(datafile.mask)
        scanParameters.append(datafile.scanParameters)
        scanParamValue.append(datafile.scanValues)
        scanParamUnit.append(datafile.scanUnits)
            
        if np.array(datafile.A3Off).shape == ():
            datafile.A3Off = 0.0
        a3.append(datafile.A3-datafile.A3Off)
        if np.array(datafile.A4Off).shape == ():
            datafile.A4Off = [0.0]
        a4.append(datafile.A4-datafile.A4Off)
    
    if files[0].type=='nxs':
        return I,qx,qy,energy,Norm,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit,temperature,H,K,L
    else:
        return I,Monitor,a3,a3Off,a4,a4Off,instrumentCalibrationEf,\
        instrumentCalibrationA4,instrumentCalibrationEdges,Ei,scanParameters,scanParamValue,scanParamUnit,temperature

def assertFile(file):
    """Make sure that file exists for methods to work"""
    if not os.path.isfile(file):
        df = DataFile(file.replace('.nxs','.hdf'))
        con = df.convert(binning=8)
        con.saveNXsqom(file)

def checkNICOS(f):
    """Check if open hdf file is NICOS"""
    if not f.get('/entry/data/data') is None: # we have a NICOS file
        return True
    elif not f.get('/entry/data/counts') is None:
        return False
    else:
        return None # Corresponding to an error message

def isMultiFLEXX(fileLocation):
    with open(fileLocation,'r') as f:
        line = f.readline()
        while line:
            if line.find("### Instrument")!=-1:
                break
            line = f.readline()
        
        line = f.readline()[1:]
        line = line.strip().split(':')[0].split('_')[0]
    
    return line == 'mira'
