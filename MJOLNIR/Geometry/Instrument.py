import sys
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import numpy as np
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector,Wedge
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

class Instrument(GeometryConcept.GeometryConcept):
    def __init__(self, position=(0,0,0),wedges=[],filename='',**kwargs):
        """Instrument object used to calculated analytic scattering coverage. 
        Based on the GeometryConcept object it contains all needed information about the setup used in further calculations.

        Kwargs:
            
            - position (float 3d): Position of the instrument alwasy at origin(?) (default (0,0,0))

            - wedges (list of wedges or single wedge): Wedge or list of wedges which the instrument consists of (default empty)

            - filename (string): Filename of xml file (ending in xml). To load binary files use self.load(filename).

        Raises:
            
            - AttributeError
        
        """
        
        self._wedges = []

        
        self._settings = {}
        if filename !='':
            if(filename.split('.')[-1]=='xml'):
                parseXML(self,filename)
            else:
                raise ValueError('File not of type XML.')
        else:
            super(Instrument,self).__init__(position)
            
            for key in kwargs:
                self.settings[key]=kwargs[key]
            self._settings['Initialized']=False
            self.append(wedges)

    @property
    def wedges(self):
        return self._wedges

    @wedges.getter
    def wedges(self):
        return self._wedges

    @wedges.setter
    def wedges(self,wedges):
        if len(self.wedges)!=0:
            warnings.warn('The list of wedges is not empty! Appending new wedges(s)')
        if isinstance(wedges, list):
            for ana in wedges:
                if not issubclass(type(ana),Wedge.Wedge):
                    raise AttributeError('Object is not an wedge or a simple list of these')
                self._wedges.append(ana)
        else:
            if not issubclass(type(wedges),Wedge.Wedge):
                raise AttributeError('Object is not an analyser or a simple list of these')
            self._wedges.append(wedges)
    
    def append(self,wedge):
        """Append wedge(s) to instrument.

        Args
            
            - wedge (Wedge(s)): Single wedge or list of wedges
        """
        if isinstance(wedge,list):
            for obj in wedge:
                if issubclass(type(obj),Wedge.Wedge):
                    self._wedges.append(obj)
                else:
                    raise AttributeError('Object not wedge or a simple list of wedges')
        else:
            if issubclass(type(wedge),Wedge.Wedge):
                    self._wedges.append(wedge)
            else:
                raise AttributeError('Object not wedge or a simple list of wedges')

    def plot(self,ax):
        """Recursive plotting routine."""
        for wedge in self.wedges:
            wedge.plot(ax,offset=self.position)

    @property
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')

    def __str__(self):
        string = '{} with settings:\n'.format(self.__class__)
        for attrib in self.settings:
            string+='{}:\t{}\n'.format(attrib,self.settings[attrib])
        string+='\n'    
        string+='Containing the following wedges:\n'
        for wedge in self.wedges:
            string+=str(wedge)+'\n'
        return string

    def initialize(self):
        """Method to initialize and perform analytical calulations of scattering quantities. 
        Initializes:

            -  A4: Matrix holding pixel A4. Shape (len(Wedges),len(detectors),pixels)
            
            -  Ef: Matrix holding pixel Ef. Shape (len(Wedges),len(detectors),pixels)
        """
        factorLambdasqrtE = 9.0445678

        if(len(self.wedges)==0):
            raise ValueError('Instrument does not contain any wedges and can thus not be initialized.')
        self._A4 = []
        self._Ef = []
        beamDirection = np.array([0.0,1.0,0.0])

        for wedge in self.wedges:
            detectorPixelPositions,analyserPixelPositions = wedge.calculateDetectorAnalyserPositions()

            A4 = [np.arccos(np.divide(np.dot(AnalyserPos,beamDirection),
                np.linalg.norm(AnalyserPos,axis=1))) for AnalyserPos in analyserPixelPositions]
            relPos = [detectorPixelPositions[i]-analyserPixelPositions[i] for i in range(len(analyserPixelPositions))]

            A6 = [np.arccos(np.divide(np.einsum('ij,ij->i',analyserPixelPositions[i],relPos[i]),
                np.linalg.norm(analyserPixelPositions[i],axis=1)*np.linalg.norm(relPos[i],axis=1))) for i in range(len(analyserPixelPositions))]
            Ef = [np.power(factorLambdasqrtE/(wedge.analysers[0].d_spacing*2.0*np.sin(A6Sub/2.0)),2.0) for A6Sub in A6] ## TODO: d_spacing Make generic
            self._A4.append(A4)
            self._Ef.append(Ef)





        self.settings['Initialized']=True

    @property
    def A4(self):
        return self._A4
    @A4.getter
    def A4(self):
        if(self.settings['Initialized']==False):
            raise RuntimeError('Instrument is not initialized.')
        return self._A4
    @A4.setter
    def A4(self,*args,**kwargs):
        raise NotImplementedError('A4 cannot be overwritten.')

    @property
    def Ef(self):
        return self._Ef
    @Ef.getter
    def Ef(self):
        if(self.settings['Initialized']==False):
            raise RuntimeError('Instrument is not initialized.')
        return self._Ef
    @Ef.setter
    def Ef(self,*args,**kwargs):
        raise NotImplementedError('Ef cannot be overwritten.')

    
    def saveXML(self,filename):
        """Method for saving current file as XML in filename."""
        XMLString = '<?xml version="1.0"?>\n'
        XMLString+= '<Instrument '
        for attrib in self.settings:
            XMLString+="{}='{}' ".format(attrib,self.settings[attrib])

        XMLString+="position='"+','.join([str(x) for x in self.position])+"'"
        XMLString+='>\n'
            
        for wedge in self.wedges:
            XMLString+="\t<Wedge "
            for attrib in wedge.settings:
                XMLString+="{}='{}' ".format(attrib,wedge.settings[attrib])

            XMLString+="position='"+','.join([str(x) for x in wedge.position])+"'"
            XMLString+='>\n'

            for item in wedge.analysers + wedge.detectors:
                itemClass = str(item.__class__).split('.')[-1][:-2]
                XMLString+="\t\t<{}".format(itemClass)
                for key in item.__dict__:
                    value = item.__getattribute__(key)
                    if isinstance(value,type(np.array([0,0,0]))):
                        valueStr = ','.join([str(x) for x in item.__getattribute__(key)])
                    else:
                        valueStr = str(value)
                    XMLString+=" {}='{}'".format(str(key)[1:],valueStr)
                XMLString+="></{}>\n".format(itemClass)
                
            
            XMLString+="\t</Wedge>\n"
        XMLString+="</Instrument>\n"
    
        f = open(filename,'w')
        f.write(XMLString)
        f.close()

def parseXML(Instr,filename):
    import xml.etree.ElementTree as ET


    tree = ET.parse(filename)
    instr_root = tree.getroot()

    
    
        
    for attrib in instr_root.keys():
        if attrib=='position':
            Instr.position = np.array(instr_root.attrib[attrib].split(','))
        Instr.settings[attrib]=instr_root.attrib[attrib]
    
    
    Instr._wedges=[]
    for wedge in instr_root.getchildren():
        
        if wedge.tag in dir(Wedge):
            Wedgeclass_ = getattr(Wedge, wedge.tag)
        else:
            raise ValueError("Element is supposed to be a Wedge, but got '{}'.".format(wedge.tag))
        wedgeSettings = {}
        
        for attrib in wedge.keys():
            if attrib=='concept':
                wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=str)
            else:        
                wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=float)
            
        temp_wedge = Wedgeclass_(**wedgeSettings)
        
        
        
        
        for item in wedge.getchildren():
            if item.tag in dir(Detector):
                class_ = getattr(Detector, item.tag)
            elif item.tag in dir(Analyser):
                class_ = getattr(Analyser,item.tag)
            else:
                raise ValueError("Item '{}' not recognized as MJOLNIR detector or analyser.".format(item.tag))
            
            itemSettings = {}
            for attrib in item.keys():
                attribVal = item.get(attrib).strip().split(',')
                if len(attribVal)==1:
                    if(attrib=='split'):
                        try:
                            val=float(attribVal[0])
                        except ValueError:
                            val=[]
                        itemSettings[attrib]=val
                    else:
                        itemSettings[attrib]=float(attribVal[0])
                else:
                    if(attrib=='split'):
                        #print(type(attribVal))
                        itemSettings[attrib]=attribVal
                    else:
                        itemSettings[attrib]=np.array(attribVal,dtype=float)    
            try:
                temp_item = class_(**itemSettings)
            except TypeError as e:
                print(e.args[0])
                raise ValueError('Item {} misses argument(s):{}'.format(class_,e.args[0].split(':')[0]))
            except AttributeError as e:
                raise AttributeError('Error in passing {} with attributes {}'.format(class_,itemSettings))
            except ValueError:
                raise ValueError('Item {} not initialized due to error.'.format(class_))
            #print(temp_item)
            temp_wedge.append(temp_item)
            #print()

        #print(str(temp_wedge))
        Instr.append(temp_wedge)
    #print(str(Instr))
   


def test_Instrument_init():
    Instr = Instrument()

    assert(np.all(Instr.position==(0,0,0)))

    Det = Detector.Detector(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.Analyser(position=(0.5,0,0),direction=(1,0,1))
    
    wedge = Wedge.Wedge(detectors=[Det,Det],analysers=Ana)

    Instr.wedges=[wedge,wedge]

    assert(Instr.settings['Initialized']==False)



def test_Instrument_error():
    
    try:
        Instr = Instrument(filename='wrongDummyFile.bin')
        assert False
    except ValueError:
        assert True

    Instr = Instrument()

    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    try:
        Instr.wedges=Ana
        assert False
    except AttributeError:
        assert True

    try:
        Instr.wedges=[Ana,Ana]
        assert False
    except AttributeError:
        assert True

    try:
        Instr.append("Wrong object type")
        assert False
    except AttributeError:
        assert True
    
    try:
        Instr.append(["List of",3.0,"wrong objects"])
        assert False
    except AttributeError:
        assert True

    try:
        Instr.settings = {'Name','New dictionary'}
        assert False
    except NotImplementedError:
        return True


def test_Instrument_warnings():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.wedges = wedge

    with warnings.catch_warnings(record=True) as w: # From https://docs.python.org/3.1/library/warnings.html
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        Instr.wedges = wedge
        # Verify some things
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert 'The list of wedges is not empty! Appending new wedges(s)' in str(w[0].message)


def test_Instrument_append():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.append([wedge,wedge])
    Instr.append(wedge)

    assert(len(Instr.wedges)==3)


def test_Instrument_plot():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append(wedge)
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    Instr.plot(ax)

def test_Instrument_Setting(): 
    Instr = Instrument()
    Instr.settings['SettingVersion']=1.0
    assert(Instr.settings['SettingVersion']==1.0)


def test_Instrument_Initialization():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0),concept='ManyToMany')
    pixels=33
    split = [12]
    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0),pixels=pixels,split=split)
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    

    wedge.append([Det,Det,Ana,Ana,Ana])

    try:
        Instr.initialize()
        assert False
    except ValueError:
        assert True

    try:
        print(Instr.A4)
        assert False
    except RuntimeError:
        assert True

    try:
        print(Instr.Ef)
        assert False
    except RuntimeError:
        assert True

    Instr.append(wedge)
    try:
        Instr.initialize()
        assert False
    except ValueError:
        assert True
    Instr.wedges[0].detectors[0].split = [12,20]
    Instr.initialize()

    assert(len(Instr.A4)==1)
    assert(len(Instr.A4[0])==2)
    assert(len(Instr.A4[0][0])==pixels)
    assert(len(Instr.A4)==len(Instr.Ef))
    assert(len(Instr.A4[0])==len(Instr.Ef[0]))
    assert(len(Instr.A4[0][0])==len(Instr.Ef[0][0]))
    assert(Instr.settings['Initialized']==True)

    try:
        Instr.A4 = []
        assert False
    except NotImplementedError:
        assert True

    try:
        Instr.Ef = []
        assert False
    except NotImplementedError:
        assert True



def test_Istrument_saveload():
    import os
    Instr = Instrument(position=(0,1,0))
    Instr2 = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append(wedge)

    tempFile = 'temp.bin'
    Instr.save(tempFile)
    Instr2.load(tempFile)
    os.remove(tempFile)
    

    assert(Instr==Instr2)



def test_parseXML(): # Improve this test!

    tempFileName = '__temp__.xml'
        
    Instr = Instrument()
    Instr.settings['Author'] = 'Jakob Lass'

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append([wedge,wedge])
    Instr.append(wedge)
    Instr.saveXML(tempFileName)
        
    InstrLoaded = Instrument(filename=tempFileName)
    os.remove(tempFileName)

    assert(Instr==InstrLoaded)


def test_XML_errors():

    fileString = ""
    fileString+="<?xml version='1.0'?>"
    fileString+="<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>"
    fileString+="<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>"
    fileString+="<FlatAnalyser direction='0.707,0.0,0.707' d_spacing='3.35' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>"
    fileString+="<TubeDetector1D position='1.198,0.0580,0.71' direction='0.998,0.04841,0.0' pixels='456' length='0.883' diameter='0.02' split='57, 114, 171, 228, 285, 342, 399'></TubeDetector1D>"
    fileString+="</Wedge>"
    fileString+="</Instrument>"

    temp_file = 'Tempfile.xml'
    f = open(temp_file,'w')
    f.write(fileString)
    f.close()

    try:
        Instr = Instrument(filename=temp_file)
        assert False
    except ValueError:
        assert True

    fileString = ""
    fileString+="<?xml version='1.0'?>"
    fileString+="<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>"
    fileString+="<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>"
    fileString+="<FlatAnalyser position='0.0580,0.71' direction='0.707,0.0,0.707' d_spacing='3.35' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>"
    fileString+="<TubeDetector1D position='1.198,0.0580,0.71' direction='0.998,0.04841,0.0' pixels='456' length='0.883' diameter='0.02' split='57, 114, 171, 228, 285, 342, 399'></TubeDetector1D>"
    fileString+="</Wedge>"
    fileString+="</Instrument>"
    f = open(temp_file,'w')
    f.write(fileString)
    f.close()
    try:
        Instr = Instrument(filename=temp_file)
        assert False
    except AttributeError:
        assert True

    fileString = ""
    fileString+="<?xml version='1.0'?>"
    fileString+="<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>"
    fileString+="<FlatAnalyser position='0.0,0.0,0.0' concept='ManyToMany'>"
    fileString+="<FlatAnalyser position='0.0580,0.71' direction='0.707,0.0,0.707' d_spacing='3.35' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>"
    fileString+="<TubeDetector1D position='1.198,0.0580,0.71' direction='0.998,0.04841,0.0' pixels='456' length='0.883' diameter='0.02' split='57, 114, 171, 228, 285, 342, 399'></TubeDetector1D>"
    fileString+="</FlatAnalyser>"
    fileString+="</Instrument>"
    f = open(temp_file,'w')
    f.write(fileString)
    f.close()
    try:
        Instr = Instrument(filename=temp_file)
        assert False
    except ValueError:
        assert True
    os.remove(temp_file)

def test_instrument_string_dummy(): # Todo: Improve test!
    Instr = Instrument()

    string = str(Instr)
    assert True
