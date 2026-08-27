from MJOLNIR.Geometry.Instrument import Instrument,prediction
import MJOLNIR.Geometry.Analyser as Analyser
import MJOLNIR.Geometry.Detector as Detector
import MJOLNIR.Geometry.Wedge as Wedge
from MJOLNIR.Data import Sample
import pytest
import numpy as np
import warnings
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import os
import pytest

dataPath = 'samlpedata'

@pytest.mark.integration
def test_Instrument_init():
    Instr = Instrument()

    assert(np.all(Instr.position==(0,0,0)))

    Det = Detector.Detector(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.Analyser(position=(0.5,0,0),direction=(1,0,1))
    
    wedge = Wedge.Wedge(detectors=[Det,Det],analysers=Ana)

    Instr.wedges=[wedge,wedge]

    assert(Instr.settings['Initialized']==False)


@pytest.mark.integration
def test_Instrument_error():
    
    with pytest.raises(ValueError):
        Instr = Instrument(fileName='wrongDummyFile.bin')

    Instr = Instrument()

    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    with pytest.raises(AttributeError):
        Instr.wedges=Ana

    with pytest.raises(AttributeError):
        Instr.wedges=[Ana,Ana]

    with pytest.raises(AttributeError):
        Instr.append("Wrong object type")
    
    with pytest.raises(AttributeError):
        Instr.append(["List of",3.0,"wrong objects"])

    with pytest.raises(NotImplementedError):
        Instr.settings = {'Name','New dictionary'}

@pytest.mark.integration
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

@pytest.mark.integration
def test_Instrument_append():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.append([wedge,wedge])
    Instr.append(wedge)

    assert(len(Instr.wedges)==3)

@pytest.mark.integration
@pytest.mark.gui
def test_Instrument_plot():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append(wedge)
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    Instr.plot(ax)

@pytest.mark.unit
def test_Instrument_Setting(): 
    Instr = Instrument()
    Instr.settings['SettingVersion']=1.0
    assert(Instr.settings['SettingVersion']==1.0)

@pytest.mark.integration
def test_Instrument_Initialization():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0),concept='ManyToMany')
    pixels=33
    split = [12]
    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0),pixels=pixels,split=split)
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    

    wedge.append([Det,Det,Ana,Ana,Ana])

    with pytest.raises(ValueError):
        Instr.initialize()
    
    with pytest.raises(RuntimeError):
        print(Instr.A4)
    
    with pytest.raises(RuntimeError):
        print(Instr.Ef)
    
    Instr.append(wedge)
    with pytest.raises(ValueError):
        Instr.initialize()
    Instr.wedges[0].detectors[0].split = [0,12,20,pixels]
    Instr.initialize()

    assert(len(Instr.A4)==1)
    assert(len(Instr.A4[0])==2)
    assert(len(Instr.A4[0][0])==pixels)
    assert(len(Instr.A4)==len(Instr.Ef))
    assert(len(Instr.A4[0])==len(Instr.Ef[0]))
    assert(len(Instr.A4[0][0])==len(Instr.Ef[0][0]))
    assert(Instr.settings['Initialized']==True)

    with pytest.raises(NotImplementedError):
        Instr.A4 = []


    with pytest.raises(NotImplementedError):
        Instr.Ef = []


@pytest.mark.integration
def test_Instrument_saveload():
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


@pytest.mark.integration
@pytest.mark.data
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
        
    InstrLoaded = Instrument(fileName=tempFileName)
    os.remove(tempFileName)

    assert(Instr==InstrLoaded) 

@pytest.mark.integration
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

    with pytest.raises(ValueError):
        Instr = Instrument(fileName=temp_file)
        del Instr


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
    with pytest.raises(AttributeError):
        Instr = Instrument(fileName=temp_file)


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
    with pytest.raises(ValueError):
        Instr = Instrument(fileName=temp_file)

    os.remove(temp_file)

@pytest.mark.unit
def test_instrument_string_dummy(): # Todo: Improve test!
    Instr = Instrument()

    string = str(Instr)
    del string
    assert True

@pytest.mark.integration
def test_instrument_create_xml():

    Instr = Instrument()
    filename = 'temp'
    Instr.generateCAMEAXML(filename)

    Instr2 = Instrument(fileName=filename+'.xml')
    os.remove(filename+'.xml')
    assert(len(Instr2.wedges)==8)


@pytest.mark.unit
def test_Normalization_tables(quick):

    Instr = Instrument(fileName=os.path.join('Data','CAMEA_Updated.xml'))
    Instr.initialize()

    NF = os.path.join(dataPath,'camea2023n000083.hdf')
    #AF = 'TestData/1024/A4Normalization.h5'

    with pytest.raises(AttributeError):
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation=os.path.join(dataPath,''),plot=False,tables=[]) # No binning specified 

    with pytest.raises(AttributeError):
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation=os.path.join(dataPath,''),plot=False,tables=['Nothing?']) # Wrong binning

    if not quick==True:
        Instr.generateCalibration(Vanadiumdatafile=NF,  savelocation=os.path.join(dataPath,''),plot=False,tables=[1,3,8],sampleMass=4.7) 
    else:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation=os.path.join(dataPath,''),plot=False,tables=[1],sampleMass=4.7) 

@pytest.mark.integration
@pytest.mark.gui
def test_Prediction():
    A3Start = 0.0
    A3Stop = 100
    A3Steps = 101 
    Ei = 5.0
    A4 = [-36,-40]
    points = False
    # [H,K,L,A3,A4,0.0,0.0,Ei,Ef]
    HKL1 = np.array([1,0,0])
    HKL2 = np.array([0,0,1])
    A3R1 = 25.0
    A3R2 = 115.0
    #r1 = np.array([1,0,0,25.0,-24,0.0,0.0,Ei,Ei])
    #r2 = np.array([0,0,1,115.0,-24,0.0,0.0,Ei,Ei])
    
    cell = np.array([6.0,6.0,6.0,90.0,90.0,90.0])

    sample = Sample.calculateSample(cell,HKL1,HKL2,A3R1=A3R1,A3R2=A3R2, Ei=Ei,Ef=Ei)

    plt.ion()

    ax = prediction(A3Start=A3Start,A3Stop=A3Stop,A3Steps=A3Steps,A4Positions=A4,Ei=Ei,sample=sample,
    points=points, instrument='CAMEA')

    ax = prediction(A3Start=A3Start,A3Stop=A3Stop,A3Steps=A3Steps,A4Positions=A4,Ei=Ei,sample=sample,
    points=points, instrument='MultiFLEXX')

    ax = prediction(A3Start=A3Start,A3Stop=A3Stop,A3Steps=A3Steps,A4Positions=A4,Ei=Ei,sample=sample,
    points=points, instrument='Bambus')
    

