from MJOLNIR.Geometry.Instrument import Instrument
import MJOLNIR.Geometry.Analyser as Analyser
import MJOLNIR.Geometry.Detector as Detector
import MJOLNIR.Geometry.Wedge as Wedge
import pytest
import numpy as np
import warnings
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import os

dataPath = 'Data'

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
        Instr = Instrument(fileName='wrongDummyFile.bin')
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
    Instr.wedges[0].detectors[0].split = [0,12,20,pixels]
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
        Instr = Instrument(fileName=temp_file)
        del Instr
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
        Instr = Instrument(fileName=temp_file)
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
        Instr = Instrument(fileName=temp_file)
        assert False
    except ValueError:
        assert True
    os.remove(temp_file)

def test_instrument_string_dummy(): # Todo: Improve test!
    Instr = Instrument()

    string = str(Instr)
    del string
    assert True
    
def test_instrument_create_xml():

    Instr = Instrument()
    filename = 'temp'
    Instr.generateCAMEAXML(filename)

    Instr2 = Instrument(fileName=filename+'.xml')
    os.remove(filename+'.xml')
    assert(len(Instr2.wedges)==8)


@pytest.mark.unit
def test_Normalization_tables(quick):

    Instr = Instrument(fileName=os.path.join(dataPath,'CAMEA_Updated.xml'))
    Instr.initialize()

    NF = os.path.join(dataPath,'camea2018n000038.hdf')
    #AF = 'TestData/1024/A4Normalization.h5'

    try:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation=os.path.join(dataPath,''),plot=False,tables=[]) # No binning specified 
        assert False
    except AttributeError:
        assert True

    try:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation=os.path.join(dataPath,''),plot=False,tables=['Nothing?']) # Wrong binning
        assert False
    except AttributeError:
        assert True

    if not quick==True:
        Instr.generateCalibration(Vanadiumdatafile=NF,  savelocation=os.path.join(dataPath,''),plot=False,tables=[1,3,8]) 
    else:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation=os.path.join(dataPath,''),plot=False,tables=[1]) 


