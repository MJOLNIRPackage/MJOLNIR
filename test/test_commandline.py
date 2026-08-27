from MJOLNIR.Data import DataFile
import subprocess
import os,sys
import pytest

dataFiles = [os.path.join('samlpedata',f) for f in ['camea2018n000136.hdf',
'camea2018n000178.hdf']]

returnText = [b'camea2018n000136.hdf: sc a3 0 da3 0.5 np 121 mn 150000	A3 scan around 1 0 0 YMnO3 T=10, 2T= -20\n',
b'camea2018n000178.hdf: sc a3 0 da3 1 np 181 mn 100000	PbTi T=1.5K Ei=5.5 2t=-10 HHL plane around 1 1 0\n'
]

if sys.platform == 'win32':
    returnText = [t[:-1] + b'\r\n' for t in returnText]


## Calibration inspector
@pytest.mark.commandline
@pytest.mark.unit
def test_CalibrationInspector_Help_Text():
    result = subprocess.check_output(['MJOLNIRCalibrationInspector','-h'], text=True)
    assert("usage: MJOLNIRCalibrationInspector [-h] [-s SAVE]" in result or "usage: MJOLNIRCalibrationInspector.exe [-h] [-s SAVE]" in result)
    assert("Inspection tool to visualize calibration tables in a data file." in result)

    
@pytest.mark.commandline
@pytest.mark.unit
@pytest.mark.gui
@pytest.mark.skip(reason="Fails on the headless Travis-ci")
def test_CalibrationInsepctor_Run():
    try:
      os.makedirs('_temp')
    except:
      pass
    subprocess.check_output(['MJOLNIRCalibrationInspector', '-s _temp/','-b 1',dataFiles[0]])
    # Creates 4 files in _temp
    filesCreated = ['_temp/Final_Energy_Individual_1.png',
    '_temp/Final_Energy_Overview_1.png',
    '_temp/Instrument_calibration_1.png',
    '_temp/Normalization_1.png']
    for file in filesCreated:
        assert(os.path.exists(file))
        os.remove(file)
    os.rmdir('_temp')
    


@pytest.mark.commandline
@pytest.mark.unit
def test_History_Help_Text():
    result = subprocess.check_output(['MJOLNIRHistory', '-h'], text=True)

    assert("usage: MJOLNIRHistory [-h] [-s SAVE] [-r] [DataFile [DataFile ...]]" in result or "usage: MJOLNIRHistory.exe [-h] [-s SAVE] [-r]" in result)
    assert("History tool for displaying files and command for selected data files." in result)


@pytest.mark.commandline
@pytest.mark.unit
@pytest.mark.data
def test_History_SingleFile():
    result = subprocess.check_output(['MJOLNIRHistory', dataFiles[0]])
    print(result)
    assert(returnText[0]==result)

@pytest.mark.commandline
@pytest.mark.unit
@pytest.mark.data
def test_History_MultipleFiles():
    call = ['MJOLNIRHistory'] + dataFiles
    
    result = subprocess.check_output(call)
    print(result)
    assert(b''.join(returnText)==result)

@pytest.mark.commandline
@pytest.mark.unit
@pytest.mark.data
def test_History_MultipleFiles_repeat():
    call = ['MJOLNIRHistory'] + dataFiles[-2:]
    results = subprocess.check_output(call)
    results2 = subprocess.check_output(['MJOLNIRHistory', '-r'])
    assert(results2 == results)
    
    
    
    

### 3DView
@pytest.mark.commandline
@pytest.mark.unit
def test_3DView_Help_Text():
    result = subprocess.check_output(['MJOLNIR3DView','-h'], text=True)
    assert("usage: MJOLNIR3DView [-h] [-r] [-b BINNING] [-d DQXDQYDE DQXDQYDE DQXDQYDE]" in result or "usage: MJOLNIR3DView.exe [-h] [-r] [-b BINNING]" in result)
    assert("Conversion tool for quick visualization using the viewer3D." in result)
    

  
@pytest.mark.skip(reason="Fails on the headless Travis-ci")    
@pytest.mark.commandline
@pytest.mark.unit
@pytest.mark.data
@pytest.mark.gui
def test_3DView_Run_Through():
    call = ['MJOLNIR3DView'] + [dataFiles[0],dataFiles[0]]
    subprocess.check_output(call)
    subprocess.check_output(['MJOLNIR3DView','-r','-b 1','-d', '0.1','0.1','0.2'])
    subprocess.check_output(['MJOLNIR3DView','-r','-b 1','-m 0','-M 1e-5'])
    assert(True)
    