from MJOLNIR.Data import DataFile
import subprocess
import os,sys
import pytest

dataFiles = [os.path.join('samlpedata',f) for f in ['camea2018n000136.hdf','camea2018n000136.nxs',
'camea2018n000137.nxs',
'camea2018n000178.hdf']]

returnText = [b'camea2018n000136.hdf: sc a3 0 da3 0.5 np 121 mn 150000	A3 scan around 1 0 0 YMnO3 T=10, 2T= -20\n',
b'camea2018n000136.nxs: sc a3 0 da3 0.5 np 121 mn 150000	A3 scan around 1 0 0 YMnO3 T=10, 2T= -20\n',
b'camea2018n000137.nxs: sc a3 0 da3 0.5 np 121 mn 150000	A3 scan around 1 0 0 YMnO3 T=10, 2T= -24\n',
b'camea2018n000178.hdf: sc a3 0 da3 1 np 181 mn 100000	PbTi T=1.5K Ei=5.5 2t=-10 HHL plane around 1 1 0\n'
]

if sys.platform == 'win32':
    returnText = [t[:-1] + b'\r\n' for t in returnText]


## Calibration inspector

def test_CalibrationInspector_Help_Text():
    result = subprocess.check_output(['MJOLNIRCalibrationInspector','-h'])
    print(result)
    helpTextWin32 = b'usage: MJOLNIRCalibrationInspector [-h] [-s SAVE]\r\n                                   [-p [PLOTLIST [PLOTLIST ...]]] [-b BINNING]\r\n                                   [DataFile]\r\n\r\nInspection tool to visualize calibration tables in a data file.\r\n\r\npositional arguments:\r\n  DataFile              Data file from which calibration table is to be\r\n                        plotted. If none provided file dialogue will appear.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -p [PLOTLIST [PLOTLIST ...]], --plot [PLOTLIST [PLOTLIST ...]]\r\n                        List of wanted plots to be generated. Should be\r\n                        "A4","Normalization","Ef","EfOverview". Default all of\r\n                        them.\r\n  -b BINNING, --binning BINNING\r\n                        Binning to be inspected. Default \'8\'\r\n'
    helpTextLinux = b'usage: MJOLNIRCalibrationInspector [-h] [-s SAVE]\n                                   [-p [PLOTLIST [PLOTLIST ...]]] [-b BINNING]\n                                   [DataFile]\n\nInspection tool to visualize calibration tables in a data file.\n\npositional arguments:\n  DataFile              Data file from which calibration table is to be\n                        plotted. If none provided file dialogue will appear.\n\noptional arguments:\n  -h, --help            show this help message and exit\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\n  -p [PLOTLIST [PLOTLIST ...]], --plot [PLOTLIST [PLOTLIST ...]]\n                        List of wanted plots to be generated. Should be\n                        "A4","Normalization","Ef","EfOverview". Default all of\n                        them.\n  -b BINNING, --binning BINNING\n                        Binning to be inspected. Default \'8\'\n'
    helpTextWin10 = b'usage: MJOLNIRCalibrationInspector.exe [-h] [-s SAVE]\r\n                                       [-p [PLOTLIST [PLOTLIST ...]]]\r\n                                       [-b BINNING]\r\n                                       [DataFile]\r\n\r\nInspection tool to visualize calibration tables in a data file.\r\n\r\npositional arguments:\r\n  DataFile              Data file from which calibration table is to be\r\n                        plotted. If none provided file dialogue will appear.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -p [PLOTLIST [PLOTLIST ...]], --plot [PLOTLIST [PLOTLIST ...]]\r\n                        List of wanted plots to be generated. Should be\r\n                        "A4","Normalization","Ef","EfOverview". Default all of\r\n                        them.\r\n  -b BINNING, --binning BINNING\r\n                        Binning to be inspected. Default \'8\'\r\n'
    helpTextPython39 = b'usage: MJOLNIRCalibrationInspector.exe [-h] [-s SAVE] [-p [PLOTLIST ...]]\r\n                                       [-b BINNING]\r\n                                       [DataFile]\r\n\r\nInspection tool to visualize calibration tables in a data file.\r\n\r\npositional arguments:\r\n  DataFile              Data file from which calibration table is to be\r\n                        plotted. If none provided file dialogue will appear.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -p [PLOTLIST ...], --plot [PLOTLIST ...]\r\n                        List of wanted plots to be generated. Should be\r\n                        "A4","Normalization","Ef","EfOverview". Default all of\r\n                        them.\r\n  -b BINNING, --binning BINNING\r\n                        Binning to be inspected. Default \'8\'\r\n'
    helpTextTravisPython39 = b'usage: MJOLNIRCalibrationInspector [-h] [-s SAVE] [-p [PLOTLIST ...]]\n                                   [-b BINNING]\n                                   [DataFile]\n\nInspection tool to visualize calibration tables in a data file.\n\npositional arguments:\n  DataFile              Data file from which calibration table is to be\n                        plotted. If none provided file dialogue will appear.\n\noptional arguments:\n  -h, --help            show this help message and exit\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\n  -p [PLOTLIST ...], --plot [PLOTLIST ...]\n                        List of wanted plots to be generated. Should be\n                        "A4","Normalization","Ef","EfOverview". Default all of\n                        them.\n  -b BINNING, --binning BINNING\n                        Binning to be inspected. Default \'8\'\n'

    helpTextPython39Windows =  b'usage: MJOLNIRCalibrationInspector [-h] [-s SAVE] [-p [PLOTLIST ...]]\r\n                                   [-b BINNING]\r\n                                   [DataFile]\r\n\r\nInspection tool to visualize calibration tables in a data file.\r\n\r\npositional arguments:\r\n  DataFile              Data file from which calibration table is to be\r\n                        plotted. If none provided file dialogue will appear.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -p [PLOTLIST ...], --plot [PLOTLIST ...]\r\n                        List of wanted plots to be generated. Should be\r\n                        "A4","Normalization","Ef","EfOverview". Default all of\r\n                        them.\r\n  -b BINNING, --binning BINNING\r\n                        Binning to be inspected. Default \'8\'\r\n'
    helpTextWin32Bool = result == helpTextWin32
    helpTextLinuxBool = result == helpTextLinux
    helpTextWin10Bool = result == helpTextWin10
    helpTextPython39Bool = result == helpTextPython39
    helpTextTravisPython39Bool = result == helpTextTravisPython39
    helpTextPython39WindowsBool = result == helpTextPython39Windows
    assert(helpTextWin32Bool or helpTextLinuxBool or helpTextWin10Bool or helpTextPython39Bool or helpTextTravisPython39Bool or helpTextPython39WindowsBool)


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
    




def test_History_Help_Text():
    result = subprocess.check_output(['MJOLNIRHistory', '-h'])
    print(result)
    helpTextWin32 = b'usage: MJOLNIRHistory [-h] [-s SAVE] [-r] [DataFile [DataFile ...]]\r\n\r\nHistory tool for displaying files and command for selected data files.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated history will be saved.\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n'
    helpTextLinux = b'usage: MJOLNIRHistory [-h] [-s SAVE] [-r] [DataFile [DataFile ...]]\n\nHistory tool for displaying files and command for selected data files.\n\npositional arguments:\n  DataFile              Data file(s) to be used. If none provided file\n                        dialogue will appear. Using string format, directory\n                        and year is also possible. See documentation.\n\noptional arguments:\n  -h, --help            show this help message and exit\n  -s SAVE, --save SAVE  Location to which the generated history will be saved.\n  -r, --reuse           Set flag to reuse files from previous usage. Default\n                        false.\n'
    helpTextWin10 = b'usage: MJOLNIRHistory.exe [-h] [-s SAVE] [-r] [DataFile [DataFile ...]]\r\n\r\nHistory tool for displaying files and command for selected data files.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated history will be saved.\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n'
    helpTextPython39Windows = b'usage: MJOLNIRHistory.exe [-h] [-s SAVE] [-r] [DataFile ...]\r\n\r\nHistory tool for displaying files and command for selected data files.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated history will be saved.\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n'
    helpTextTravisPython39Windows = helpTextPython39Windows.decode().replace('\r','').encode()
    helpTextTravisPython39 = helpTextTravisPython39Windows.decode().replace('.exe','').encode()
    helpTextTravisMoreTest = helpTextPython39Windows.decode().replace('.exe','').encode()

    helpTextWin32Bool = result == helpTextWin32
    helpTextLinuxBool = result == helpTextLinux
    helpTextWin10Bool = result == helpTextWin10
    helpTextPython39Bool = result == helpTextPython39Windows
    helpTextTravisPython39WindowsBool = result == helpTextTravisPython39Windows
    helpTextTravisPython39Bool = result == helpTextTravisPython39
    helpTextTravisMoreTestBool = result == helpTextTravisMoreTest
    assert(helpTextWin32Bool or helpTextLinuxBool or helpTextWin10Bool or helpTextPython39Bool or helpTextTravisPython39WindowsBool or helpTextTravisPython39Bool or helpTextTravisMoreTestBool)


def test_History_SingleFile():
    result = subprocess.check_output(['MJOLNIRHistory', dataFiles[0]])
    print(result)
    assert(returnText[0]==result)


def test_History_MultipleFiles():
    call = ['MJOLNIRHistory'] + dataFiles
    
    result = subprocess.check_output(call)
    print(result)
    assert(b''.join(returnText)==result)

def test_History_MultipleFiles_repeat():
    call = ['MJOLNIRHistory'] + dataFiles[-2:]
    results = subprocess.check_output(call)
    results2 = subprocess.check_output(['MJOLNIRHistory', '-r'])
    assert(results2 == results)
    
    
    
    
### Converter ###


def test_Convert_Help_Text():
    result = subprocess.check_output(['MJOLNIRConvert', '-h'])
    print(result)
    helpTextWin32 = b"usage: MJOLNIRConvert [-h] [-s SAVE] [-b BINNING] [-r]\r\n                      [DataFile [DataFile ...]]\r\n\r\nConversion tool for converting output h5 files to nxs files.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -b BINNING, --binning BINNING\r\n                        Binning performed. Default '8'\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n"
    helpTextLinux = b"usage: MJOLNIRConvert [-h] [-s SAVE] [-b BINNING] [-r]\n                      [DataFile [DataFile ...]]\n\nConversion tool for converting output h5 files to nxs files.\n\npositional arguments:\n  DataFile              Data file(s) to be used. If none provided file\n                        dialogue will appear. Using string format, directory\n                        and year is also possible. See documentation.\n\noptional arguments:\n  -h, --help            show this help message and exit\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\n  -b BINNING, --binning BINNING\n                        Binning performed. Default '8'\n  -r, --reuse           Set flag to reuse files from previous usage. Default\n                        false.\n"
    helpTextWin10 = b'usage: MJOLNIRConvert.exe [-h] [-s SAVE] [-b BINNING] [-r]\r\n                          [DataFile [DataFile ...]]\r\n\r\nConversion tool for converting output h5 files to nxs files.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -b BINNING, --binning BINNING\r\n                        Binning performed. Default \'8\'\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n'
    helpTextPython39Windows = b"usage: MJOLNIRConvert.exe [-h] [-s SAVE] [-b BINNING] [-r] [DataFile ...]\r\n\r\nConversion tool for converting output h5 files to nxs files.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -s SAVE, --save SAVE  Location to which the generated file will be saved.\r\n  -b BINNING, --binning BINNING\r\n                        Binning performed. Default '8'\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n"
    helpTextTravisPython39Windows = helpTextPython39Windows.decode().replace('\r','').encode()
    helpTextTravisPython39 = helpTextTravisPython39Windows.decode().replace('.exe','').encode()
    helpTextTravisMoreTest = helpTextPython39Windows.decode().replace('.exe','').encode()

    helpTextWin32Bool = result == helpTextWin32
    helpTextLinuxBool = result == helpTextLinux
    helpTextWin10Bool = result == helpTextWin10
    helpTextPython39Bool = result == helpTextPython39Windows
    helpTextTravisPython39WindowsBool = result == helpTextTravisPython39Windows
    helpTextTravisPython39Bool = result == helpTextTravisPython39
    helpTextTravisMoreTestBool = result == helpTextTravisMoreTest
    assert(helpTextWin32Bool or helpTextLinuxBool or helpTextWin10Bool or helpTextPython39Bool or helpTextTravisPython39WindowsBool or helpTextTravisPython39Bool or helpTextTravisMoreTestBool)


    
def test_Convert_binning():
    subprocess.check_output(['MJOLNIRConvert',dataFiles[0]])
    f = DataFile.DataFile(dataFiles[0].replace('hdf','nxs'))
    assert(f.binning == 8)
    subprocess.check_output(['MJOLNIRConvert',dataFiles[0],'-b 1'])
    f = DataFile.DataFile(dataFiles[0].replace('hdf','nxs'))
    assert(f.name == dataFiles[0].split(os.path.sep)[1].replace('hdf','nxs'))
    assert(f.binning == 1)
    
def test_Convert_SaveLocation():
    if not os.path.isdir(os.path.join('Data','Data')):
          os.makedirs(os.path.join('Data','Data'))
    subprocess.check_output(['MJOLNIRConvert',dataFiles[0],'-s Data'+os.path.sep])
    fileName = dataFiles[0].split(os.path.sep)[1].replace('hdf','nxs')
    f = DataFile.DataFile('Data'+os.path.sep+fileName)
    
    assert(f.binning == 8)
    print(fileName)
    print(f.name)
    assert(f.name == fileName)
    
    
def test_Convert_Reuse():
    subSet = [dataFiles[0],dataFiles[3]]
    call = ['MJOLNIRConvert'] + subSet
    subprocess.check_output(call)
    print('calling: "'+' '.join(call)+'"')
    for file in subSet:
        os.remove(file.replace('hdf','nxs'))
    subprocess.check_output(['MJOLNIRConvert','-r'])
    for file in subSet:
        assert(os.path.isfile(file.replace('hdf','nxs')))


### 3DView
def test_3DView_Help_Text():
    result = subprocess.check_output(['MJOLNIR3DView','-h'])
    print(result)
    helpTextWin32 = b"usage: MJOLNIR3DView [-h] [-r] [-b BINNING] [-d DQXDQYDE DQXDQYDE DQXDQYDE]\r\n                     [-M VMAX] [-m VMIN]\r\n                     [DataFile [DataFile ...]]\r\n\r\nConversion tool for quick visualization using the viewer3D.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n  -b BINNING, --binning BINNING\r\n                        Binning performed. Default '8'\r\n  -d DQXDQYDE DQXDQYDE DQXDQYDE, --dQxdQydE DQXDQYDE DQXDQYDE DQXDQYDE\r\n                        Binning used to plot in 3D, Default [0.03,0.03,0.08]\r\n  -M VMAX, --VMax VMAX  Maximal value for plotting, default max of data\r\n  -m VMIN, --VMin VMIN  Minimal value for plotting, default min of data\r\n"
    helpTextLinux = b"usage: MJOLNIR3DView [-h] [-r] [-b BINNING] [-d DQXDQYDE DQXDQYDE DQXDQYDE]\n                     [-M VMAX] [-m VMIN]\n                     [DataFile [DataFile ...]]\n\nConversion tool for quick visualization using the viewer3D.\n\npositional arguments:\n  DataFile              Data file(s) to be used. If none provided file\n                        dialogue will appear. Using string format, directory\n                        and year is also possible. See documentation.\n\noptional arguments:\n  -h, --help            show this help message and exit\n  -r, --reuse           Set flag to reuse files from previous usage. Default\n                        false.\n  -b BINNING, --binning BINNING\n                        Binning performed. Default \'8\'\n  -d DQXDQYDE DQXDQYDE DQXDQYDE, --dQxdQydE DQXDQYDE DQXDQYDE DQXDQYDE\n                        Binning used to plot in 3D, Default [0.03,0.03,0.08]\n  -M VMAX, --VMax VMAX  Maximal value for plotting, default max of data\n  -m VMIN, --VMin VMIN  Minimal value for plotting, default min of data\n"
    helpTextWin10 = b'usage: MJOLNIR3DView.exe [-h] [-r] [-b BINNING]\r\n                         [-d DQXDQYDE DQXDQYDE DQXDQYDE] [-M VMAX] [-m VMIN]\r\n                         [DataFile [DataFile ...]]\r\n\r\nConversion tool for quick visualization using the viewer3D.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n  -b BINNING, --binning BINNING\r\n                        Binning performed. Default \'8\'\r\n  -d DQXDQYDE DQXDQYDE DQXDQYDE, --dQxdQydE DQXDQYDE DQXDQYDE DQXDQYDE\r\n                        Binning used to plot in 3D, Default [0.03,0.03,0.08]\r\n  -M VMAX, --VMax VMAX  Maximal value for plotting, default max of data\r\n  -m VMIN, --VMin VMIN  Minimal value for plotting, default min of data\r\n'
    helpTextPython39 = b"usage: MJOLNIR3DView [-h] [-r] [-b BINNING] [-d DQXDQYDE DQXDQYDE DQXDQYDE]\n                     [-M VMAX] [-m VMIN]\n                     [DataFile ...]\n\nConversion tool for quick visualization using the viewer3D.\n\npositional arguments:\n  DataFile              Data file(s) to be used. If none provided file\n                        dialogue will appear. Using string format, directory\n                        and year is also possible. See documentation.\n\noptional arguments:\n  -h, --help            show this help message and exit\n  -r, --reuse           Set flag to reuse files from previous usage. Default\n                        false.\n  -b BINNING, --binning BINNING\n                        Binning performed. Default '8'\n  -d DQXDQYDE DQXDQYDE DQXDQYDE, --dQxdQydE DQXDQYDE DQXDQYDE DQXDQYDE\n                        Binning used to plot in 3D, Default [0.03,0.03,0.08]\n  -M VMAX, --VMax VMAX  Maximal value for plotting, default max of data\n  -m VMIN, --VMin VMIN  Minimal value for plotting, default min of data\n"
    helpTextTravisPython39 = helpTextPython39.decode().replace('\r','').encode()
    helpTextTravisMoreTest = b"usage: MJOLNIR3DView [-h] [-r] [-b BINNING] [-d DQXDQYDE DQXDQYDE DQXDQYDE]\r\n                     [-M VMAX] [-m VMIN]\r\n                     [DataFile ...]\r\n\r\nConversion tool for quick visualization using the viewer3D.\r\n\r\npositional arguments:\r\n  DataFile              Data file(s) to be used. If none provided file\r\n                        dialogue will appear. Using string format, directory\r\n                        and year is also possible. See documentation.\r\n\r\noptional arguments:\r\n  -h, --help            show this help message and exit\r\n  -r, --reuse           Set flag to reuse files from previous usage. Default\r\n                        false.\r\n  -b BINNING, --binning BINNING\r\n                        Binning performed. Default '8'\r\n  -d DQXDQYDE DQXDQYDE DQXDQYDE, --dQxdQydE DQXDQYDE DQXDQYDE DQXDQYDE\r\n                        Binning used to plot in 3D, Default [0.03,0.03,0.08]\r\n  -M VMAX, --VMax VMAX  Maximal value for plotting, default max of data\r\n  -m VMIN, --VMin VMIN  Minimal value for plotting, default min of data\r\n"


    helpTextWin32Bool = result == helpTextWin32
    helpTextLinuxBool = result == helpTextLinux
    helpTextWin10Bool = result == helpTextWin10
    helpTextPython39Bool = result == helpTextPython39
    helpTextTravisPython39Bool = result == helpTextTravisPython39
    helpTextTravisMoreTestBool = result == helpTextTravisMoreTest
    assert(helpTextWin32Bool or helpTextLinuxBool or helpTextWin10Bool or helpTextPython39Bool or helpTextTravisPython39Bool or helpTextTravisMoreTestBool)


  
@pytest.mark.skip(reason="Fails on the headless Travis-ci")    
def test_3DView_Run_Through():
    call = ['MJOLNIR3DView'] + [dataFiles[0],dataFiles[0]]
    subprocess.check_output(call)
    subprocess.check_output(['MJOLNIR3DView','-r','-b 1','-d', '0.1','0.1','0.2'])
    subprocess.check_output(['MJOLNIR3DView','-r','-b 1','-m 0','-M 1e-5'])
    assert(True)
    