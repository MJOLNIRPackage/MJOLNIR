import numpy as np
from MJOLNIR.Data.DataFile import DataFile,decodeStr,createEmptyDataFile,assertFile
from MJOLNIR import _tools
import MJOLNIR.Data.Sample
import matplotlib.pyplot as plt
import os

dataPath = 'Data'


def test_DataFile():
    try:
        DF = DataFile('/nope.txt')
        assert False
    except:
        assert True

    try:
        DF= DataFile(os.path.join(dataPath,'CAMEA_Full.xml')) # Wrong file
        assert False
    except:
        assert True

    files = [os.path.join(dataPath,'camea2018n000137.hdf'),
             os.path.join(dataPath,'camea2018n000137.nxs')]
    DF1 = DataFile(files[0])
    assertFile(files[1])
    DF2 = DataFile(files[1])
    s = str(DF2)
    sampleS = str(DF2.sample)
    print(str(DF1.sample))
    print(str(DF2.sample))
    assert(DF1.sample == DF2.sample)

def test_DataFile_equility():
    f1 = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    print('----------')
    f2 = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    assert(f1==f2)
    print('----------')
    f3 = DataFile(f2)
    assert(f1==f3)
    print('----------')
    f4 = DataFile(os.path.join(dataPath,'camea2018n000137.hdf'))
    assert(f1!=f4)
    f5 = f2.convert(binning=8)
    f5.saveNXsqom(os.path.join(dataPath,'camea2018n000136.nxs'))
    f3 = DataFile(os.path.join(dataPath,'camea2018n000136.nxs'))
    assert(f1==f3)
    print('----------')

def test_DataFile_plotA4():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2018n000136.hdf')
    fileName2= os.path.join(dataPath,'camea2018n000136.nxs')
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
    import matplotlib
    matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2018n000136.hdf')
    fileName2= os.path.join(dataPath,'camea2018n000136.nxs')
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
    import matplotlib
    matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2018n000136.hdf')
    fileName2= os.path.join(dataPath,'camea2018n000136.nxs')
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
    import matplotlib
    matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2018n000136.hdf')
    fileName2= os.path.join(dataPath,'camea2018n000136.nxs')
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

    files = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000136.nxs')]
    assertFile(files[1])
    for file in files:
        dfile = DataFile(file)
        assert(dfile.scanParameters[0]=='rotation_angle')
        assert(len(dfile.scanParameters)==len(dfile.scanUnits))
        assert(len(dfile.scanParameters)==len(dfile.scanValues))
        assert(len(dfile.scanParameters)==1)
        assert(dfile.scanUnits[0]=='degree')
        ##assert(np.all(dfile.scanValues==np.arange(0,150,1)))


def test_DataFile_Error():
    df = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))

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
    
    df2 = DataFile(os.path.join(dataPath,'camea2018n000136.nxs'))
    try:
        df.saveNXsqom(os.path.join(dataPath,'saving.nxs'))
        assert False
    except AttributeError: # File does not have original_file attribute
        assert True

    df2.type = 'WrongType'
    try:
        df2.saveNXsqom(os.path.join(dataPath,'saving.nxs'))
        assert False
    except AttributeError: # Manually change type to wrong
        assert True


def test_DataFile_SaveLoad():
    df = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    df.saveHDF(os.path.join(dataPath,'camea2018n000136_2.hdf'))
    df2= DataFile(os.path.join(dataPath,'camea2018n000136_2.hdf'))
    failed = []

    for att,val in df.__dict__.items():
        if att in ['name','fileLocation']: # Name and location are to be different
            continue
        if isinstance(val,np.ndarray):
            if val.dtype == 'O':
                continue
            try:
                test = np.all(np.isclose(val,getattr(df2,att)))
            except:
                test = np.all(val==getattr(df2,att))
        else:
            test = val == getattr(df2,att)
        if not test:
            failed.append(att)
    print(failed)
    assert(len(failed)==0)
    os.remove(df2.fileLocation)


def test_DataFile_CreateEmpty(): # TODO: Make this test!!!
    nf = np.array([os.path.join(dataPath,'Normalization_1.calib'),
    os.path.join(dataPath,'Normalization_3.calib'),os.path.join(dataPath,'Normalization_8.calib')])

    A3 = np.linspace(0,180,181)
    A4 = -16
    Ei = 5.5
    Monitor = 1e5
    sample = MJOLNIR.Data.Sample.Sample(a=6.0,b=6.0,c=12.2,projectionVector2=[1,0,0],projectionVector1=[0,2,1],gamma=120.,beta=80.,alpha=90.)

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


def test_updateCalibration():
    calibFiles = [os.path.join(dataPath,'Normalization80_1.calib'),
                    os.path.join(dataPath,'Normalization80_3.calib'),
                    os.path.join(dataPath,'Normalization80_5.calib')]


    df = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    print(df.I)
    print('----------------------')
    df.loadBinning(1)

    binnings = df.possibleBinnings # is 1,3,8
    edges = df.instrumentCalibrationEdges

    df.updateCalibration(calibFiles)

    df.loadBinning(1)
    newBinnings = df.possibleBinnings # is 1,3,8,5
    newEdges = df.instrumentCalibrationEdges
    assert(len(newBinnings)!=len(binnings)) # Addition of binning 5
    assert(not np.any(newEdges!=edges)) # Check if all elemenst are equal


    df.updateCalibration(calibFiles,overwrite=True)
    df.loadBinning(1)

    newEdges = df.instrumentCalibrationEdges
    assert(np.any(newEdges!=edges)) # Check if all elemenst are equal

#
#def test_DataFile_BoundaryCalculation(quick):
#    if quick==True:
#        binning = [1,3,8]
#    else:
#        binning = [1]
#    for B in binning:
#        print('Using binning {}'.format(B))
#        df = DataFile(os.path.join(dataPath,'camea2018n000017.hdf'))
#        converted = df.convert(binning=B)
#        EP,EBins = converted.calculateEdgePolygons()
#        areas = np.array([e.area for e in EP])
#        assert(np.all(areas>2.0)) # Value found by running algorithm
#        assert(len(EBins)==B*8+1)
        

