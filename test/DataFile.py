import numpy as np
from MJOLNIR.Data.DataFile import DataFile,decodeStr,createEmptyDataFile, possibleAttributes, shallowRead
import MJOLNIR.Data.Sample
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pytest

dataPath = 'samlpedata'


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

    files = [os.path.join(dataPath,f) for f in ['camea2018n000137.hdf','camea2018n000137.nxs','camea2022n000894.hdf','camea2022n000894.nxs']]
    DF1 = DataFile(files[0])

    size = DF1.size
    shape = DF1.shape
    assert(size==DF1.I.size)
    assert(np.all(shape == DF1.I.shape))

    try:
        DF1.size = 0
        assert False
    except AttributeError:
        assert True
    try:
        DF1.shape = 0
        assert False
    except AttributeError:
        assert True


def test_DataFile_equility():
    f1 = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    print('----------')
    f2 = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    assert(f1==f2)
    
def test_DataFile_plotA4():
    plt.ioff()
    import matplotlib
    #matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2022n000894.hdf')
    fileName2= os.path.join(dataPath,'camea2022n000894.nxs')
    file = DataFile(fileName)
    

    try:
        file.plotA4(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotA4(1)
    fig2 = file.plotA4()
    
    plt.close('all')

    
def test_DataFile_plotEf():
    plt.ioff()
    import matplotlib
    #matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2022n000894.hdf')
    
    file = DataFile(fileName)

    try:
        file.plotEf(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotEf(1)
    fig2 = file.plotEf()
    
    
    plt.close('all')

def test_DataFile_plotEfOverview():
    plt.ioff()
    import matplotlib
    #matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2022n000894.hdf')

    file = DataFile(fileName)

    try:
        file.plotEfOverview(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotEfOverview(1)
    fig2 = file.plotEfOverview()

    plt.close('all')

def test_DataFile_plotNormalization():
    plt.ioff()
    import matplotlib
    #matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2022n000894.hdf')
    
    file = DataFile(fileName)
    
    try:
        file.plotNormalization(binning=20) # Binning not found in data file
        assert False
    except AttributeError:
        assert True

    fig = file.plotNormalization(1)
    fig2 = file.plotNormalization()

    plt.close('all')

def test_DataFile_decodeString():
    a = b'String'
    b = 'String'

    c =1.1 # Float

    assert(decodeStr(a)==decodeStr(b))
    assert(c == decodeStr(c))

def test_DataFile_ScanParameter():

    files = [os.path.join(dataPath,'camea2022n000894.hdf')]
    for file in files:
        dfile = DataFile(file)
        assert(dfile.scanParameters[0]=='A3')
        assert(len(dfile.scanParameters)==len(dfile.scanUnits))
        assert(len(dfile.scanParameters)==len(dfile.scanValues))
        assert(len(dfile.scanParameters)==1)
        assert(dfile.scanUnits[0]=='degree')
        ##assert(np.all(dfile.scanValues==np.arange(0,150,1)))


def test_DataFile_Error():
    df = DataFile(os.path.join(dataPath,'camea2022n000894.hdf'))

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
    

@pytest.mark.skip(reason="Save/load feature currently under investigation")
def test_DataFile_SaveLoad():
    df = DataFile(os.path.join(dataPath,'camea2022n000894.hdf'))

    df2= DataFile(os.path.join(dataPath,'camea2022n000894_2.hdf'))
    failed = []

    for att,val in df.__dict__.items():
        if att in ['name','fileLocation','fromNICOS']: # Name and location are to be different
            continue
        if att == 'monitor_1':
            print("monitor_1:", repr(val))
            print("type:", type(val), type(val))
            print("shape:", getattr(val, 'shape', None),
                getattr(val, 'shape', None))
            val2 = getattr(df2,att)
            print("monitor_1:", repr(val2))
            print("type:", type(val2), type(val2))
            print("shape:", getattr(val2, 'shape', None),
                getattr(val2, 'shape', None))
        if isinstance(val,np.ndarray):
            if val.dtype == 'O':
                continue
            try:
                test = np.all(np.isclose(val,getattr(df2,att)))
            except:
                test = np.all(val==getattr(df2,att))
        else:
            val2 = getattr(df2,att)
            test = val == val2
        if not test:
            failed.append(att)
    print(failed)

    assert(len(failed)==0)
    os.remove(df2.fileLocation)


def test_DataFile_CreateEmpty(): # TODO: Make this test!!!
    nf = np.array([os.path.join('Data','Normalization_1.calib'),
    os.path.join('Data','Normalization_3.calib'),os.path.join('Data','Normalization_8.calib')])

    A3 = np.linspace(0,180,181)
    A3Position = 30.0
    A4 = -16
    Ei = 5.5
    Monitor = 1e5
    projectionVector1 = np.array([0,2,1,A3Position+30,A4,0.0,0.0,Ei,Ei])
    projectionVector2 = np.array([1,0,0,A3Position,A4,0.0,0.0,Ei,Ei])
    sample = MJOLNIR.Data.Sample.Sample(a=6.0,b=6.0,c=12.2,projectionVector2=projectionVector2,projectionVector1=projectionVector1,gamma=120.,beta=80.,alpha=90.)

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
    calibFiles = [os.path.join('Data','Normalization80_1.calib'),
                    os.path.join('Data','Normalization80_3.calib'),
                    os.path.join('Data','Normalization80_5.calib')]


    df = DataFile(os.path.join(dataPath,'camea2022n000894.hdf'))
    print(df.I)
    print('----------------------')
    df.loadBinning(1)

    binnings = df.possibleBinnings # is 1,3,8
    edges = df.instrumentCalibrationEdges

    df.updateCalibration(calibFiles)

    df.loadBinning(1)
    newBinnings = df.possibleBinnings # is 1,3,8,5
    newEdges = df.instrumentCalibrationEdges
    
    assert(not np.any(newEdges!=edges)) # Check if all elements are equal


    df.updateCalibration(calibFiles,overwrite=True)
    df.loadBinning(1)

    newEdges = df.instrumentCalibrationEdges
    
    assert(np.any(newEdges!=edges)) # Check if all elements are equal

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
        

def test_shallowRead():
    # read out all possible things
    parameters = possibleAttributes
    files = [os.path.join(dataPath,'camea2018n000136.hdf'),
             os.path.join(dataPath,'camea2018n000137.hdf'),
             os.path.join(dataPath,'camea2022n000894.hdf')]
    result = shallowRead(files,parameters)

    assert(len(result) == len(files))
    assert(len(result[0]) == len(parameters))
    assert(len(result[0]) == len(result[1]))

    try:
        shallowRead(files,['NotExisting!'])
        assert False
    except AttributeError:
        assert True