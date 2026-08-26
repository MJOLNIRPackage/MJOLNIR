import numpy as np
from MJOLNIR.Data.DataFile import DataFile,decodeStr,createEmptyDataFile, possibleAttributes, shallowRead
import MJOLNIR.Data.Sample
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pytest

dataPath = 'samlpedata'

@pytest.mark.data
@pytest.mark.unit
def test_DataFile():
    with pytest.raises(FileNotFoundError):
        DF = DataFile('/nope.txt')

    with pytest.raises(FileNotFoundError):
        DF= DataFile(os.path.join(dataPath,'CAMEA_Full.xml')) # Wrong file


    files = [os.path.join(dataPath,f) for f in ['camea2018n000137.hdf','camea2018n000137.nxs','camea2022n000894.hdf','camea2022n000894.nxs']]
    DF1 = DataFile(files[0])

    size = DF1.size
    shape = DF1.shape
    assert(size==DF1.I.size)
    assert(np.all(shape == DF1.I.shape))

    with pytest.raises(AttributeError):
        DF1.size = 0
    with pytest.raises(AttributeError):
        DF1.shape = 0


@pytest.mark.data
@pytest.mark.unit
def test_DataFile_equility():
    f1 = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    print('----------')
    f2 = DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    assert(f1==f2)


@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataFile_plotA4():
    plt.ioff()

    fileName = os.path.join(dataPath,'camea2022n000894.hdf')
    file = DataFile(fileName)
    

    with pytest.raises(AttributeError):
        file.plotA4(binning=20) # Binning not found in data file

    fig = file.plotA4(1)
    fig2 = file.plotA4()
    
    plt.close('all')


@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataFile_plotEf():
    plt.ioff()
    import matplotlib
    #matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2022n000894.hdf')
    
    file = DataFile(fileName)

    with pytest.raises(AttributeError):
        file.plotEf(binning=20) # Binning not found in data file

    fig = file.plotEf(1)
    fig2 = file.plotEf()
    
    
    plt.close('all')

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataFile_plotEfOverview():
    plt.ioff()

    fileName = os.path.join(dataPath,'camea2022n000894.hdf')

    file = DataFile(fileName)

    with pytest.raises(AttributeError):
        file.plotEfOverview(binning=20) # Binning not found in data file

    fig = file.plotEfOverview(1)
    fig2 = file.plotEfOverview()

    plt.close('all')

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataFile_plotNormalization():
    plt.ioff()
    import matplotlib
    #matplotlib.use('Agg')
    fileName = os.path.join(dataPath,'camea2022n000894.hdf')
    
    file = DataFile(fileName)
    
    with pytest.raises(AttributeError):
        file.plotNormalization(binning=20) # Binning not found in data file


    fig = file.plotNormalization(1)
    fig2 = file.plotNormalization()

    plt.close('all')

@pytest.mark.unit
def test_DataFile_decodeString():
    a = b'String'
    b = 'String'

    c =1.1 # Float

    assert(decodeStr(a)==decodeStr(b))
    assert(c == decodeStr(c))


@pytest.mark.data
@pytest.mark.unit
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

@pytest.mark.data
@pytest.mark.unit
def test_DataFile_Error():
    df = DataFile(os.path.join(dataPath,'camea2022n000894.hdf'))

    # Not implimented
    with pytest.raises(NotImplementedError):
        df+df

    df.instrument = 'WrongInstrument'
    with pytest.raises(AttributeError):
        df.convert(binning=1)
    

@pytest.mark.skip(reason="Save/load feature currently under investigation")
@pytest.mark.data
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


@pytest.mark.integration
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

    with pytest.raises(AttributeError):
        _ = createEmptyDataFile(A3=10,A4=10,Ei=10,sample=sample) # No change in any parameter
    
    with pytest.raises(AttributeError):
        _ = createEmptyDataFile(A3=[10,11],A4=[10,11,12],Ei=10,sample=sample) # Two parameters change but not with the same shape
    

    df = createEmptyDataFile(A3=A3,A4=A4,Ei=Ei,sample=sample,Monitor=Monitor,normalizationFiles = nf)
    
    # Check the contents of df
    assert(df.sample == sample)
    assert(len(df.possibleBinnings)==len(nf))

@pytest.mark.data
@pytest.mark.integration
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

       
@pytest.mark.data
@pytest.mark.unit
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

    with pytest.raises(AttributeError):
        shallowRead(files,['NotExisting!'])