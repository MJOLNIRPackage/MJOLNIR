from importlib.util import decode_source
import numpy as np
import MJOLNIR.Data.DataFile
from MJOLNIR.Data.DataSet import DataSet,calculateGrid3D,binData3D,cut1DE,fmt,figureRowColumns,centeroidnp,compareNones,OxfordList, load
from MJOLNIR import _tools
import MJOLNIR.Data.Sample
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os,sys
import warnings
import pytest
from MJOLNIR.Data import Mask

pythonVersion = sys.version_info[0]

dataPath = 'samlpedata'

def test_DataSet_Creation():

    try:
        dataset = DataSet(OtherSetting=10.0,YetAnotherWrongSetting=20.0)
        assert False
    except AttributeError:
        assert True
    dataset = DataSet(Author='Jakob Lass')
    if(dataset.settings['Author']!='Jakob Lass'):
        assert False


def test_Dataset_Initialization():

    emptyDataset = DataSet()
    del emptyDataset
    MJOLNIR.Data.DataFile.assertFile(os.path.join(dataPath,'camea2018n000136.nxs'))
    dataset = DataSet(dataFiles=[os.path.join(dataPath,'camea2018n000136.hdf')],convertedFiles=os.path.join(dataPath,'camea2018n000137.nxs'),calibrationfiles=[])
    
    assert(dataset.dataFiles[0].name=='camea2018n000136.hdf')
    assert(dataset.convertedFiles[0].name=='camea2018n000137.nxs')
    assert(dataset.normalizationfiles == [])
    Str = str(dataset)

                                                                                                                 
def test_DataSet_Error():
    
    try:
        ds = DataSet(normalizationfiles=[10,11])
        assert False
    except:
        assert True

    ds = DataSet()
    
    try: # No data files
        ds.convertDataFile()
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.binData3D(0.1,0.1,0.1)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.cut1D([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.plotCut1D([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.cutQE([0,0],[1,1],0.1,0.01,5.5,6.0)
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.cutPowder(np.linspace(0,4,5),np.linspace(0,4,5))
        assert False
    except AttributeError:
        assert True

    try: # No data files
        ds.plotCutPowder(np.linspace(0,4,5),np.linspace(0,4,5))
        assert False
    except AttributeError:
        assert True
           
    

    try: # Wrong data file type
        ds.dataFiles = 100
        assert False
    except AttributeError:
        assert True


    try: # Can't overwrite settings
        ds.settings={}
        assert False
    except NotImplementedError:
        assert True

    try:# Wrong data file type
        ds.convertedFiles = 10
        assert False
    except AttributeError:
        assert True


    ds.dataFiles = os.path.join(dataPath,'camea2018n000136.hdf')

def test_LoadBambusData():
    ds = DataSet(dataFiles=[os.path.join(dataPath,'BambusTest.dat')])


    ## Set up values to check
    Ei = 3.684
    A4 = 38.194
    A3 = 20.0
    
    assert(ds[0].instrument == 'Bambus')
    assert(ds[0].binning == 1)
    assert(ds[0].dasel == (0,0))
    assert(ds.instrumentCalibrationEf.shape == [(100,4)])

    scanParameters = ds[0].scanParameters
    scanValues = ds[0].scanValues

    assert(len(scanParameters) == len(scanValues))
    

    for param,value in zip(['Ei','A4','A3'],[Ei,A4,A3]):
        assert(np.isclose(getattr(ds[0],param)[0],value))


    try: # Only binning 1 can be used
        ds.convertDataFile(binning = 3)
        assert False
    except AttributeError:
        assert True

    
    assert(ds[0].I.shape == (scanValues.shape[1],100,1))

    ds.convertDataFile()


def test_DataSet_Pythonic():
    dataFiles = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(dataFiles=dataFiles)
    assert(len(dataset)==2)
    for df in dataset:
        print(df)
    initShape = dataset.I.shape
    names = [dataset[i].name for i in range(len(dataset))]
    names.reverse()
    for i,df in enumerate(list(reversed(dataset))):
        names[i]==df.name

    dataset.append(dataFiles)
    assert(len(dataset)==4)
    secondShape = dataset.Monitor.shape
    assert(np.all(secondShape!=initShape))
    del dataset[3]
    del dataset[2]
    try:
        dataset[10]
        assert False
    except IndexError:
        assert True
    
    try:
        del dataset[10]
        assert False
    except IndexError:
        assert True

    try:
        dataset.append('NoFile')
    except:
        assert True
    
    dataset.append(MJOLNIR.Data.DataFile.DataFile(dataFiles[0]))
    assert(len(dataset)==3)
    assert(dataset.I.shape!=secondShape)



def test_DataSet_Equality():
    D1 = DataSet(dataFiles=os.path.join(dataPath,'camea2018n000136.hdf'))#,convertedFiles=['TestData/VanNormalization.nxs')])
    assert(D1==D1)

def test_DataSet_SaveLoad():
    
    D1 = DataSet(dataFiles=os.path.join(dataPath,'camea2018n000136.hdf'))#,convertedFiles = 'TestData/VanNormalization.nxs'))

    temp = 'temporary.bin'

    D1.save(temp)
    D2 = load(temp)
    os.remove(temp)
    assert(D1==D2) 

def test_DataSet_str():
    D1 = DataSet(dataFiles=os.path.join(dataPath,'camea2018n000136.hdf'))#,normalizationfiles = 'TestData/VanNormalization.hdf'))
    string = str(D1)
    print(string)


def test_DataSet_Convert_Data():
    dataFiles = os.path.join(dataPath,'camea2018n000136.hdf')
    dataset = DataSet(dataFiles=dataFiles)
    

    try:
        dataset.convertDataFile(dataFiles=dataFiles,binning=100)
        assert False
    except AttributeError: # Cant find normalization table
        assert True

    try:
        dataset.convertDataFile(dataFiles='FileDoesNotExist',binning=1)
        assert False
    except AttributeError: # FileDoesNotExist
        assert True

    try:
        os.remove(os.path.join(dataPath,'camea2018n000136.nxs'))
    except:
        pass
    dataset.convertDataFile(dataFiles=dataFiles,binning=8,saveLocation=os.path.join(dataPath,''),saveFile=True)
    convertedFile = dataset.convertedFiles[0]
    
    otherFile = MJOLNIR.Data.DataFile.DataFile(dataFiles.replace('.hdf','.nxs'))
    assert(otherFile.difference(convertedFile)==['counts'])
    #assert(otherFile==convertedFile)
    #assert(convertedFile==otherFile)
    os.remove(os.path.join(dataPath,'camea2018n000136.nxs'))
    


def test_DataSet_3DMesh():
    
    x = np.linspace(0,1,2)
    y = np.linspace(0,1,5)
    z = np.linspace(1,2,5)

    X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
    XT1,YT1,ZT1 = calculateGrid3D(X,Y,Z)

    assert(XT1.shape==(3,6,6))
    assert(np.all(XT1[:,0,0]==np.array([-0.5,0.5,1.5])))
    assert(np.all(YT1[0,:,0]==np.array([-0.125,0.125,0.375,0.625,0.875,1.125])))
    assert(np.all(YT1[0,:,0]==ZT1[0,0,:]-1.0))



def test_DataSet_BinData():
    I = np.random.randint(0,100,(10,20,30))
    Norm = np.random.rand(10,20,30)
    Posx = np.linspace(0,1,10)
    Posy = np.linspace(0,1,20)
    Posz = np.linspace(1,2,30)
    PosX,PosY,PosZ = np.meshgrid(Posx,Posy,Posz,indexing='ij')



    pos = [PosX.flatten(),PosY.flatten(),PosZ.flatten()]
    Data,bins = binData3D(0.5,0.25,0.25,pos,I,norm=Norm)

    ReBinnedI = Data[0]
    RebinnedNorm = Data[1]
    RebinnedNormCount = Data[2]


    assert(ReBinnedI.shape==(3,5,5))
    #assert(np.all(bins[0].shape=(4,6,6)))
    assert(RebinnedNorm.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.dtype==int)
    assert(RebinnedNorm.dtype==Norm.dtype)
    assert(ReBinnedI.dtype==I.dtype)


def test_DataSet_full_test():
    import MJOLNIR.Data.Viewer3D
    
    import matplotlib.pyplot as plt
    import os
    plt.ioff()
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf')]

    dataset = DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation=os.path.join(dataPath,''),saveFile=True)
    import matplotlib
    matplotlib.use('Agg')
    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    
    viewer = MJOLNIR.Data.Viewer3D.Viewer3D(Intensity,bins)
    viewer = dataset.View3D(0.08,0.08,0.25,CurratAxeBraggList=[[1,0,0]])
    
    if pythonVersion == 3: # Only possible in python 3
        viewer.ax.set_xticks_base(0.5)
        viewer.ax.set_yticks_base(0.5)

    viewer.setProjection(0)
    viewer.setPlane(4)
    del viewer 
    viewer = dataset.View3D(0.08,0.08,0.25,rlu=False,log=True)
    
    viewer.ax.get_figure().savefig('View3D.png')

    os.remove(os.path.join(dataPath,'camea2018n000136.nxs'))
    os.remove('View3D.png')
    del viewer
    plt.close('all')

def test_DataSet_Visualization():
    import warnings
    from MJOLNIR.Data import Viewer3D,DataFile
    DataFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]

    dataset = DataSet(dataFiles=DataFiles)
    dataset.convertDataFile(saveLocation='Data\\',saveFile=True)

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    Data,bins = dataset.binData3D(0.08,0.08,0.25,dataFiles = [MJOLNIR.Data.DataFile.DataFile(os.path.join('Data','camea2018n000136.nxs'))])
    
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    
    viewer = Viewer3D.Viewer3D(Intensity,bins)
    viewer.caxis = (0,100)
    try:
        viewer.caxis = 'Wrong type'
        assert False
    except AttributeError:
        assert True
    
    try:
        viewer.caxis = [0,1,2,3,4] # Too long input
        assert False
    except AttributeError:
        assert True
    
    try:
        viewer.setAxis(20) # Must bee 0,1, or 2
        assert False
    except AttributeError:
        assert True

    plt.plot()
    plt.close('all')
    os.remove(os.path.join('Data','camea2018n000136.nxs'))
    

def test_DataSet_binEdges():
    X = np.random.rand(100)*3 # array between 0 and 3 -ish
    X.sort()
    tolerance = 0.01
    Bins = _tools.binEdges(X,tolerance=tolerance)

    assert(Bins[0]==X[0]-0.1*tolerance)
    assert(np.isclose(Bins[-1],X[-1],atol=5) or Bins[-1]>X[-1])
    assert(len(Bins)<=3.0/tolerance)
    assert(np.all(np.diff(Bins[:-1])>tolerance*0.99))

def test_DataSet_1Dcut():
    q1 =  np.array([1.23,-1.51])
    q2 =  np.array([1.54, -1.25])
    width = 0.1

    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile(saveFile=False)

    ax,Data,bins = ds.plotCut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,fmt='.')
    Data2,bins2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5)
    
    # Check that the two data sets have the same values (except for Data2 also having 'BinDistance')
    #print(Data2)
    #print(Data.loc[:, Data.columns != 'BinDistance'])
    assert(Data2.equals(Data.loc[:, Data.columns != 'BinDistance']))

    Data,bins = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,extend=False)
    assert(np.all(np.logical_and(bins[0][:,0]>=q1[0]-0.1,bins[0][:,0]<=q2[0]+0.1))) 
    # x-values should be between 1.1 and 2.0 correpsonding to q points given (add some extra space due to way bins are created (binEdges))

    #q3 = np.array([1.1,1.1])
    #q4 = np.array([2.0,2.0])
    Data,bins = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,extend=False)
    
    assert(np.all(bins[0][:,0]>=q1[0]-0.1))
    assert(np.all(bins[0][:,0]<=q2[0]+0.1))
    assert(np.all(bins[0][:,1]>=q1[1]-0.1))
    assert(np.all(bins[0][:,1]<=q2[1]+0.1))
    # x and y-values should be between 1.1 and 2.0 correpsonding to q points given (add some extra space due to way bins are created (binEdges))

    Q1 = np.array([1,0,0])
    Q2 = np.array([0.5,1,0])

    ax,Data,bins = ds.plotCut1D(Q1,Q2,width,rlu=True,minPixel=0.01,Emin=2.0,Emax=2.5,fmt='.')
    Data2,bins2 = ds.cut1D(Q1,Q2,width,rlu=True,minPixel=0.01,Emin=2.0,Emax=2.5)

    assert(Data2.equals(Data.loc[:,Data.columns!='BinDistance']))
    assert(np.all(np.array([np.all(np.isclose(bins[i],bins2[i])) for i in range(len(bins))]).flatten()))

    q1,q2 = ds.convertToQxQy([Q1,Q2])
    D1,b1 = ds.cut1D(Q1,Q2,width,rlu=True,minPixel=0.01,Emin=2.0,Emax=2.5)
    D2,b2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5)

    # Convert b1 to HKL in order to check if conversion works
    BinPos,OrthoPos,E = b1
    BinPos = np.concatenate([ds.convertToQxQy(BinPos[:,:3]),BinPos[:,-1].reshape(-1,1)],axis=1)
    OrthoPos = ds.convertToQxQy(OrthoPos)
    b1 = [BinPos,OrthoPos,E]

    
    assert(np.all(np.isclose(D1,D2)))
    
    assert(np.all(np.array([np.all(np.isclose(b1[i],b2[i])) for i in range(len(b1))]).flatten()))

    # Check that generating a plot from previous data is equivalent to directly plotting

    
    ax1,cut,bins = ds.plotCut1D(Q1,Q2,width=0.1,minPixel=0.04,Emin=2.0,Emax=2.5,ufit=False)

    ax2,*_ = ds.plotCut1D(Q1,Q2,Emin=2.0,Emax=2.5,width=0.1,minPixel=0.04,data=[cut,bins])



    for key in ['legend','title','xlabel','ylabel']:
        val1 = ax1.properties()[key]
        val2 = ax2.properties()[key]
        assert(val2==val1)

    # check plotted data in line plot 
    # ax1, find index of line2D
    id1 = np.arange(len(ax1.properties()['children']))[np.array([isinstance(child,mpl.lines.Line2D) for child in ax1.properties()['children']])][0]
    id2 = np.arange(len(ax2.properties()['children']))[np.array([isinstance(child,mpl.lines.Line2D) for child in ax2.properties()['children']])][0]

    np.testing.assert_array_equal(ax2.properties()['children'][id1]._xy,ax1.properties()['children'][id2]._xy)


def test_DataSet_1Dcut_ufit():
    q1 =  np.array([1.23,-1.51])
    q2 =  np.array([1.54, -1.25])
    width = 0.1

    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile(saveFile=False)

    ax,dataset = ds.plotCut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,fmt='.',ufit=True)
    dataset2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,Emin=2.0,Emax=2.5,ufit=True)
    

    files = ', '.join([x.replace('hdf','nxs').split(os.path.sep)[-1] for x in convertFiles])

    assert(np.all([np.all(x==y) for x,y in zip(dataset.fit_columns,dataset2.fit_columns)]))
    assert(dataset.meta == dataset2.meta)

    assert(dataset.meta['instrument'] == 'CAMEA')
    assert(dataset.meta['datafilename'] == files)

    

def test_DataSet_1DcutE():
    q =  np.array([1.23,-1.25]).reshape(2,1)
    width = 0.1
    Emin = 1.5
    Emax = 2.5
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000137.hdf')]
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile(saveFile=True)
    Datset._getData()
    I,qx,qy,energy,Norm,Monitor = Datset.I.extractData(),Datset.qx.extractData(),Datset.qy.extractData(),Datset.energy.extractData(),Datset.Norm.extractData(),Datset.Monitor.extractData()

    [intensity,MonitorCount,Normalization,normcounts],[bins] = cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=Emin,E2=Emax,q=q,width=width,minPixel=0.01)
    Q = Datset.convertToHKL(q.reshape(2))
    
    Data,[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=Q,width=width,minPixel=0.01)
    assert(np.min(bins)>=Emin-0.01) # Check that bins do not include data outside of cut
    assert(np.max(bins)<=Emax+0.01)
    assert(len(bins)==len(intensity)+1)# Bins denotes edges and must then be 1 more than intensity

    assert(intensity.shape==MonitorCount.shape) # Check that all matrices are cut equally
    assert(intensity.shape==Normalization.shape)
    assert(intensity.shape==normcounts.shape)

    Data,[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=q,width=width,minPixel=0.01,rlu=False)
    
    Data,[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    
    assert(np.all(np.isclose(np.diff(bins),0.01)))
    assert(bins.min()>=Emin)
    assert(bins.max()<=Emax)

    try: # no points inside energy interval
        cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=500,E2=700,q=q,width=width,minPixel=0.01)
        assert False
    except AttributeError:
        assert True

    try: # no points inside q
        cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=5,E2=7,q=np.array([20.0,0]).reshape(2,1),width=width,minPixel=0.01)
        assert False
    except AttributeError:
        assert True

    # Check the data of plot to be the same as cut
    Data,[bins] = Datset.cut1DE(E1=Emin,E2=Emax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    ax,Data2,[bins2] = Datset.plotCut1DE(E1=Emin,E2=Emax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    assert(Data.equals(Data2))
    assert(np.all(np.isclose(bins,bins2)))

    ufitData = Datset.cut1DE(E1=Emin,E2=Emax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True,ufit=True)
    ax,ufitData2 = Datset.plotCut1DE(E1=Emin,E2=Emax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True,ufit=True)

    files = ', '.join([x.replace('hdf','nxs').split(os.path.sep)[-1] for x in convertFiles])
    
    assert(np.all([np.all(np.isclose(x,y,equal_nan=True)) for x,y in zip(ufitData.fit_columns,ufitData2.fit_columns)]))
    assert(ufitData.meta == ufitData2.meta)

    assert(ufitData.meta['instrument'] == 'CAMEA')
    assert(ufitData.meta['datafilename'] == files)

    ax,Data3,[bins3] = Datset.plotCut1DE(E1=Emin,E2=Emax,q=Q,width=0.1,minPixel=0.01,constantBins=True)



def test_DataSet_2Dcut():
    q1 =  np.array([1.23,-1.25])
    q2 =  np.array([1.54, -1.51])
    width = 0.1
    minPixel=0.02
    EnergyBins = np.linspace(2,3,4)
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000137.hdf')]

    Datset = DataSet(dataFiles = convertFiles)

    Datset.convertDataFile(saveFile=False)
    ax,Data,bins = Datset.plotCutQE(q1,q2,width=width,minPixel=minPixel,EnergyBins=EnergyBins,rlu=False)# Remove to improve test coverage ,vmin=0.0 , vmax= 5e-06)
    Data2,bins2,_ = Datset.cutQE(q1,q2,width=width,minPixel=minPixel,EnergyBins=EnergyBins,rlu=False)

    comparisons = list(Data.columns)
    del comparisons[comparisons.index('BinDistance')]

    np.all([np.allclose(Data[p],Data2[p],equal_nan=True) for p in comparisons])

    Q1 = np.array([1,0,0])
    Q2 = np.array([0.5,1,0])

    q1,q2 = Datset.convertToQxQy([Q1,Q2])

    Data1,bins,_ = Datset.cutQE(Q1,Q2,width,minPixel,EnergyBins=EnergyBins,rlu=True)
    Data2,bins2,_ = Datset.cutQE(q1,q2,width,minPixel,EnergyBins=EnergyBins,rlu=False)

    comparisons = list(Data.columns)
    del comparisons[comparisons.index('BinDistance')]
    
    assert(np.all(comparisons == Data2.columns))

    
    np.all([np.allclose(Data1[p],Data2[p],equal_nan=True) for p in comparisons])

def test_DataSet_cutPowder():

    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile(saveFile=True)
    mask = Datset.mask#np.ones_like(Datset.mask)[0]

    Datset.mask = mask
    Datset.mask = np.logical_not(mask)
    

    eBins = _tools.binEdges(Datset.energy,0.25)
    qBins = np.linspace(0.0,4,201)

    ax,D,QCentres,ECentres = Datset.plotCutPowder(eBins,qBins)# Remove to improve test ,vmin=0,vmax=1e-6)
    D2 = Datset.cutPowder(eBins,qBins)
    assert(np.all(D.equals(D2)))
    #assert(np.all(np.isclose(QCentres2,QCentres)))
    #assert(np.all(np.isclose(ECentres2,ECentres)))

    
def test_DataSet_createRLUAxes():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    fig = plt.figure()
    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile(saveFile=True)

    ax = ds.createQAxis()
    ax = ds.createQAxis(basex=0.5,figure=fig)
    ax = ds.createQAxis(basey=0.5)

    if pythonVersion == 3: # Only possible in python 3
        ax.set_xticks_base(0.2)
        ax.set_yticks_base(0.5)

    V1,V2,V3 = [2,0,0],[-2,3,0],[2,-3,0]
    ax.set_axis(V1,V2)
    ax.set_axis(V1,V2,V3)

    plt.close('all')


def test_DataSet_createQEAxes():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile(saveFile=True)

    ax = ds.createQEAxes(projectionVector1=ds.sample[0].projectionVector1,projectionVector2=ds.sample[0].projectionVector2)

    try:
        ax = ds.createQEAxes(axis=2) # Axis only allowed to be 0 or 1
    except AttributeError:
        assert True
    
    try:
        ax = ds.createQEAxes(projectionVector1=[1,0,0],projectionVector2=[1,2,3,4,5]) # Wrong shape of vector
    except AttributeError:
        assert True
    plt.close('all')



def test_DataSet_plotQPlane():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000137.hdf')]#'TestData/ManuallyChangedData/A3.hdf')]
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile(saveFile=True)

    EmptyDS = DataSet()
    try:
        Datset.plotQPlane() # No Bins, Emin or Emax
        assert False
    except AttributeError:
        assert True
    try:
        Datset.plotQPlane(EBins=[10]) # Length of bins is 1
        assert False
    except AttributeError:
        assert True
    
    try:
        Datset.plotQPlane(EMin=20,EMax=10) # EMin>EMax
        assert False
    except AttributeError:
        assert True
    
    try:
        EmptyDS.plotQPlane(EMin=2,EMax=3) # Empty DataSet
        assert False
    except AttributeError:
        assert True


    EMin = Datset.energy.min()
    EMax = EMin+0.5
    Data,ax1 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=True,log=False,rlu=True)
    Data,ax2 = Datset.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=True,rlu=True)
    AX = Datset.createQAxis()
    Data,ax3 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,ax=AX,colorbar=True,vmin=0,vmax=1e-6,zorder=10)
    
    ax1.set_clim(-20,-15)
    ax2.set_clim(0,1e-6)
    Data,ax3 = Datset.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.05,yBinTolerance=0.05)
    
    cmap = plt.cm.coolwarm

    Dataset = DataSet(dataFiles=convertFiles)
    for d in Dataset.dataFiles:
        d.A3Off +=90 # rotate data to fall into problem of arctan2
    Data,ax2 = Datset.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.05,yBinTolerance=0.05,enlargen=False,log=True,rlu=True,cmap=cmap)
    
    try:
        Datset.plotQPlane(EMin,EMax,binning='notABinningMethod')
        assert False
    except:
        assert True

    # 3D
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.colors import ListedColormap
    cmap = plt.cm.coolwarm
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)

    fig = plt.figure(figsize=(10,11))
    ax = fig.add_subplot(111, projection='3d')

    Energies = np.concatenate(Datset.energy,axis=0)
    E = np.arange(Energies.min()+0.35,Energies.max(),0.35)


    PandasData,ax = \
    Datset.plotQPlane(EBins=E,ax = ax,xBinTolerance=0.03,yBinTolerance=0.03,
            binning='polar',vmin=7.5e-7,vmax=7e-6,antialiased=True,cmap=cmap,rlu=True,extend='max')
    plt.close('all')

#@pytest.mark.unit
#def test_DataSet_plotA3A4(quick):
#    plt.ioff()
#
#    File1 = os.path.join(dataPath,'camea2018n000136.hdf')
#    File2 = os.path.join(dataPath,'camea2018n000137.hdf')
#
#    DS = DataSet(dataFiles=[File1,File2])
#    DS.convertDataFile(saveFile=True)
#
#    F1 = DS.convertedFiles[0]
#    F2 = DS.convertedFiles[1]
#
#    files = [F1,F2]
#    axes = [plt.figure().gca(),plt.figure().gca()]
#    try:
#        plotA3A4(files,planes=[],ax=axes) # 64 planes and only 2 axes
#        assert False
#    except AttributeError:
#        assert True
#
#    try:
#        plotA3A4(files,planes=None,ax=[]) # 64 planes and only 2 axes
#        assert False
#    except AttributeError:
#        assert True 
#
#    try:
#        plotA3A4(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
#        assert False
#    except AttributeError:
#        assert True
#    
#    try:
#        ei = F1.Ei
#        F2.Ei = F1.Ei*10
#        plotA3A4(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
#        assert False
#    except AttributeError:
#        F2.Ei = ei
#        assert True
#
#    try:
#        plotA3A4(files,planes=[10,[22]],ax=axes,singleFigure=True) # 2 axes and singleFigure true
#        assert False
#    except AttributeError:
#        assert True
#    if not quick==True:
#        print('___________')
#        plotA3A4(files,planes=[10,[22,23]],ax=axes) # Plot plane 10 and 22+23 in the provided axes
#        print('___________')
#        DS.plotA3A4(planes=[19,[22,25]]) # Plot planes in new axes
#        print('___________')
#        DS.plotA3A4([F1,F1],planes=[19,[22,25]]) # Plot planes in new axes
#        print('___________')
#        patches,energies=DS.plotA3A4([F1],planes=[10,25],returnPatches=True)
#        print('___________')
#        assert(len(patches)==2)
#        assert(len(energies)==2)
#    plt.close('all')

@pytest.mark.unit
def test_DataSet_plotQPatches(quick):
    assert True
    #     plt.ioff()

#     File1 = 'TestData/T0Phonon10meV.nxs')
#     File2 = 'TestData/T0Phonon10meV93_5A4.nxs')

#     DS = DataSet(convertedFiles=[File1,File2])

#     F1 = DS.convertedFiles[0]
#     F2 = DS.convertedFiles[1]

#     files = [F1,F2]
#     axes = [plt.figure().gca(),plt.figure().gca()]
#     try:
#         plotQPatches(files,planes=[],ax=axes) # 64 planes and only 2 axes
#         assert False
#     except AttributeError:
#         assert True
        
#     try:
#         plotQPatches(files,planes=[[0,2,3],23,44],ax=axes) # 3 planes and 2 axes
#         assert False
#     except AttributeError:
#         assert True

#     try:
#         plotQPatches(files,planes=[10,[22]],ax=axes,singleFigure=True) # 2 axes and singleFigure true
#         assert False
#     except AttributeError:
#         assert True

#     if not quick==True:
#         plotQPatches(files,planes=[10,[22,23]],ax=axes) # Plot plane 10 and 22+23 in the provided axes
#         DS.plotQPatches(planes=[19,[22,25]],A4Extend=0.5,A3Extend=1) # Plot planes in new axes
#         DS.plotQPatches(dataFiles=[files[0],files[0]],planes=[19,[22,25]],A4Extend=0.5,A3Extend=1) # Plot planes in new axes and only one file
#     plt.close('all')
    

def test_DataSet_fmt():
    assert('$1.00 \\times 10^{1}$' == fmt(10,'Unused'))
    assert('$1.00 \\times 10^{-10}$' == fmt(1e-10,'Unused'))
    assert('$2.55 \\times 10^{-2}$' == fmt(0.0255,'Unused'))
    assert('$2.56 \\times 10^{-2}$' == fmt(0.02556,'Unused'))
    

def test_DataSet_figureRowColumns():
    assert(np.all(np.array([3,4])==np.array(figureRowColumns(10)))) # 10 -> 3,4
    assert(np.all(np.array([3,3])==np.array(figureRowColumns(9)))) # 9 -> 3,3
    assert(np.all(np.array([1,1])==np.array(figureRowColumns(1)))) # 1 -> 1,1
    try:
        figureRowColumns(0) # No figures
        assert False
    except AttributeError:
        assert True

    assert(np.all(np.array([8,8])==np.array(figureRowColumns(63)))) # 63 -> 8,8
    

def test_DataSet_centeroidnp():
    pos = np.array([[0,0],[1,0],[0,1],[1,1]],dtype=float)
    assert(np.all(np.isclose(np.array([0.5,0.5]),centeroidnp(pos))))

    pos2 = np.array([[1.2,2.2],[7.5,1.0],[11.0,0.0],[4.0,-1.0],[2.0,2.0]],dtype=float)
    assert(np.all(np.isclose(np.array([5.14,0.84]),centeroidnp(pos2))))
    
def test_DataSet_compareNones():
    assert(compareNones(np.array([None]),np.array([None]),0.1))
    assert(not compareNones(np.array([None]),np.array([0.5]),0.1))
    assert(not compareNones(np.array([0.5]),np.array([None]),0.1))
    assert(compareNones(np.array([0.4]),np.array([0.5]),0.2))
    assert(not compareNones(np.array([0.4]),np.array([0.5]),0.001))

    assert(not np.all(compareNones(np.array([0.4,10.2,10.0]),np.array([0.5]),0.001)))
    assert(np.all(compareNones(np.array([0.4,10.2,10.0]),np.array([0.4,10.2,10.0]),0.001)))


def test_DataSet_cutQELine():
    QPoints = np.array([[0.3,-1],[0.7,-1.4],[1.6,-0.9],[0.3,-0.9]],dtype=float)
    QPointsHKL=np.array([[1.0,0.0,0.0],
                        [0.5,1.5,0.0],
                        [1.7,-0.1,0.0],
                        [1.0,1.0,0.0]])


    EnergyBins = np.linspace(1.7,2.7,5)
    minPixel = 0.001
    width=0.1
    DataFile = [os.path.join(dataPath,'camea2018n000137.hdf')]

    dataset = DataSet(convertedFiles=DataFile)
    dataset.convertDataFile(saveFile=False)

    try: # No Q-points
        dataset.cutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu=True)
        assert False
    except AttributeError:
        assert True
    try: # Wrong RLU-input
        dataset.cutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu='42') # Wrong RLU-input
        assert False
    except AttributeError:
        assert True


    DataList,BinList=dataset.cutQELine(QPointsHKL,EnergyBins,width=width,minPixel=minPixel,rlu=True)
    DataList2,BinList2=dataset.cutQELine(QPoints,EnergyBins,width=width,minPixel=minPixel,rlu=False)

    assert(len(DataList['qCut'][0])==(len(QPointsHKL)-1)*(len(EnergyBins)-1))# Assert that there are 3 cuts with 4 energies

    assert(len(DataList2['qCut'][0])==(len(QPoints)-1)*(len(EnergyBins)-1))

    
    

def test_DataSet_plotCutQELine():
    
    Points = np.array([[0.7140393034102988,-0.4959224853328328],
                        [1.128363301356428,-1.6520150761601147],
                        [1.9002545852012716,-0.9393552598967219],
                        [1.0432282332853056,-0.12375569239528339]],dtype=float)
    QPoints = np.zeros((Points.shape[0],3))
    QPoints[:,:2]=Points
    EnergyBins = np.linspace(1.7,2.7,11)
    minPixel = 0.001
    width=0.1
    import matplotlib
    matplotlib.use('Agg')

    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(convertedFiles=DataFile)
    dataset.convertDataFile(saveFile=False)
    
    try: # No Q-points
        dataset.plotCutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu=False)
        assert False
    except AttributeError:
        assert True

    try: # No wrong dim of QPonts
        dataset.plotCutQELine(QPoints,EnergyBins,width=width,minPixel=minPixel,rlu=False)
        assert False
    except AttributeError:
        assert True

    try: # No wrong dim of QPonts
        dataset.plotCutQELine(QPoints[:,:2],EnergyBins,width=width,minPixel=minPixel,rlu=True)
        assert False
    except AttributeError:
        assert True


    ax,Data,Bins = dataset.plotCutQELine(
        QPoints[:,:2],EnergyBins,width=width,minPixel=minPixel,rlu=False,vmin=0.0,vmax=1.5e-6,log=True,seperatorWidth=3)


    HKLPoints = np.array([[1.0,0.0,0.0],
                        [0.5,1.5,0.0],
                        [1.7,-0.1,0.0],
                        [1.0,1.0,0.0]])



    ax,Data,Bins = dataset.plotCutQELine(
        HKLPoints,EnergyBins,width=width,minPixel=minPixel,rlu=True,plotSeperator = False,colorbar=True,log=True)

    
    if True:
        # 3D
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib.colors import ListedColormap
        cmap = plt.cm.coolwarm
        my_cmap = cmap(np.arange(cmap.N))
        my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
        my_cmap = ListedColormap(my_cmap)

        fig = plt.figure(figsize=(10,11))
        ax = fig.add_subplot(111, projection='3d')

        Energies = np.concatenate(dataset.energy,axis=0)
        E = np.arange(Energies.min()+0.35,Energies.max(),0.35)
        
        try:
            ax,Data,Bins = \
            dataset.plotCutQELine(QPoints=HKLPoints,EnergyBins=E,ax = ax,width=0.05,minPixel=0.01,
                    vmin=7.5e-7,vmax=7e-6,cmap=cmap,rlu=True)
            assert False
        except NotImplementedError:
            assert True



def test_DataSet_extractDetectorData():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]#['TestData/ManuallyChangedData/A3.nxs'),'TestData/ManuallyChangedData/A3.nxs')]
    dataset = DataSet(DataFile)

    binning = 1
    dataset.convertDataFile(binning=binning,saveFile=True)

    
    try:
        dataset.extractDetectorData(A4=10000.0) # A4 outside of detector
        assert False
    except AttributeError:
        assert True

    try:
        dataset.extractDetectorData(Ef=10000.0) # Ef outside of detector
        assert False
    except AttributeError:
        assert True
    

    Efs = dataset.convertedFiles[0].instrumentCalibrationEf[:,1].reshape(104,8*binning)
    AnalyserSelection = 5
    Ef = np.mean(Efs[:,AnalyserSelection])

    A4s = dataset.convertedFiles[0].instrumentCalibrationA4.reshape(104,8*binning)
    DetectorSelection = 19
    A4 = np.mean(A4s[DetectorSelection])-dataset.convertedFiles[0].A4Off


    DatBoth = dataset.extractData(A4=A4,Ef=Ef)
    DatBothId = dataset.extractData(A4=A4,EfId=AnalyserSelection)
    DatOne = dataset.extractData(A4=A4)
    DatOne2= dataset.extractData(Ef=Ef)
    DatAll = dataset.extractData()
    DatAllRaw = dataset.extractData(raw=True)


    # Independent of number of files:
    assert(len(DatAllRaw)==3) # Check that 3 lists are returned
    assert(len(DatAllRaw[0])==len(DatAllRaw[1]) and len(DatAllRaw[0])==len(DatAllRaw[2])) # Check that 3 list have same number of files

    assert(np.all(DatBothId[0]==DatBoth[0])) # Check that ID and value gives the same.

    # The shape of raw is the same as non-raw
    assert(len(DatAllRaw[0])==len(DatAll)) # Have same number of files

    for i in range(len(dataset.convertedFiles)):
        assert(DatAllRaw[0][i].shape==DatAllRaw[1][i].shape and DatAllRaw[0][i].shape==DatAllRaw[2][i].shape) # Check that 3 list have same shape
        assert(DatAllRaw[0][i].shape==DatAll[i].shape) 
        

def test_DataSet_subract():
    #Simple test of subtracting the same data file from it-self
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)

    dataset2 = DataSet(DataFile)

    try:
        subtracted = dataset-dataset2
    except AttributeError: # Data set need to be converted
        assert True
    
    dataset.convertDataFile()
    dataset2.convertDataFile()
    subtracted = dataset-dataset2
    assert(np.sum(subtracted.I.data)==0)
    assert(np.all(np.isclose(subtracted.Norm.data,dataset.Norm.data)))
    assert(np.all(np.isclose(subtracted.Monitor.data,dataset.Monitor.data)))

def test_DataSet_OxfordList():
    l = ['Apples','Pears']
    S = OxfordList(l)
    assert(S=='Apples and Pears')

    l.append('Oranges')
    S = OxfordList(l)
    assert(S=='Apples, Pears, and Oranges')

    assert(OxfordList([]) is None)
    assert(OxfordList(['Apples'])=='Apples')


def test_DataSet_MultiFLEXX():
    fileLocation = _tools.fileListGenerator('65059',folder=os.path.join('Data',''),instrument='MultiFLEXX')

    ds = DataSet(fileLocation)
    ds.convertDataFile(saveFile = False)
    import matplotlib
    matplotlib.use('Agg')

    V = ds.View3D(0.05,0.05,0.5,grid=True)

def test_DataSet_ELine():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)

    dataset.convertDataFile()
    import matplotlib
    matplotlib.use('Agg')

    Q1 = [1.0,-0.185,0.0]
    Q2 = [0.5,1.5,0.0]
    Emin = 1.65
    Emax = 3.3

    CutData,Bins = dataset.cutELine(Q1, Q2, Emin=Emin, Emax=Emax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=True, dataFiles=None, constantBins=False)
    ax, CutDataPlot, BinsPlot = dataset.plotCutELine(Q1, Q2, Emin=Emin, Emax=Emax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=True, dataFiles=None, constantBins=False)

    assert(np.all([np.all(np.isclose(B,B2)) for B,B2 in zip(Bins,BinsPlot)]))
    assert(np.all(CutDataPlot.equals(CutData)))

    assert(np.all(np.logical_and(CutData['Energy']>=Emin,CutData['Energy']<=Emax*1.01))) # Allow for slightly heigher energy
    assert(np.logical_and(np.all(CutData['H']<=Q1[0]*1.01),np.all(CutData['H']>=Q2[0]*0.99)))
    assert(np.logical_and(np.all(CutData['K']>=Q1[1]*1.01),np.all(CutData['K']<=Q2[1]*1.01)))
    assert(np.all(np.isclose(CutData['L'],0.0,atol=1e-6)))
    assert(np.all([np.all(np.logical_and(B[0]>=Emin*0.99,B[0]<=Emax*1.05)) for B in Bins])) # Allow for slightly heigher energy


    if sys.version[0] == '3':
        ax.set_xticks_base(0.01)
        ax.set_xticks_number(10)
        ax.set_xticks_base()
        

    Q1raw = dataset.convertToQxQy(Q1)
    Q2raw = dataset.convertToQxQy(Q2)

    CutData,Bins = dataset.cutELine(Q1raw, Q2raw, Emin=Emin, Emax=Emax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=False)
    ax, CutDataPlot, BinsPlot = dataset.plotCutELine(Q1raw, Q2raw, Emin=Emin, Emax=Emax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=False)


    assert(np.all([np.all(np.isclose(B,B2)) for B,B2 in zip(Bins,BinsPlot)]))
    assert(np.all(CutDataPlot.equals(CutData)))

    assert(np.all(np.logical_and(CutData['Energy']>=Emin,CutData['Energy']<=Emax*1.01))) # Allow for slightly heigher energy
    assert(np.logical_and(np.all(CutData['Qx']<=Q1raw[0]*1.01),np.all(CutData['Qx']>=Q2raw[0]*0.99)))
    assert(np.logical_and(np.all(CutData['Qy']<=Q1raw[1]*0.99),np.all(CutData['Qy']>=Q2raw[1]*1.01)))

    assert(np.all([np.all(np.logical_and(B[0]>=Emin*0.99,B[0]<=Emax*1.05)) for B in Bins])) # Allow for slightly heigher energy


def test_updateCalibration():
    calibFiles = [os.path.join('Data','Normalization80_1.calib'),
                    os.path.join('Data','Normalization80_3.calib'),
                    os.path.join('Data','Normalization80_5.calib')]


    ds = DataSet(os.path.join(dataPath,'camea2018n000136.hdf'))
    
    df = ds[0]

    df.loadBinning(1)

    binnings = df.possibleBinnings # is 1,3,8
    edges = df.instrumentCalibrationEdges

    ds.updateCalibration(calibFiles)

    df.loadBinning(1)
    newBinnings = df.possibleBinnings # is 1,3,8,5
    newEdges = df.instrumentCalibrationEdges
    assert(len(newBinnings)!=len(binnings)) # Addition of binning 5
    assert(not np.any(newEdges!=edges)) # Check if all elemenst are equal


    ds.updateCalibration(calibFiles,overwrite=True)
    df.loadBinning(1)

    newEdges = df.instrumentCalibrationEdges
    assert(np.any(newEdges!=edges)) # Check if all elemenst are equal


def testplotRaw1D_Error():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
     # Scan variables are A3 and A3+A4
    dataset = DataSet(DataFile)
    dataset[0].scanParameters = ['Ei']
    import matplotlib
    matplotlib.use('Agg')

    try:
        dataset.plotRaw1D() # Two different scan types
        assert False
    except AttributeError:
        assert True
        
    DataFile = [os.path.join(dataPath,'camea2018n000137.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)
    try:
        dataset.plotRaw1D(analyzerSelection=[0]) # Not enough analyzers provided
        assert False
    except AttributeError:
        assert True
    try:
        dataset.plotRaw1D(detectorSelection=[0]) # Not enough detectors provided
        assert False
    except AttributeError:
        assert True
        
    try:
        dataset.plotRaw1D(legend=[0]) # Not enough legend labels provided
        assert False
    except AttributeError:
        assert True
        
def testplotRaw1D():
    DataFile = [os.path.join(dataPath,'camea2018n000137.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)
    import matplotlib
    matplotlib.use('Agg')

    ax = dataset.plotRaw1D()
    
    ax = dataset.plotRaw1D(legend=['1','2'],detectorSelection=[0,0],analyzerSelection=[5,5],grid=True)
    assert True

def testMasking():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
     # Scan variables are A3 and A3+A4
    ds = DataSet(DataFile)
    ds.convertDataFile()

    circ = Mask.circleMask(center=np.array([1.0,0.0]),radiusPoint=np.array([1.1,0.0]),coordinates =['h','k'])
    Emask = Mask.lineMask(1.7,end=2.0,coordinates='energy')
    rect = Mask.rectangleMask(corner1=np.array([1.0,0.0]),corner2=np.array([1.5,0.5]),coordinates=['h','k'])

    mask = circ*Emask+Emask*rect
    calcMask = (circ*Emask)(ds)
    ds.mask = mask(ds)
    masks = [mask(df) for df in ds]
    
    assert(np.sum(ds.I.mask)==np.sum([np.sum(m) for m in masks]))


    # Test failing mask
    try:
        ds.mask = 'Error'
        assert False
    except AttributeError:
        assert True
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mask = np.array([np.ones_like(I) for I in ds.I]) # Mask all
        ds.mask = mask

    ds = DataSet([DataFile[0],DataFile[0]])
    ds.convertDataFile()

    
    print([x.shape for x in mask])
    print(mask[0].shape)
    
    ds.mask = mask[0] # Only provide 1 mask to be applied to both data files
    
    

def testupdateSampleParameters():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    ds = DataSet(DataFile)

    unitCell = ds.sample[0].unitCell
    unitCell[:2]*=2.0
    ds.updateSampleParameters(unitCell=unitCell)
    newCell = ds.sample[0].cell
    newUB = ds.sample[0].UB

    for d in ds:
        assert(np.all(np.isclose(d.sample.cell,newCell)))
        assert(np.all(np.isclose(d.sample.UB,newUB)))

def test_CurratAxeMasking():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    ds = DataSet(DataFile)
    
    # Correct wrong alignment of about 3.5 degrees
    theta=np.rad2deg(0.05979407243252655)
    for df in ds:
        df.A3-=theta
    ds.convertDataFile(binning=3)

    averageSum = np.mean(ds.I.extractData())

    braggpeaks = [[1,0,0]]
    # Calcualte mask directly from dataset
    curratAxeMask = ds.calculateCurratAxeMask(braggpeaks,dH=0.1,dK=0.1,dL=0.001)
    ds.mask = curratAxeMask

    # Calcualte mask using a mask object and on the dataset
    CAMask = Mask.CurratAxeMask(braggpeaks,dH=0.1,dK=0.1,dL=0.001)

    evaluatedCAMask = CAMask(ds)
    assert(np.all(np.equal(evaluatedCAMask,curratAxeMask)))

    # Calcualte the two different contributions
    CAMaskMono = Mask.CurratAxeMask(braggpeaks,dH=0.1,dK=0.1,dL=0.001,spurionType='Monochromator')
    CAMaskAna = Mask.CurratAxeMask(braggpeaks,dH=0.1,dK=0.1,dL=0.001,spurionType='Analyser')
    CAMaskMonoEvaluated = CAMaskMono(ds)
    CAMaskAnaEvaluated = CAMaskAna(ds)
    CAMaskCombi = [np.logical_or(CAMM,CAMA) for CAMM,CAMA in zip(CAMaskMonoEvaluated,CAMaskAnaEvaluated)]
    assert(np.all(np.equal(evaluatedCAMask,CAMaskCombi)))


    averageSumMasked = np.mean(ds.I.extractData())
    assert(averageSum-averageSumMasked>0.004)

    curratAxeMask = ds.calculateCurratAxeMask(braggpeaks,dqx=0.15,dqy=0.15)
    ds.mask = curratAxeMask


    averageSumMasked = np.mean(ds.I.extractData())
    assert(averageSum-averageSumMasked>0.004)
    


def test_absoluteNormalziation():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf')]
    ds = DataSet(DataFile)

    try:
        ds.absoluteNormalize(10.0,'MnF2')
        assert False
    except  AttributeError: # Must be converted first!
        assert True

    ds.convertDataFile()
    norm = np.mean(ds.Norm.extractData())
    
    # Use value for MnF2 to check 
    ds.absoluteNormalize(sampleMass=6.2,sampleChemicalFormula='MnF2',formulaUnitsPerUnitCell=2,
                                      correctVanadium=False)
    
    factor = 0.06088201383247563 # Factor calculated for MnF2

    assert(np.isclose(factor,ds.absoluteNormalized))

    # Redo normalization to retrive the same factor
    ds.absoluteNormalize(sampleMass=6.2,sampleChemicalFormula='MnF2',formulaUnitsPerUnitCell=2,
                                      correctVanadium=False)

    assert(np.isclose(factor,ds.absoluteNormalized))

def test_CustomAxisInput():
    # dset object generation
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),
                os.path.join(dataPath,'camea2018n000136.hdf')]
    ds = DataSet(DataFile)
    ds.convertDataFile()
    # parameter definition for plots
    cutqe_params = {
        "q1":np.array([0.75,0,0], float),
        "q2":np.array([1,0,0], float),
        "EMin":1.8, "EMax":3.6, "dE":0.1, "width":0.1, "minPixel":0.05,
        "rlu":True, "smoothing":0.0, "vmin":0.0, "vmax":1e-07
    }

    # custom axes generation
    fig, axes = plt.subplots(2,1) # two axes stacked on top of each other in same fig
    ax_direction1, ax_direction2 = axes.flatten()

    # plotting and cutting with MJOLNIR - fixed now for QE cuts
    q2 = [1,1,0]
    cutqe_params.update({"q2":q2, "ax":ax_direction1})
    ax_direction1, data_dir1, bins_dir1 = ds.plotCutQE(**cutqe_params)

    q2 = np.array([0.7,0.7,0],float)
    cutqe_params.update({"q2":q2, "ax":ax_direction2})
    ax_direction2, data_dir2, bins_dir2 = ds.plotCutQE(**cutqe_params)

    ## Same for plotCutQELine
    cutqe_params = {
        "QPoints":np.array([[0.75,0,0],[1,0,0]], float),
        'EnergyBins':np.arange(1.8,3.6,0.1), "width":0.1, "minPixel":0.05,
        "rlu":True, "vmin":0.0, "vmax":1e-07
    }

    # custom axes generation
    fig, axes = plt.subplots(2,1) # two axes stacked on top of each other in same fig
    ax_direction1, ax_direction2 = axes.flatten()

    # plotting and cutting with MJOLNIR - fixed now for QE cuts

    cutqe_params.update({"ax":ax_direction1})
    ax_direction1, data_dir1, bins_dir1 = ds.plotCutQELine(**cutqe_params)

    QPoints = np.array([[1,0,0],[1,1,0]], float)
    cutqe_params.update({'QPoints':QPoints, "ax":ax_direction2})
    ax_direction2, data_dir2, bins_dir2 = ds.plotCutQELine(**cutqe_params)

    ## cut1D
    cutqe_params = {
            "q1":np.array([0,0,0], float),
            "q2":np.array([1,0,0], float),
            "Emin":1.8, "Emax":2.0, "width":0.1, "minPixel":0.05,
            "rlu":True, 
        }

    # custom axes generation
    fig, axes = plt.subplots(2,1) # two axes stacked on top of each other in same fig
    ax_direction1, ax_direction2 = axes.flatten()

    # plotting and cutting with MJOLNIR - fixed now for QE cuts

    cutqe_params.update({"ax":ax_direction1})
    ax_direction1, data_dir1, bins_dir1 = ds.plotCut1D(**cutqe_params)

    q2 = np.array([0.7,0.7,0],float)
    cutqe_params.update({'q2':q2, "ax":ax_direction2})
    ax_direction2, data_dir2, bins_dir2 = ds.plotCut1D(**cutqe_params)


    # cut1DE
    cut1de_params = {
            "q":np.array([1,0,0], float),
            "E1":1.8, "E2":3.6, "width":0.5, "minPixel":0.05,
            "rlu":True, 
        }

    # custom axes generation
    fig, axes = plt.subplots(2,1) # two axes stacked on top of each other in same fig
    ax_direction1, ax_direction2 = axes.flatten()

    # plotting and cutting with MJOLNIR - fixed now for QE cuts

    cut1de_params.update({"ax":ax_direction1})
    ax_direction1, data_dir1, bins_dir1 = ds.plotCut1DE(**cut1de_params)

    q = np.array([0.7,0.7,0],float)
    cut1de_params.update({'q':q, "ax":ax_direction2})
    ax_direction2, data_dir2, bins_dir2 = ds.plotCut1DE(**cut1de_params)