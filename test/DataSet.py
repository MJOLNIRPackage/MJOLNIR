from importlib.util import decode_source
import numpy as np
import MJOLNIR.Data.DataFile
from MJOLNIR.Data.DataSet import DataSet,calculateGrid3D,binData3D,cut1DE,fmt,figureRowColumns,centeroidnp,compareNones,OxfordList, load, Dataset
from MJOLNIR import _tools

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os,sys
import warnings
import pytest
from MJOLNIR.Data import Mask



pythonVersion = sys.version_info[0]

dataPath = 'samlpedata'

@pytest.mark.unit
@pytest.mark.data
def test_Dataset_Initialization():

    emptyDataset = DataSet()
    del emptyDataset

    dataset = DataSet(dataFiles=[os.path.join(dataPath,'camea2018n000136.hdf')],calibrationfiles=[])
    
    assert(dataset.dataFiles[0].name=='camea2018n000136.hdf')
    assert(dataset.normalizationfiles == [])
    

@pytest.mark.unit
@pytest.mark.data                                                                                                        
def test_DataSet_Error():
    
    with pytest.raises(AttributeError):
        ds = DataSet(normalizationfiles=[10,11])


    ds = DataSet()
    
    with pytest.raises(AttributeError):
        ds.convertDataFile()

    with pytest.raises(AttributeError):
        ds.binData3D(0.1,0.1,0.1)

    with pytest.raises(AttributeError): 
        ds.cut1D([0,0],[1,1],0.1,0.01,5.5,6.0)

    with pytest.raises(AttributeError):
        ds.plotCut1D([0,0],[1,1],0.1,0.01,5.5,6.0)

    with pytest.raises(AttributeError):
        ds.cutQE([0,0],[1,1],0.1,0.01,5.5,6.0)


    with pytest.raises(AttributeError):
        ds.cutPowder(np.linspace(0,4,5),np.linspace(0,4,5))
        
    with pytest.raises(AttributeError): # No data files
        ds.plotCutPowder(np.linspace(0,4,5),np.linspace(0,4,5))
           
    

    with pytest.raises(AttributeError):
        ds.dataFiles = 100

    with pytest.raises(AttributeError):
        ds.convertedFiles = 10


    ds.dataFiles = os.path.join(dataPath,'camea2018n000136.hdf')


@pytest.mark.skip(reason="Bambus Data file currently not available")
@pytest.mark.integration
@pytest.mark.data
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
        np.testing.assert_allclose(getattr(ds[0],param)[0],value)


    with pytest.raises(AttributeError):
        ds.convertDataFile(binning = 3)
    
    assert(ds[0].I.shape == (scanValues.shape[1],100,1))

    ds.convertDataFile()


@pytest.mark.unit
@pytest.mark.data
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
        assert names[i]==df.name

    dataset.append(dataFiles)
    assert(len(dataset)==4)
    secondShape = dataset.Monitor.shape
    assert(np.all(secondShape!=initShape))
    del dataset[3]
    del dataset[2]
    with pytest.raises(IndexError):
        dataset[10]
    
    with pytest.raises(IndexError):
        del dataset[10]

    with pytest.raises(FileNotFoundError):
        dataset.append('NoFile')
    
    dataset.append(MJOLNIR.Data.DataFile.DataFile(dataFiles[0]))
    assert(len(dataset)==3)
    assert(dataset.I.shape!=secondShape)


@pytest.mark.integration
@pytest.mark.data
def test_DataSet_Equality():
    D1 = DataSet(dataFiles=os.path.join(dataPath,'camea2018n000136.hdf'))#,convertedFiles=['TestData/VanNormalization.nxs')])
    assert(D1==D1)


@pytest.mark.unit
@pytest.mark.data
def test_DataSet_str():
    D1 = DataSet(dataFiles=os.path.join(dataPath,'camea2018n000136.hdf'))#,normalizationfiles = 'TestData/VanNormalization.hdf'))
    string = str(D1)
    print(string)

@pytest.mark.integration
@pytest.mark.data
def test_DataSet_Convert_Data(): # TODO: redo test!
    dataFiles = os.path.join(dataPath,'camea2018n000136.hdf')
    dataset = DataSet(dataFiles=dataFiles)
    

    with pytest.raises(AttributeError):
        dataset.convertDataFile(dataFiles=dataFiles,binning=100)

    with pytest.raises(FileNotFoundError):
        dataset.convertDataFile(dataFiles='FileDoesNotExist',binning=1)

    with pytest.raises(FileNotFoundError):
        os.remove(os.path.join(dataPath,'camea2018n000136.nxs'))
    
    dataset.convertDataFile(dataFiles=dataFiles,binning=8,saveLocation=os.path.join(dataPath,''))
    convertedFile = dataset.convertedFiles[0]

    " other checks go here... "

    qx1 = np.asarray([1.81623706, 1.82088352, 1.82532258, 1.8304964,  1.83547466, 1.84185814,  1.83128525, 1.83757527])
    qx2 = np.asarray([1.14197299, 1.14631861, 1.15025765, 1.1461355,  1.15032864, 1.15418509,
           1.15816701, 1.16200206, 1.16592697, 1.16982047, 1.17596058, 1.16550963,
           1.17006764, 1.17433565, 1.17788609, 1.18175158, 1.18588861, 1.18948165,
           1.19378013, 1.18823503, 1.19288716, 1.19699238, 1.20105022])
    qy1 = np.asarray([-0.46624969, -0.47152159, -0.46554817, -0.45894287, -0.45253364])
    I = np.asarray([0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0])

    np.testing.assert_allclose(dataset.qx[0][5,99,2:10],qx1)
    np.testing.assert_allclose(dataset.qx[0][35,49,5:28],qx2)
    np.testing.assert_allclose(dataset.qy[0][35,39,55:60],qy1)
    np.testing.assert_allclose(dataset.I[0][10,4,50:64],I)


@pytest.mark.unit
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


@pytest.mark.unit
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
    assert(RebinnedNormCount.dtype==float)
    assert(RebinnedNorm.dtype==Norm.dtype)
    assert(ReBinnedI.dtype==I.dtype)

@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def test_DataSet_full_test():
    import MJOLNIR.Data.Viewer3D
    
    import matplotlib.pyplot as plt
    import os
    plt.ioff()
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf')]

    dataset = DataSet(dataFiles=DataFile)
    dataset.convertDataFile()
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

    os.remove('View3D.png')
    del viewer
    plt.close('all')

@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def test_DataSet_Visualization():
    import warnings
    from MJOLNIR.Data import Viewer3D
    DataFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]

    dataset = DataSet(dataFiles=DataFiles)
    dataset.convertDataFile()

    Data,bins = dataset.binData3D(0.08,0.08,0.25)
    
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    
    viewer = Viewer3D.Viewer3D(Intensity,bins)
    viewer.caxis = (0,100)
    with pytest.raises(AttributeError):
        viewer.caxis = 'Wrong type'
    
    with pytest.raises(AttributeError):
        viewer.caxis = [0,1,2,3,4] # Too long input
    with pytest.raises(AttributeError):
        viewer.setAxis(20) # Must bee 0,1, or 2

    plt.plot()
    plt.close('all')
    
@pytest.mark.unit
def test_DataSet_binEdges():
    X = np.random.rand(100)*3 # array between 0 and 3 -ish
    X.sort()
    tolerance = 0.01
    Bins = _tools.binEdges(X,tolerance=tolerance)

    assert(Bins[0]==X[0]-0.1*tolerance)
    assert(np.isclose(Bins[-1],X[-1],atol=5) or Bins[-1]>X[-1])
    assert(len(Bins)<=3.0/tolerance)
    assert(np.all(np.diff(Bins[:-1])>tolerance*0.99))

@pytest.mark.integration
@pytest.mark.data  
@pytest.mark.gui
def test_DataSet_1Dcut():
    q1 =  np.array([1.23,-1.51])
    q2 =  np.array([1.54, -1.25])
    width = 0.1

    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile()

    ax,Data,bins = ds.plotCut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5,fmt='.')
    Data2,bins2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5)
    
    # Check that the two data sets have the same values (except for Data2 also having 'BinDistance')
    assert(Data2.equals(Data.loc[:, Data.columns != 'BinDistance']))

    Data,bins = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5,extend=False)
    assert(np.all(np.logical_and(bins[0][:,0]>=q1[0]-0.1,bins[0][:,0]<=q2[0]+0.1))) 
    # x-values should be between 1.1 and 2.0 correpsonding to q points given (add some extra space due to way bins are created (binEdges))

    #q3 = np.array([1.1,1.1])
    #q4 = np.array([2.0,2.0])
    Data,bins = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5,extend=False)
    
    assert(np.all(bins[0][:,0]>=q1[0]-0.1))
    assert(np.all(bins[0][:,0]<=q2[0]+0.1))
    assert(np.all(bins[0][:,1]>=q1[1]-0.1))
    assert(np.all(bins[0][:,1]<=q2[1]+0.1))
    # x and y-values should be between 1.1 and 2.0 correpsonding to q points given (add some extra space due to way bins are created (binEdges))

    Q1 = np.array([1,0,0])
    Q2 = np.array([0.5,1,0])

    ax,Data,bins = ds.plotCut1D(Q1,Q2,width,rlu=True,minPixel=0.01,EMin=2.0,EMax=2.5,fmt='.')
    Data2,bins2 = ds.cut1D(Q1,Q2,width,rlu=True,minPixel=0.01,EMin=2.0,EMax=2.5)

    assert(Data2.equals(Data.loc[:,Data.columns!='BinDistance']))
    assert(np.all(np.array([np.all(np.isclose(bins[i],bins2[i])) for i in range(len(bins))]).flatten()))

    q1,q2 = ds.convertToQxQy([Q1,Q2])
    D1,b1 = ds.cut1D(Q1,Q2,width,rlu=True,minPixel=0.01,EMin=2.0,EMax=2.5)
    D2,b2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5)

    # Convert b1 to HKL in order to check if conversion works
    BinPos,OrthoPos,E = b1
    BinPos = np.concatenate([ds.convertToQxQy(BinPos[:,:3]),BinPos[:,-1].reshape(-1,1)],axis=1)
    OrthoPos = ds.convertToQxQy(OrthoPos)
    b1 = [BinPos,OrthoPos,E]

    
    assert(np.all(np.isclose(D1,D2)))
    
    assert(np.all(np.array([np.all(np.isclose(b1[i],b2[i])) for i in range(len(b1))]).flatten()))

    # Check that generating a plot from previous data is equivalent to directly plotting

    
    ax1,cut,bins = ds.plotCut1D(Q1,Q2,width=0.1,minPixel=0.04,EMin=2.0,EMax=2.5,ufit=False)

    ax2,*_ = ds.plotCut1D(Q1,Q2,EMin=2.0,EMax=2.5,width=0.1,minPixel=0.04,data=[cut,bins])



    for key in ['legend','title','xlabel','ylabel']:
        val1 = ax1.properties()[key]
        val2 = ax2.properties()[key]
        assert(val2==val1)

    # check plotted data in line plot 
    # ax1, find index of line2D
    id1 = np.arange(len(ax1.properties()['children']))[np.array([isinstance(child,mpl.lines.Line2D) for child in ax1.properties()['children']])][0]
    id2 = np.arange(len(ax2.properties()['children']))[np.array([isinstance(child,mpl.lines.Line2D) for child in ax2.properties()['children']])][0]

    np.testing.assert_array_equal(ax2.properties()['children'][id1]._xy,ax1.properties()['children'][id2]._xy)


@pytest.mark.skipif(Dataset is None, reason="ufit not installed")
@pytest.mark.ufit
@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def test_DataSet_1Dcut_ufit():
    q1 =  np.array([1.23,-1.51])
    q2 =  np.array([1.54, -1.25])
    width = 0.1

    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile()

    ax,dataset = ds.plotCut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5,fmt='.',ufit=True)
    dataset2 = ds.cut1D(q1,q2,width,rlu=False,minPixel=0.01,EMin=2.0,EMax=2.5,ufit=True)
    

    files = ', '.join([x.replace('hdf','nxs').split(os.path.sep)[-1] for x in convertFiles])

    assert(np.all([np.all(x==y) for x,y in zip(dataset.fit_columns,dataset2.fit_columns)]))
    assert(dataset.meta == dataset2.meta)

    assert(dataset.meta['instrument'] == 'CAMEA')
    assert(dataset.meta['datafilename'] == files)

    
@pytest.mark.skipif(Dataset is None, reason="ufit not installed")
@pytest.mark.ufit
@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def test_DataSet_1DcutE():
    q =  np.array([1.23,-1.25]).reshape(2,1)
    width = 0.1
    EMin = 1.5
    EMax = 2.5
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000137.hdf')]
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    Datset._getData()
    I,qx,qy,energy,Norm,Monitor = Datset.I.extractData(),Datset.qx.extractData(),Datset.qy.extractData(),Datset.energy.extractData(),Datset.Norm.extractData(),Datset.Monitor.extractData()

    [intensity,MonitorCount,Normalization,normcounts],[bins] = cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=EMin,E2=EMax,q=q,width=width,minPixel=0.01)
    Q = Datset.convertToHKL(q.reshape(2))
    
    Data,[bins] = Datset.cut1DE(E1=EMin,E2=EMax,q=Q,width=width,minPixel=0.01)
    assert(np.min(bins)>=EMin-0.01) # Check that bins do not include data outside of cut
    assert(np.max(bins)<=EMax+0.01)
    assert(len(bins[0])==len(intensity)+1)# Bins denotes edges and must then be 1 more than intensity

    assert(intensity.shape==MonitorCount.shape) # Check that all matrices are cut equally
    assert(intensity.shape==Normalization.shape)
    assert(intensity.shape==normcounts.shape)

    Data,[bins] = Datset.cut1DE(E1=EMin,E2=EMax,q=q,width=width,minPixel=0.01,rlu=False)
    
    Data,[bins] = Datset.cut1DE(E1=EMin,E2=EMax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    
    np.testing.assert_allclose(np.diff(bins),0.01)
    assert(bins[0].min()>=EMin)
    assert(bins[0].max()<=EMax)

    with pytest.raises(AttributeError):
        cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=500,E2=700,q=q,width=width,minPixel=0.01)

    with pytest.raises(AttributeError):
        cut1DE(positions=[qx,qy,energy],I=I,Norm=Norm,Monitor=Monitor,E1=5,E2=7,q=np.array([20.0,0]).reshape(2,1),width=width,minPixel=0.01)

    # Check the data of plot to be the same as cut
    Data,[bins] = Datset.cut1DE(E1=EMin,E2=EMax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    ax,Data2,[bins2] = Datset.plotCut1DE(E1=EMin,E2=EMax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True)
    assert(Data.equals(Data2))
    np.testing.assert_allclose(bins[0],bins2[0])

    ufitData = Datset.cut1DE(E1=EMin,E2=EMax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True,ufit=True)
    ax,ufitData2 = Datset.plotCut1DE(E1=EMin,E2=EMax,q=q,width=0.1,minPixel=0.01,rlu=False,constantBins=True,ufit=True)

    files = ', '.join([x.replace('hdf','nxs').split(os.path.sep)[-1] for x in convertFiles])
    
    assert(np.all([np.all(np.isclose(x,y,equal_nan=True)) for x,y in zip(ufitData.fit_columns,ufitData2.fit_columns)]))
    assert(ufitData.meta == ufitData2.meta)

    assert(ufitData.meta['instrument'] == 'CAMEA')
    assert(ufitData.meta['datafilename'] == files)

    ax,Data3,[bins3] = Datset.plotCut1DE(E1=EMin,E2=EMax,q=Q,width=0.1,minPixel=0.01,constantBins=True)


@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
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

    Datset.convertDataFile()
    ax,Data,bins = Datset.plotCutQE(q1,q2,width=width,minPixel=minPixel,EnergyBins=EnergyBins,rlu=False)# Remove to improve test coverage ,vmin=0.0 , vmax= 5e-06)
    Data2,bins2,_ = Datset.cutQE(q1,q2,width=width,minPixel=minPixel,EnergyBins=EnergyBins,rlu=False)

    comparisons = list(Data.columns)
    del comparisons[comparisons.index('BinDistance')]

    assert(np.all([np.allclose(Data[p],Data2[p],equal_nan=True) for p in comparisons]))

    Q1 = np.array([1,0,0])
    Q2 = np.array([0.5,1,0])

    q1,q2 = Datset.convertToQxQy([Q1,Q2])

    Data1,bins,_ = Datset.cutQE(Q1,Q2,width,minPixel,EnergyBins=EnergyBins,rlu=True)
    Data2,bins2,_ = Datset.cutQE(q1,q2,width,minPixel,EnergyBins=EnergyBins,rlu=False)

    comparisons = list(Data.columns)
    del comparisons[comparisons.index('BinDistance')]
    
    assert(np.all(comparisons == Data2.columns))

    
    assert(np.all([np.allclose(Data1[p],Data2[p],equal_nan=True) for p in comparisons]))

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataSet_cutPowder():

    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()
    mask = Datset.mask#np.ones_like(Datset.mask)[0]

    Datset.mask = mask
    Datset.mask = np.logical_not(mask)
    

    eBins = _tools.binEdges(Datset.energy,0.25)
    qBins = np.linspace(0.0,4,201)

    ax,D,QCentres,ECentres = Datset.plotCutPowder(eBins,qBins)# Remove to improve test ,vmin=0,vmax=1e-6)
    D2 = Datset.cutPowder(eBins,qBins)
    assert(np.all(D.equals(D2)))


@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataSet_createRLUAxes():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    fig = plt.figure()
    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile()

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

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataSet_createQEAxes():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000136.hdf')]
    
    ds = DataSet(dataFiles = convertFiles)
    ds.convertDataFile()

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


@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
def test_DataSet_plotQPlane():
    plt.ioff()
    import matplotlib
    matplotlib.use('Agg')

    convertFiles = [os.path.join(dataPath,'camea2018n000137.hdf')]#'TestData/ManuallyChangedData/A3.hdf')]
    
    Datset = DataSet(dataFiles = convertFiles)
    Datset.convertDataFile()

    EmptyDS = DataSet()
    with pytest.raises(AttributeError):
        Datset.plotQPlane() # No Bins, EMin or EMax
        
    with pytest.raises(AttributeError):
        Datset.plotQPlane(EBins=[10]) # Length of bins is 1

    with pytest.raises(AttributeError):
        Datset.plotQPlane(EMin=20,EMax=10) # EMin>EMax
    
    with pytest.raises(AttributeError):
        EmptyDS.plotQPlane(EMin=2,EMax=3) # Empty DataSet
    
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
    
    with pytest.raises(AttributeError):
        Datset.plotQPlane(EMin,EMax,binning='notABinningMethod')

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


    
@pytest.mark.unit
def test_DataSet_fmt():
    assert('$1.00 \\times 10^{1}$' == fmt(10,'Unused'))
    assert('$1.00 \\times 10^{-10}$' == fmt(1e-10,'Unused'))
    assert('$2.55 \\times 10^{-2}$' == fmt(0.0255,'Unused'))
    assert('$2.56 \\times 10^{-2}$' == fmt(0.02556,'Unused'))
    
@pytest.mark.unit
def test_DataSet_figureRowColumns():
    assert(np.all(np.array([3,4])==np.array(figureRowColumns(10)))) # 10 -> 3,4
    assert(np.all(np.array([3,3])==np.array(figureRowColumns(9)))) # 9 -> 3,3
    assert(np.all(np.array([1,1])==np.array(figureRowColumns(1)))) # 1 -> 1,1

    with pytest.raises(AttributeError):
        figureRowColumns(0) # No figures
        
    assert(np.all(np.array([8,8])==np.array(figureRowColumns(63)))) # 63 -> 8,8
    
@pytest.mark.unit
def test_DataSet_centeroidnp():
    pos = np.array([[0,0],[1,0],[0,1],[1,1]],dtype=float)
    np.testing.assert_allclose(np.array([0.5,0.5]),centeroidnp(pos))

    pos2 = np.array([[1.2,2.2],[7.5,1.0],[11.0,0.0],[4.0,-1.0],[2.0,2.0]],dtype=float)
    np.testing.assert_allclose(np.array([5.14,0.84]),centeroidnp(pos2))

@pytest.mark.unit
def test_DataSet_compareNones():
    assert(compareNones(np.array([None]),np.array([None]),0.1))
    assert(not compareNones(np.array([None]),np.array([0.5]),0.1))
    assert(not compareNones(np.array([0.5]),np.array([None]),0.1))
    assert(compareNones(np.array([0.4]),np.array([0.5]),0.2))
    assert(not compareNones(np.array([0.4]),np.array([0.5]),0.001))

    assert(not np.all(compareNones(np.array([0.4,10.2,10.0]),np.array([0.5]),0.001)))
    assert(np.all(compareNones(np.array([0.4,10.2,10.0]),np.array([0.4,10.2,10.0]),0.001)))

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
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
    dataset.convertDataFile()

    with pytest.raises(AttributeError):
        dataset.cutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu=True)
    with pytest.raises(AttributeError):
        dataset.cutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu='42') # Wrong RLU-input

    DataList,BinList=dataset.cutQELine(QPointsHKL,EnergyBins,width=width,minPixel=minPixel,rlu=True)
    DataList2,BinList2=dataset.cutQELine(QPoints,EnergyBins,width=width,minPixel=minPixel,rlu=False)
    print(DataList)
    assert(len(DataList)==(len(QPointsHKL)-1))# Assert that there are 3 cuts 

    assert(len(DataList2)==(len(QPoints)-1))

    
    
@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
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
    dataset.convertDataFile()
    
    with pytest.raises(AttributeError):
        dataset.plotCutQELine([],EnergyBins,width=width,minPixel=minPixel,rlu=False)


    with pytest.raises(AttributeError):
        dataset.plotCutQELine(QPoints,EnergyBins,width=width,minPixel=minPixel,rlu=False)

    with pytest.raises(AttributeError):
        dataset.plotCutQELine(QPoints[:,:2],EnergyBins,width=width,minPixel=minPixel,rlu=True)

    ax,Data = dataset.plotCutQELine(
        QPoints=QPoints[:,:2],EnergyBins=EnergyBins,width=width,minPixel=minPixel,rlu=False,vmin=0.0,vmax=1.5e-6,log=True,seperatorWidth=3)


    HKLPoints = np.array([[1.0,0.0,0.0],
                        [0.5,1.5,0.0],
                        [1.7,-0.1,0.0],
                        [1.0,1.0,0.0]])



    ax,Data = dataset.plotCutQELine(
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
        
        with pytest.raises(NotImplementedError):
            ax,Data,Bins = \
            dataset.plotCutQELine(QPoints=HKLPoints,EnergyBins=E,ax = ax,width=0.05,minPixel=0.01,
                    vmin=7.5e-7,vmax=7e-6,cmap=cmap,rlu=True)


@pytest.mark.data
@pytest.mark.integration
@pytest.mark.slow
def test_DataSet_extractDetectorData():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]#['TestData/ManuallyChangedData/A3.nxs'),'TestData/ManuallyChangedData/A3.nxs')]
    dataset = DataSet(DataFile)

    binning = 1
    dataset.convertDataFile(binning=binning)

    with pytest.raises(AttributeError):
        dataset.extractDetectorData(A4=10000.0) # A4 outside of detector

    with pytest.raises(AttributeError):
        dataset.extractDetectorData(Ef=10000.0) # Ef outside of detector
    
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
        
@pytest.mark.data
@pytest.mark.integration
@pytest.mark.slow
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
    np.testing.assert_allclose(subtracted.Norm.data,dataset.Norm.data)
    np.testing.assert_allclose(subtracted.Monitor.data,dataset.Monitor.data)

@pytest.mark.unit
def test_DataSet_OxfordList():
    l = ['Apples','Pears']
    S = OxfordList(l)
    assert(S=='Apples and Pears')

    l.append('Oranges')
    S = OxfordList(l)
    assert(S=='Apples, Pears, and Oranges')

    assert(OxfordList([]) is None)
    assert(OxfordList(['Apples'])=='Apples')

@pytest.mark.skip(reason="Current data structure not defined")
@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def test_DataSet_MultiFLEXX():
    fileLocation = _tools.fileListGenerator('65059',folder=os.path.join('Data',''),instrument='MultiFLEXX')

    ds = DataSet(fileLocation)
    ds.convertDataFile()
    import matplotlib
    matplotlib.use('Agg')

    V = ds.View3D(0.05,0.05,0.5,grid=True)

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def test_DataSet_ELine():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)

    dataset.convertDataFile()
    import matplotlib
    matplotlib.use('Agg')

    Q1 = [1.0,-0.185,0.0]
    Q2 = [0.5,1.5,0.0]
    EMin = 1.65
    EMax = 3.3

    CutData,Bins = dataset.cutELine(Q1, Q2, EMin=EMin, EMax=EMax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=True, dataFiles=None, constantBins=False)
    ax, CutDataPlot, BinsPlot = dataset.plotCutELine(Q1, Q2, EMin=EMin, EMax=EMax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=True, dataFiles=None, constantBins=False)

    assert(np.all([np.all(np.isclose(B,B2)) for B,B2 in zip(Bins,BinsPlot)]))
    np.testing.assert_allclose(CutDataPlot, CutData)

    assert(np.all(np.logical_and(CutData['Energy']>=EMin,CutData['Energy']<=EMax*1.01))) # Allow for slightly heigher energy
    assert(np.logical_and(np.all(CutData['H']<=Q1[0]*1.01),np.all(CutData['H']>=Q2[0]*0.99)))
    assert(np.logical_and(np.all(CutData['K']>=Q1[1]*1.01),np.all(CutData['K']<=Q2[1]*1.01)))
    np.testing.assert_allclose(CutData['L'],0.0,atol=1e-6)
    assert(np.all([np.all(np.logical_and(B[0]>=EMin*0.99,B[0]<=EMax*1.05)) for B in Bins])) # Allow for slightly heigher energy


    if sys.version[0] == '3':
        ax.set_xticks_base(0.01)
        ax.set_xticks_number(10)
        ax.set_xticks_base()
        

    Q1raw = dataset.convertToQxQy(Q1)
    Q2raw = dataset.convertToQxQy(Q2)

    CutData,Bins = dataset.cutELine(Q1raw, Q2raw, EMin=EMin, EMax=EMax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=False)
    ax, CutDataPlot, BinsPlot = dataset.plotCutELine(Q1raw, Q2raw, EMin=EMin, EMax=EMax, energyWidth = 0.05, minPixel = 0.02, width = 0.02, rlu=False)


    assert(np.all([np.all(np.isclose(B,B2)) for B,B2 in zip(Bins,BinsPlot)]))
    assert(np.all(CutDataPlot.equals(CutData)))

    assert(np.all(np.logical_and(CutData['Energy']>=EMin,CutData['Energy']<=EMax*1.01))) # Allow for slightly heigher energy
    assert(np.logical_and(np.all(CutData['Qx']<=Q1raw[0]*1.01),np.all(CutData['Qx']>=Q2raw[0]*0.99)))
    assert(np.logical_and(np.all(CutData['Qy']<=Q1raw[1]*0.99),np.all(CutData['Qy']>=Q2raw[1]*1.01)))

    assert(np.all([np.all(np.logical_and(B[0]>=EMin*0.99,B[0]<=EMax*1.05)) for B in Bins])) # Allow for slightly heigher energy

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.slow
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

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def testplotRaw1D_Error():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
     # Scan variables are A3 and A3+A4
    dataset = DataSet(DataFile)
    dataset[0].scanParameters = ['Ei']
    import matplotlib
    matplotlib.use('Agg')

    with pytest.raises(AttributeError):
        dataset.plotRaw1D() # Two different scan types
        
    DataFile = [os.path.join(dataPath,'camea2018n000137.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)
    with pytest.raises(AttributeError):
        dataset.plotRaw1D(analyzerSelection=[0]) # Not enough analyzers provided

    with pytest.raises(AttributeError):
        dataset.plotRaw1D(detectorSelection=[0]) # Not enough detectors provided

        
    with pytest.raises(AttributeError):
        dataset.plotRaw1D(legend=[0]) # Not enough legend labels provided
    

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
@pytest.mark.slow
def testplotRaw1D():
    DataFile = [os.path.join(dataPath,'camea2018n000137.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    dataset = DataSet(DataFile)
    import matplotlib
    matplotlib.use('Agg')

    ax = dataset.plotRaw1D()
    
    ax = dataset.plotRaw1D(legend=['1','2'],detectorSelection=[0,0],analyzerSelection=[5,5],grid=True)
    assert True

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.slow
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
    with pytest.raises(AttributeError):
        ds.mask = 'Error'
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mask = np.array([np.ones_like(I) for I in ds.I]) # Mask all
        ds.mask = mask

    ds = DataSet([DataFile[0],DataFile[0]])
    ds.convertDataFile()

    
    print([x.shape for x in mask])
    print(mask[0].shape)
    
    ds.mask = mask[0] # Only provide 1 mask to be applied to both data files
    
    
@pytest.mark.data
@pytest.mark.integration
def testupdateSampleParameters():
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    ds = DataSet(DataFile)

    unitCell = ds.sample[0].unitCell
    unitCell[:2]*=2.0
    ds.updateSampleParameters(unitCell=unitCell)
    newCell = ds.sample[0].cell
    newUB = ds.sample[0].UB

    for d in ds:
        np.testing.assert_allclose(d.sample.cell, newCell)
        np.testing.assert_allclose(d.sample.UB, newUB)

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.slow
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
    

@pytest.mark.data
@pytest.mark.integration
def test_absoluteNormalziation(): # TODO: Improve and recalculate the expected factor for MnF2
    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf')]
    ds = DataSet(DataFile)

    with pytest.raises(AttributeError): # Must be converted first!
        ds.absoluteNormalize(10.0,'MnF2')

    ds.convertDataFile()
    norm = np.mean(ds.Norm.extractData())
    
    # Use value for MnF2 to check 
    ds.absoluteNormalize(sampleMass=6.2,sampleChemicalFormula='MnF2',formulaUnitsPerUnitCell=2,
                                      correctVanadium=False)
    
    factor = 2*0.06088201383247563 # Factor calculated for MnF2

    np.testing.assert_allclose(factor,ds.absoluteNormalized)

    # Redo normalization to retrive the same factor
    ds.absoluteNormalize(sampleMass=6.2,sampleChemicalFormula='MnF2',formulaUnitsPerUnitCell=2,
                                      correctVanadium=False)

    np.testing.assert_allclose(factor,ds.absoluteNormalized)

@pytest.mark.data
@pytest.mark.integration
@pytest.mark.gui
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
    ax_direction1, data_dir1 = ds.plotCutQELine(**cutqe_params)

    QPoints = np.array([[1,0,0],[1,1,0]], float)
    cutqe_params.update({'QPoints':QPoints, "ax":ax_direction2})
    ax_direction2, data_dir2 = ds.plotCutQELine(**cutqe_params)

    ## cut1D
    cutqe_params = {
            "q1":np.array([0,0,0], float),
            "q2":np.array([1,0,0], float),
            "EMin":1.8, "EMax":2.0, "width":0.1, "minPixel":0.05,
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

@pytest.mark.data
@pytest.mark.integration
def test_symmetrize():

    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),]
    ds = DataSet(DataFile)
    ds.convertDataFile()

    qxOriginal = ds[0].qx
    qyOriginal = ds[0].qy

    def identity(H,K,L):
        return H,K,L

    ds.symmetrize(identity)

    qx = ds[0].qx
    qy = ds[0].qy


    np.testing.assert_allclose(qx, qxOriginal, atol=1e-6)
    np.testing.assert_allclose(qy, qyOriginal, atol=1e-6)

    def absolute(H,K,L):
        return np.abs(H),np.abs(K),np.abs(L)

    ds.symmetrize(absolute)

    assert(np.all(ds[0].h.min()>=0))
    assert(np.all(ds[0].k.min()>=0))
    assert(np.all(ds[0].l.min()>=0))


@pytest.mark.data
@pytest.mark.integration
def test_sorting():

    DataFile = [os.path.join(dataPath,'camea2018n000136.hdf'),os.path.join(dataPath,'camea2018n000137.hdf')]
    ds = DataSet(DataFile)
    ds.convertDataFile()

    ds2 = DataSet([DataFile[1],DataFile[0]])
    ds2.convertDataFile()

    names = [df.name for df in ds]
    names2 = [df.name for df in ds2]

    assert(names[0] == names2[1] and names[1] == names2[0])

    idx = ds.argAutoSort()
    ds.autoSort()

    for idx1,df in zip(idx,ds):
        oldName = names[idx1]
        newName = df.name
        assert(oldName==newName)

    idxes = ds2.argSortTo(ds)
    ds2.applySorting(idxes)
    names = [df.name for df in ds]
    names2 = [df.name for df in ds2]

    assert(names[0] == names2[0] and names[1] == names2[1])

    
