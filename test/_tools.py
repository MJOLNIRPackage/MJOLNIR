import numpy as np
from MJOLNIR._tools import *

import os
import pytest

@pytest.mark.unit
def test_minMax():
    L = np.random.rand(10,3,2)
    minmax = minMax(L)
    
    assert(len(minmax)==2)
    np.testing.assert_allclose(minmax[0],np.min(L))
    np.testing.assert_allclose(minmax[1],np.max(L))

    minmax = np.array(minMax(L,axis=-1))
    assert(minmax.shape==(2,10,3))

    minmax = np.array(minMax(L,axis=0))
    assert(minmax.shape==(2,3,2))

@pytest.mark.unit
def test_unitVector():
    V = np.array([1.0,2.0,3.0])
    v = unitVector(V)
    np.testing.assert_allclose(np.linalg.norm(v),1.0)


@pytest.mark.unit
def test_rotations():
    vectors = [np.array([0,0,3.0]),np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0]),np.random.rand(3)]
    rotations = [rotate2X(v) for v in vectors]
    rotVector = [np.dot(rotations[i],vectors[i]) for i in range(len(vectors))]
    for i in range(len(rotVector)):
        np.testing.assert_allclose(np.linalg.norm(vectors[i]),np.linalg.norm(rotVector[i]))
        np.testing.assert_allclose(rotVector[i][0],np.linalg.norm(rotVector[i]))

@pytest.mark.unit
def test_vectorAngle():
    v1 = np.array([1,0,0])
    v2 = np.array([0,1,0])
    v3 = np.array([1,1,1])
    theta1 = vectorAngle(v1,v2)
    theta2 = vectorAngle(v1,v3)
    theta3 = vectorAngle(v2,v3)
    print(theta2)
    np.testing.assert_allclose(theta1,np.pi/2)
    np.testing.assert_allclose(theta2,theta3)
    np.testing.assert_allclose(theta2,0.955316618125)

@pytest.mark.unit
def test_Rotation_matrix():
    M1 = rotationMatrix(0,0,0)

    np.testing.assert_allclose(M1,np.identity(3))

    M2 = rotationMatrix(90,0,0)
    M3 = rotationMatrix(np.pi/2,np.pi/3,np.pi/4,format='rad')
    M4 = rotationMatrix(90,60,45)
    
    np.testing.assert_allclose(M2,np.array([[1,0,0],[0,0,-1],[0,1,0]]),atol=1e-8)
    np.testing.assert_allclose(M3,M4)
    
    M4_check = np.array([[3.53553391e-01,6.12372436e-01,7.07106781e-01],[3.53553391e-01,6.12372436e-01,-7.07106781e-01],[-8.66025404e-01,5.00000000e-01,3.06161700e-17]])
    np.testing.assert_allclose(M4,M4_check,atol=1e-8)



@pytest.mark.unit
def test_binEdges():
    values = np.exp(np.linspace(-0.1,1,101)) # Generate non-linear points
    #values[0]-tolerance/2.0 and ends at values[-1]+tolerance/2.0.
    minBin = 0.1
    bins = binEdges(values,minBin)

    assert(np.isclose(bins[0],values[0]-0.1*minBin)) # First bin starts at values[0]-tolerance*0.1
    assert(bins[-1]>=values[-1]+0.1*minBin) # Last bin ends at values[-1]+tolerance*0.1
    assert(np.all(np.diff(bins)>=minBin*0.99)) # Assert that all bins are at least of size minBin

    binsCut = binEdges(values,0.01,startPoint=0.0)
    assert(binsCut[0]>=0.0)

    binsCut = binEdges(values+np.min(values)+0.0099,0.01,startPoint=0.0,endPoint=0.5)
    assert(binsCut[0]>=0.0)
    assert(binsCut[-1]<=0.5)

    binsCut = binEdges(values+np.min(values)+0.0101,0.01,startPoint=0.0,endPoint=0.5)
    assert(binsCut[0]>0.01)

    binsCut = binEdges(np.linspace(0,0.95,101),0.01,startPoint=0.0,endPoint=1.0101)
    assert(binsCut[-1]<=1.0)

    binsCut = binEdges(values-np.max(values)+1.0,0.01,startPoint=0.0,endPoint=1.01)
    assert(binsCut[-1]<=1.01)

@pytest.mark.unit
def test_fileListGenerator():
    numberStr = '0,20-23-24,4000'
    year = 2018
    folder = os.path.join('','home','camea')
    with pytest.raises(AttributeError):
        files = fileListGenerator(numberStr,folder,year)

    numberStr = '0,20-23,4000'

    files = fileListGenerator(numberStr,folder,year)
    filesCorrect = np.array([
        os.path.join('','home','camea','camea2018n000000.hdf'), os.path.join('','home','camea','camea2018n000020.hdf'), 
        os.path.join('','home','camea','camea2018n000021.hdf'), os.path.join('','home','camea','camea2018n000022.hdf'), 
        os.path.join('','home','camea','camea2018n000023.hdf'), os.path.join('','home','camea','camea2018n004000.hdf')])
    assert(np.all(filesCorrect == np.array(files)))

    # opposite direction
    numberStr = '23-20'
    filesCorrect = np.array([
        os.path.join('','home','camea','camea2018n000023.hdf'),os.path.join('','home','camea','camea2018n000022.hdf'), 
        os.path.join('','home','camea','camea2018n000021.hdf'),os.path.join('','home','camea','camea2018n000020.hdf')])

    files = fileListGenerator(numberStr,folder,year)
    assert(np.all(filesCorrect == np.array(files)))

@pytest.mark.unit
def test_RoundBinning():
    
    X = np.random.rand(2,3000)
    binning = [0.05,0.1,0.2]
    try:
        RoundBinning(X,binning)
    except AttributeError:
        assert True
    
    binning = [0.05,0.05]
    BP, I, C  = RoundBinning(X,binning)

    Int = np.random.rand(3000)
    BP2,I2,C2,Data = RoundBinning(X,binning[0],Int)

    assert(np.all(BP==BP2))
    assert(np.all(I==I2))
    assert(np.all(C==C2))
    assert(np.all(Data.shape[0]==BP.shape[1]))

@pytest.mark.unit
def test_generateLabel():
    v = [1,-1,2]
    label = generateLabel(v)
    assert(label=='(H, -K, 2L)')

    v = np.array([-0.5,22.0000,-3.0])
    label = generateLabel(v)
    assert(label=='(-0.5H, 22K, -3L)')

@pytest.mark.unit
def test_generateLabelDirection():
    v = [1,-1,2]
    label = generateLabelDirection(v)
    assert(label=='(H, -H, 2H)')

    v = np.array([0.0,22.0000,-3.0])
    label = generateLabelDirection(v)
    assert(label=='(0, 22K, -3K)')

    v = np.array([0.0,0.0000,-3.333])
    label = generateLabelDirection(v)
    assert(label=='(0, 0, -3.333L)')

@pytest.mark.unit
def test_molarMassCalculation():
    water = 18.01528 # g/mol
    value,elements = calculateMolarMass('H2O',returnElements=True)
    assert(np.isclose(value,water,atol=1e-4))
    assert(elements=={'H':2,'O':1})

    glucose = 180.156 # g/mol
    assert(np.isclose(calculateMolarMass('C6H12O6'),glucose,atol=1e-4))

    YMnO3 = 191.84 # g/mol
    assert(np.isclose(calculateMolarMass('YMnO3'),YMnO3,atol=1e-4))

    YBCO = 666.19 # g/mol
    value = calculateMolarMass('YBa2Cu3O7')
    assert(np.isclose(value,YBCO,atol=1e-4))
    
    
    Crazy = 'U2.3(H3O2.2)2.2Be((He2.2Ge)2Fe)3.3'
    weight = 2724.6618342799998
    elem = {'H': 13.2,  'O': 9.68,  'He': 29.04,  'Ge': 13.2,  'Fe': 6.6,  'U': 4.6,  'Be': 2.0}
    value,elements = calculateMolarMass(Crazy,formulaUnitsPerUnitCell=2,returnElements=True)
    assert(np.isclose(value,weight,atol=1e-4))
    assert(np.all([np.isclose(elem[key],elements[key],atol=1e-8) for key in elem.keys()]))

@pytest.mark.unit
def test_absoluteNormalization():# TODO: Improve and recalculate the expected factor 
    sampleMass = 6.2
    normFactor = calculateAbsoluteNormalization(sampleChemicalFormula='MnF2',formulaUnitsPerUnitCell=2,sampleMass=sampleMass,correctVanadium=True)
    # Known value for MnF2
    assert(np.isclose(normFactor,2*2.0016916920401816e-07))

    vanadiumMass=15.25
    vanadiumMonitor=1e5
    vanadiumSigmaIncoherent = 5.08
    constants = 4*np.pi/13.77
    Van = calculateAbsoluteNormalization(sampleChemicalFormula='V',formulaUnitsPerUnitCell=1,sampleMass=vanadiumMass,correctVanadium=True)
    

    Van*=vanadiumMonitor*vanadiumSigmaIncoherent/constants
    assert(np.isclose(Van,2*1.0))

    # Same as above but use the molecular mass in stead of calculated
    sampleMass = 6.2
    normFactor2 = calculateAbsoluteNormalization(sampleMolarMass=92.9349*2,sampleMass=sampleMass,correctVanadium=True)
    # Known value calcualted above
    assert(np.isclose(normFactor,normFactor2))


@pytest.mark.unit
def test_KWavelength():
    k = np.pi
    wavelength = WavelengthK(k)
    assert(np.isclose(wavelength,2.0))
    
    k_prime = KWavelength(wavelength)
    assert(np.isclose(k,k_prime))

@pytest.mark.unit
def test_KWavelength_multidim():
    k = np.random.rand(10,20)*2*np.pi
    wavelength = WavelengthK(k)
    
    k_prime = KWavelength(wavelength)
    np.testing.assert_allclose(k,k_prime)
    
@pytest.mark.unit
def test_EnergyK():
    E = 5.0
    k = KEnergy(E)
    l = WavelengthEnergy(E)
    assert(np.isclose(l,4.044851376460777))
    
    E_prime = EnergyK(k)
    assert(np.isclose(E,E_prime))

@pytest.mark.unit
def test_WnergyK_multidim():
    E = np.random.rand(10,20)*2*np.pi
    wavelength = WavelengthEnergy(E)
    
    E_prime = EnergyWavelength(wavelength)
    np.testing.assert_allclose(E,E_prime)
    
    
@pytest.mark.unit
def test_ScatteringAngle():
    d = 3.355 # AA
    E = 5.00 # meV
    K = KEnergy(E)
    W = WavelengthEnergy(E)
    
    A4E = ScatteringAngle(d,Energy=E)
    A4K = ScatteringAngle(d,K=K)
    A4W = ScatteringAngle(d,Wavelength=W)
    A4ERad = ScatteringAngle(d,Energy=E,degrees=False)
    A4ERad = np.rad2deg(A4ERad)
    
    A4True = 74.14270982536185
    
    np.testing.assert_allclose([A4E,A4K,A4W,A4ERad],A4True,atol=1e-5)

@pytest.mark.unit
def test_DSpacing():
    twoTheta = 74.14270982536185#
    E = 5.00 # meV
    K = KEnergy(E)
    W = WavelengthEnergy(E)
    
    dE = DSpacing(twoTheta,Energy=E)
    dK = DSpacing(twoTheta,K=K)
    dW = DSpacing(twoTheta,Wavelength=W)
    dERad = DSpacing(np.deg2rad(twoTheta),Energy=E,degrees=False)
    
    
    dTrue = 3.355 # AA
    
    assert(np.all(np.isclose([dE,dK,dW,dERad],dTrue)))
    
@pytest.mark.unit
def test_ScatteringAngle_Errors():
    d = 3.355 # AA
    E = 5.00 # meV
    K = KEnergy(E)
    W = WavelengthEnergy(E)
    
    with pytest.raises(AttributeError):
        _ = ScatteringAngle(d,Energy=E,K=K)
    
    with pytest.raises(AttributeError):
        _ = ScatteringAngle(d,Energy=E,Wavelength=W)
    
    with pytest.raises(AttributeError):
        _ = ScatteringAngle(d,Energy=E,Wavelength=W,K=K)
    
    with pytest.raises(AttributeError):
        _ = ScatteringAngle(d)
        
@pytest.mark.unit
def test_DSpacing_Errors():
    twoTheta = 3.355 # AA
    E = 5.00 # meV
    K = KEnergy(E)
    W = WavelengthEnergy(E)
    
    with pytest.raises(AttributeError):
        _ = DSpacing(twoTheta,Energy=E,K=K)
    
    with pytest.raises(AttributeError):
        _ = DSpacing(twoTheta,Energy=E,Wavelength=W)
    
    with pytest.raises(AttributeError):
        _ = DSpacing(twoTheta,Energy=E,Wavelength=W,K=K)
    
    with pytest.raises(AttributeError):
        _ = DSpacing(twoTheta)
    
@pytest.mark.unit
def test_scaling():# Test standard scaling functions used for rescaling of energy bins

    x = np.random.rand(20)*20-10 # points between -10 and 10
    np.testing.assert_allclose(x,rescale(scale(x)))
    np.testing.assert_allclose(x,rescale(scale(x,f=0.1),f=0.1))

def test_getNext():
    
    data = [1,2,3,4,5,6,7,8]

    assert(len(data)==8)
        
        
    for d in getNext(data,delete=True):
        print(d)
    assert(len(data)==0)
    
    data = [1,2,3,4,5,6,7,8]
    
    for d in getNext(data,delete=False):
        print(d)
    
    assert(len(data)==8)
    