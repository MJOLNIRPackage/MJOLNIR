import numpy as np
from MJOLNIR._tools import rotate2X, minMax, unitVector, vectorAngle, rotationMatrix, binEdges, fileListGenerator, RoundBinning, generateLabel, generateLabelDirection

import os

def test_minMax():
    L = np.random.rand(10,3,2)
    minmax = minMax(L)
    
    assert(len(minmax)==2)
    assert(np.isclose(minmax[0],np.min(L)))
    assert(np.isclose(minmax[1],np.max(L)))

    minmax = np.array(minMax(L,axis=-1))
    assert(minmax.shape==(2,10,3))

    minmax = np.array(minMax(L,axis=0))
    assert(minmax.shape==(2,3,2))

def test_unitVector():
    V = np.array([1.0,2.0,3.0])
    v = unitVector(V)
    assert(np.isclose(np.linalg.norm(v),1.0))


def test_rotations():
    vectors = [np.array([0,0,3.0]),np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0]),np.random.rand(3)]
    rotations = [rotate2X(v) for v in vectors]
    rotVector = [np.dot(rotations[i],vectors[i]) for i in range(len(vectors))]
    for i in range(len(rotVector)):
        assert(np.isclose(np.linalg.norm(vectors[i]),np.linalg.norm(rotVector[i])))
        print(rotVector[i][0],np.linalg.norm(rotVector[i]))
        assert(np.isclose(rotVector[i][0],np.linalg.norm(rotVector[i])))


def test_vectorAngle():
    v1 = np.array([1,0,0])
    v2 = np.array([0,1,0])
    v3 = np.array([1,1,1])
    theta1 = vectorAngle(v1,v2)
    theta2 = vectorAngle(v1,v3)
    theta3 = vectorAngle(v2,v3)
    print(theta2)
    assert(np.isclose(theta1,np.pi/2))
    assert(np.isclose(theta2,theta3))
    assert(np.isclose(theta2,0.955316618125))
    
def test_Rotation_matrix():
    M1 = rotationMatrix(0,0,0)

    assert(np.all(np.isclose(M1,np.identity(3))))

    M2 = rotationMatrix(90,0,0)
    M3 = rotationMatrix(np.pi/2,np.pi/3,np.pi/4,format='rad')
    M4 = rotationMatrix(90,60,45)
    
    assert(np.all(np.isclose(M2,np.array([[1,0,0],[0,0,-1],[0,1,0]]))))
    assert(np.all(np.isclose(M3,M4)))
    
    M4_check = np.array([[3.53553391e-01,6.12372436e-01,7.07106781e-01],[3.53553391e-01,6.12372436e-01,-7.07106781e-01],[-8.66025404e-01,5.00000000e-01,3.06161700e-17]])
    assert(np.all(np.isclose(M4,M4_check)))




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
def test_fileListGenerator():
    numberStr = '0,20-23-24,4000'
    year = 2018
    folder = os.path.join('','home','camea')
    try:
        files = fileListGenerator(numberStr,folder,year)
        assert False # Too many dasches
    except AttributeError:
        assert True

    numberStr = '0,20-23,4000'

    files = fileListGenerator(numberStr,folder,year)
    filesCorrect = np.array([
        os.path.join('','home','camea','camea2018n000000.hdf'), os.path.join('','home','camea','camea2018n000020.hdf'), 
        os.path.join('','home','camea','camea2018n000021.hdf'), os.path.join('','home','camea','camea2018n000022.hdf'), 
        os.path.join('','home','camea','camea2018n000023.hdf'), os.path.join('','home','camea','camea2018n004000.hdf')])
    assert(np.all(filesCorrect == np.array(files)))


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


def test_generateLabel():
    v = [1,-1,2]
    label = generateLabel(v)
    assert(label=='(H, -K, 2L)')

    v = np.array([-0.5,22.0000,-3.0])
    label = generateLabel(v)
    assert(label=='(-0.5H, 22K, -3L)')


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