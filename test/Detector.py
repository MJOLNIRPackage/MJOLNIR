import numpy as np
from MJOLNIR.Geometry.Detector import Detector, TubeDetector1D
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def test_init():
    GenericDetector = Detector(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(np.all(GenericDetector.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericDetector.direction==(1.0,0.0,0.0)))

def test_Generic_plot():
    GenericDetector = Detector(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    try:
        GenericDetector.plot(ax)
        assert False
    except NotImplementedError:
        assert True



def test_TubeDetector_init():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0),pixels=20,length=0.3,diameter=0.025,split=[0,57,57*2])
    assert(np.all(TubeDetector.position==np.array([0.0,1.0,0.0])))
    assert(np.all(TubeDetector.direction==(1.0,0.0,0.0)))
    assert(TubeDetector.pixels==20)
    assert(TubeDetector.length==0.3)
    assert(TubeDetector.diameter==0.025)
    assert(np.all(TubeDetector.split==np.array([0,57,57*2])))

def test_TubeDetector_pixels():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        TubeDetector.pixels=0
        assert False
    except AttributeError:
        assert True

def test_TubeDetector_length():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        TubeDetector.length=-0.1
        assert False
    except AttributeError:
        assert True
    
def test_TubeDetector_diameter():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        TubeDetector.diameter=-0.1
        assert False
    except AttributeError:
        assert True

def test_TubeDetector_split():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0),pixels=100)
    try:
        TubeDetector.split=-0.1
        assert False
    except AttributeError:
        assert True

    TubeDetector.split=[50,60,100]
    pixelPos = TubeDetector.getPixelPositions()
    assert(len(pixelPos)==2)
    assert(len(pixelPos[0])==10)


def test_TubeDetector1D_plot():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    TubeDetector.plot(ax)
    

def test_TubeDetector1D_getPixelPositions():
    TubeDetector = TubeDetector1D(position=(1.0,0.0,1.0),direction=(1.0,0,0),length=0.5,pixels=5)
    positions = TubeDetector.getPixelPositions()

    AssumedPositions = np.array([[0.8,0,1],[0.9,0,1],[1.0,0,1],[1.1,0,1],[1.2,0,1]])
    print(positions)
    print(AssumedPositions)
    assert(np.all(AssumedPositions==positions))

    
    