from MJOLNIR.Geometry.Detector import Detector
from MJOLNIR.Geometry.Wedge import Wedge
import MJOLNIR.Geometry.Detector as Detector
import MJOLNIR.Geometry.Analyser as Analyser
import warnings
import numpy as np
import matplotlib.pyplot as plt



def test_Wedge_init():
    wedge = Wedge(position=(0,0,0))

    Det = Detector.Detector(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.Analyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.detectors=Det
    wedge.analysers=Ana

    wedge2 = Wedge(position=(0,0,0))
    wedge2.detectors=[Det,Det]
    wedge2.analysers=[Ana,Ana]

    wedge3 = Wedge(detectors=[Det,Det],analysers=Ana)
    assert(len(wedge3.analysers)==1)
    assert(len(wedge3.detectors)==2)


def test_Wedge_error():
    wedge = Wedge(position=(0,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    try:
        wedge.detectors=Ana
        assert False
    except AttributeError:
        assert True

    try:
        wedge.analysers=Det
        assert False
    except AttributeError:
        assert True

    try:
        wedge.append("Wrong object type")
        assert False
    except AttributeError:
        assert True
    
    try:
        wedge.append(["List of",3.0,"wrong objects"])
        assert False
    except AttributeError:
        assert True


    try:
        wedge.settings['concept']='OneToOne'
        wedge.calculateDetectorAnalyserPositions()
        assert False
    except ValueError:
        assert True

    try:
        wedge.settings['concept']='OneToOne'
        wedge.append([Det,Det,Ana])
        wedge.calculateDetectorAnalyserPositions()
        assert False
    except RuntimeError:
        assert True

    try:
        wedge.settings['concept']='ManyToMany'
        wedge.detectors[0].split=[10,30,44]
        wedge.calculateDetectorAnalyserPositions()
        assert False
    except ValueError:
        assert True

    try:
        wedge.settings['concept']='Wrong'
        wedge.calculateDetectorAnalyserPositions()
        assert False
    except ValueError:
        assert True


def test_Wedge_warnings():
    wedge = Wedge()

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.detectors = Det
    wedge.analysers = Ana
    with warnings.catch_warnings(record=True) as w: # From https://docs.python.org/3.1/library/warnings.html
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        wedge.detectors = Det
        wedge.analysers = Ana
        # Verify some things
        assert len(w) == 2
        assert issubclass(w[0].category, UserWarning)
        assert issubclass(w[1].category, UserWarning)
        assert 'The list of detectors is not empty! Appending new detector(s)' in str(w[0].message)
        assert 'The list of analysers is not empty! Appending new analyser(s)' in str(w[1].message)

def test_Wedge_append():
    wedge = Wedge()

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Det,Ana])
    wedge.append(Det)
    wedge.append(Ana)

    assert(len(wedge.detectors)==3)
    assert(len(wedge.analysers)==2)

def test_Wedge_plot():
    wedge = Wedge()
    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    wedge.plot(ax)

def test_wedge_calculateDetectorAnalyserPositions_OneToOne():
    wedge = Wedge(concept='OneToOne')
    Det = Detector.TubeDetector1D(position=(1.0,0.0,1.0),direction=(1.0,0,0),length=0.5,pixels=5)
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    Det2 = Detector.TubeDetector1D(position=(1.5,0.1,1.0),direction=(1.0,0,0),length=0.5,pixels=5)
    Ana2 = Analyser.FlatAnalyser(position=(0.75,0,0),direction=(1,0,1))

    wedge.append([Det,Det2,Ana,Ana2])

    detectorPixelPositions,analyserPixelPositions = wedge.calculateDetectorAnalyserPositions()
    #print(detectorPixelPositions,analyserPixelPositions)

    DetPixelPos = np.array([[[0.8,0,1],[0.9,0,1],[1.0,0,1],[1.1,0,1],[1.2,0,1]]])
    assert(np.all(DetPixelPos==detectorPixelPositions[0][0]))
    assert(np.all([x==Ana.position for x in analyserPixelPositions[0]]))
    
    DetPixelPos2= np.array([[[1.3,0.1,1],[1.4,0.1,1],[1.5,0.1,1],[1.6,0.1,1],[1.7,0.1,1]]])
    #print(analyserPixelPositions[1])
    assert(np.all(DetPixelPos2==detectorPixelPositions[1][0]))
    assert(np.all([x[0]==Ana2.position[0] for x in analyserPixelPositions[1]]))

    offcenterpos = np.array([0.02225497,0.02166939,0.02105171,0.02041748,0.01977911])
    #print(np.sum([analyserPixelPositions[1][i][1]-offcenterpos[i] for i in range(5)]))
    assert(np.sum([analyserPixelPositions[1][i][1]-offcenterpos[i] for i in range(5)])<1e-8)
    
    

def test_wedge_calculateDetectorAnalyserPositions_ManyToMany():
    
    wedge = Wedge(concept='ManyToMany')

    
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    Det2 = Detector.TubeDetector1D(position=(1.5,0.1,1.0),direction=(1.0,0,0),length=0.5,pixels=5,split=[0,2,5])
    Ana2 = Analyser.FlatAnalyser(position=(0.75,0,0),direction=(1,0,1))

    wedge.append([Det2,Ana,Ana2])

    detectorPixelPositions,analyserPixelPositions = wedge.calculateDetectorAnalyserPositions()

    #print(detectorPixelPositions)
    #print(analyserPixelPositions)
    
    DetPixelPos2= np.array([[[1.3,0.1,1],[1.4,0.1,1],[1.5,0.1,1],[1.6,0.1,1],[1.7,0.1,1]]])
    assert(np.all(DetPixelPos2==detectorPixelPositions[0]))
    
    anapos = np.array([0.5,0.5,0.75,0.75,0.75])
    assert(np.all([analyserPixelPositions[0][i,0]==anapos[i] for i in range(5)]))

    #offcenterpos = np.array([0.00700467,0.00676014,0.02105171,0.02041748,0.01977911])
    #assert(np.sum([analyserPixelPositions[0][i][1]-offcenterpos[i] for i in range(5)])<1e-8)

def test_wedge_string_dummy():
    wedge = Wedge(concept='ManyToMany')

    string = str(wedge)
    assert True

def test_wedge_repr_dummy():
    wedge = Wedge(concept='ManyToMany')

    string = repr(wedge)
    assert True
    
