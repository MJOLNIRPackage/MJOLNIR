from MJOLNIR.Geometry.Analyser import Analyser, FlatAnalyser
import numpy as np
import warnings
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def test_init():
    GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0),d_spacing=3.35,mosaicity=60.0)
    assert(np.all(GenericAnalyser.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericAnalyser.direction==(1.0,0.0,0.0)))
    assert(GenericAnalyser.d_spacing==3.35)
    assert(GenericAnalyser.mosaicity==60.0)

def test_Generic_plot():
    
    GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0),d_spacing=3.35,mosaicity=60.0)
    ax = []
    try:
        GenericAnalyser.plot(ax)
        assert False
    except NotImplementedError:
        assert True

def test_warnings():
    GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0),d_spacing=3.35,mosaicity=60.0)

    with warnings.catch_warnings(record=True) as w: # From https://docs.python.org/3.1/library/warnings.html
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        GenericAnalyser.d_spacing = 20.0
        GenericAnalyser.mosaicity = 0.5
        # Verify some things
        assert len(w) == 2
        assert issubclass(w[0].category, UserWarning)
        assert issubclass(w[1].category, UserWarning)
        assert "The unit of d spacing is Angstrom" in str(w[0].message)
        assert "The unit of mosaicity is arcminutes" in str(w[1].message)

def test_Generic_errors():
    GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0),d_spacing=3.35,mosaicity=60.0)
    try:
        GenericAnalyser.d_spacing=-0.1
        assert False
    except AttributeError:
        assert True

    try:
        GenericAnalyser.mosaicity=-0.1
        assert False
    except AttributeError:
        assert True


def test_Analyser_init():
    Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0),d_spacing=3.35,mosaicity=60,width=0.05,height=0.1)
    assert(np.all(Analyser.position==np.array([0.0,1.0,0.0])))
    assert(np.all(Analyser.direction==(1.0,0.0,0.0)))
    assert(Analyser.width==0.05)
    assert(Analyser.height==0.1)
 

def test_Analyser_width():
    Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        Analyser.width=-0.1
        assert False
    except AttributeError:
        assert True

def test_Analyser_height():
    Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        Analyser.height=-0.1
        assert False
    except AttributeError:
        assert True


def test_FlatAnalyser_plot():
    Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    Analyser.plot(ax)
    
def test_FlatAnalyser_str():
    Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(str(Analyser)=='FlatAnalyser located at '+np.array2string(np.array([0.0,1.0,0.0])))


def test_FlatAnalyser_kwChecker():

    try:
        Analyser = FlatAnalyser(Position=(0.0,1.0,0.0),direction=(1.0,0,0))
        assert False
    except:
        assert True

    try:
        Analyser = FlatAnalyser(pos=(0.0,1.0,0.0),direction=(1.0,0,0))
        assert False
    except:
        assert True
