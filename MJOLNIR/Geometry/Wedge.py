import math,numpy as np
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Wedge(GeometryConcept.GeometryConcept):
    """Wedge object to keep track of analysers and detectors. To be used as a storage object and facilitate easy movement of multiple detectors and analysers as once."""
    def __init__(self,position=(0,0,0),detectors=[],analysers=[]):
        """
        Args:

            position (float 3): Position of wedge (default (0,0,0))

        Kwargs:

            detectors (list or single detector): Either a list or a single detector (default empty)

            analysers (list or single analyser): Either a list or a single analyser (default empty)

        .. note::
            A wedge does not have a direction. The direction of analysers and detectors are to be set individually.

        """
        super(Wedge,self).__init__(position)
        self._analysers = []
        self._detectors = []

        self.append(analysers)
        self.append(detectors)
        

    @property
    def analysers(self):
        return self._analysers

    @analysers.getter
    def analysers(self):
        return self._analysers

    @analysers.setter
    def analysers(self,Analysers):
        if len(self.analysers)!=0:
            warnings.warn('The list of analysers is not empty! Appending new analyser(s)')
        if isinstance(Analysers, list):
            for ana in Analysers:
                if not issubclass(type(ana),Analyser.Analyser):
                    raise AttributeError('Object is not an analyser or a simple list of these')
                self._analysers.append(ana)
        else:
            if not issubclass(type(Analysers),Analyser.Analyser):
                raise AttributeError('Object is not an analyser or a simple list of these')
            self._analysers.append(Analysers)
        
    @property
    def detectors(self):
        return self._detectors

    @detectors.getter
    def detectors(self):
        return self._detectors

    @detectors.setter
    def detectors(self,Detectors):
        if len(self.detectors)!=0:
            warnings.warn('The list of detectors is not empty! Appending new detector(s)')
        if isinstance(Detectors, list):
            for det in Detectors:
                if not issubclass(type(det),Detector.Detector):
                    raise AttributeError('Object is not a detector or a simple list of these')
                self._detectors.append(det)
        else:
            if not issubclass(type(Detectors),Detector.Detector):
                raise AttributeError('Object is not a detector or a simple list of these')
            self._detectors.append(Detectors)

    def append(self,Object):
        """Append Object(s) to corresponding list.

        Args
            object (Detector(s)/Analyser(s)): Single detector/analyser of list of detectors/analysers
        """
        if isinstance(Object,list):
            for obj in Object:
                if issubclass(type(obj),Analyser.Analyser):
                    self._analysers.append(obj)
                elif issubclass(type(obj),Detector.Detector):
                    self._detectors.append(obj)
                else:
                    raise AttributeError('Object not analyser or detector or a simple list of these')
        else:
            if issubclass(type(Object),Analyser.Analyser):
                    self._analysers.append(Object)
            elif issubclass(type(Object),Detector.Detector):
                self._detectors.append(Object)
            else:
                raise AttributeError('Object not analyser or detector or a simple list of these')

    def plot(self,ax,offset=(0,0,0)):
        """Recursive plotting routine."""
        for obj in self.analysers+self.detectors:
            obj.plot(ax,offset=self.position+offset)

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