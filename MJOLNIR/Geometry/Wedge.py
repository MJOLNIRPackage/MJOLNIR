import math,numpy as np
from MJOLNIR.Geometry import GeometryObject,Analyser,Detector
import warnings

class Wedge(GeometryObject.GeometryObject):
    """Wedge object"""
    def __init__(self,position,direction):
        """
        Args:

            position (float 3): Position of wedge

            direction (float 3): Direction of wedge

        """
        super(Wedge,self).__init__(position,direction)
        self._Analysers = []
        self._Detectors = []

    @property
    def Analysers(self):
        return self._Analysers

    @Analysers.getter
    def Analysers(self):
        return self._Analysers

    @Analysers.setter
    def Analysers(self,Analysers):
        if len(self.Analysers)!=0:
            warnings.warn('The list of analysers is not empty! Appending new analyser(s)')
        if isinstance(Analysers, list):
            for ana in Analysers:
                if not issubclass(type(ana),Analyser.Analyser):
                    raise AttributeError('Object is not an analyser or a simple list of these')
                self._Analysers.append(ana)
        else:
            if not issubclass(type(Analysers),Analyser.Analyser):
                raise AttributeError('Object is not an analyser or a simple list of these')
            self._Analysers.append(Analysers)
        
    @property
    def Detectors(self):
        return self._Analysers

    @Detectors.getter
    def Detectors(self):
        return self._Detectors

    @Detectors.setter
    def Detectors(self,Detectors):
        if len(self.Detectors)!=0:
            warnings.warn('The list of detectors is not empty! Appending new detector(s)')
        if isinstance(Detectors, list):
            for det in Detectors:
                if not issubclass(type(det),Detector.Detector):
                    raise AttributeError('Object is not a detector or a simple list of these')
                self._Detectors.append(det)
        else:
            if not issubclass(type(Detectors),Detector.Detector):
                raise AttributeError('Object is not a detector or a simple list of these')
            self._Detectors.append(Detectors)

    def append(self,Object):
        """Append Object(s) to corresponding list.

        Args
            object (Detector(s)/Analyser(s)): Single detector/analyser of list of detectors/analysers
        """
        if isinstance(Object,list):
            for obj in Object:
                if issubclass(type(obj),Analyser.Analyser):
                    self._Analysers.append(obj)
                elif issubclass(type(obj),Detector.Detector):
                    self._Detectors.append(obj)
                else:
                    raise AttributeError('Object not analyser or detector or a simple list of these')
        else:
            if issubclass(type(Object),Analyser.Analyser):
                    self._Analysers.append(Object)
            elif issubclass(type(Object),Detector.Detector):
                self._Detectors.append(Object)
            else:
                raise AttributeError('Object not analyser or detector or a simple list of these')

def test_Wedge_init():
    wedge = Wedge(position=(0,0,0),direction=(0,0,1))

    Det = Detector.Detector(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.Analyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.Detectors=Det
    wedge.Analysers=Ana

def test_Wedge_error():
    wedge = Wedge(position=(0,0,0),direction=(0,0,1))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    try:
        wedge.Detectors=Ana
        assert False
    except AttributeError:
        assert True

    try:
        wedge.Analysers=Det
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
    wedge = Wedge(position=(0,0,0),direction=(0,0,1))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.Detectors = Det
    wedge.Analysers = Ana
    with warnings.catch_warnings(record=True) as w: # From https://docs.python.org/3.1/library/warnings.html
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        wedge.Detectors = Det
        wedge.Analysers = Ana
        # Verify some things
        assert len(w) == 2
        assert issubclass(w[0].category, UserWarning)
        assert issubclass(w[1].category, UserWarning)
        assert 'The list of detectors is not empty! Appending new detector(s)' in str(w[0].message)
        assert 'The list of analysers is not empty! Appending new analyser(s)' in str(w[1].message)

def test_Wedge_append():
    wedge = Wedge(position=(0,0,0),direction=(0,0,1))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    wedge.append(Det)
    wedge.append(Ana)

    assert(len(wedge.Detectors)==2)
    assert(len(wedge.Analysers)==2)

