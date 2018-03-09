import math,numpy as np
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector,Wedge
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Instrument(GeometryConcept.GeometryConcept):
    def __init__(self, position=(0,0,0),wedges=[]):
        """Instrument object used to calculated analytic scattering coverage. Based on the GeometryConcept object it contains all needed information
        about the setup used in further calculations.

        Kwargs:
            
            position (float 3d): Position of the instrument alwasy at origin(?) (default (0,0,0))

            wedges (list of wedges or single wedge): Wedge or list of wedges which the instrument consists of (default empty)

        Raises:
            AttributeError
        
        """
        super(Instrument,self).__init__(position)
        self._wedges = []

        self.append(wedges)
        self._settings = {}
        self._settings['Initialized']=False
        

    @property
    def wedges(self):
        return self._wedges

    @wedges.getter
    def wedges(self):
        return self._wedges

    @wedges.setter
    def wedges(self,wedges):
        if len(self.wedges)!=0:
            warnings.warn('The list of wedges is not empty! Appending new wedges(s)')
        if isinstance(wedges, list):
            for ana in wedges:
                if not issubclass(type(ana),Wedge.Wedge):
                    raise AttributeError('Object is not an wedge or a simple list of these')
                self._wedges.append(ana)
        else:
            if not issubclass(type(wedges),Wedge.Wedge):
                raise AttributeError('Object is not an analyser or a simple list of these')
            self._wedges.append(wedges)
    
    def append(self,wedge):
        """Append wedge(s) to instrument.

        Args
            wedge (Wedge(s)): Single wedge or list of wedges
        """
        if isinstance(wedge,list):
            for obj in wedge:
                if issubclass(type(obj),Wedge.Wedge):
                    self._wedges.append(wedge)
                else:
                    raise AttributeError('Object not wedge or a simple list of wedges')
        else:
            if issubclass(type(wedge),Wedge.Wedge):
                    self._wedges.append(wedge)
            else:
                raise AttributeError('Object not wedge or a simple list of wedges')

    def plot(self,ax):
        """Recursive plotting routine."""
        for wedge in self.wedges:
            wedge.plot(ax,offset=self.position)

    @property
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')

            


def test_Instrument_init():
    Instr = Instrument()

    assert(np.all(Instr.position==(0,0,0)))

    Det = Detector.Detector(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.Analyser(position=(0.5,0,0),direction=(1,0,1))
    
    wedge = Wedge.Wedge(detectors=[Det,Det],analysers=Ana)

    Instr.wedges=[wedge,wedge]

    assert(Instr.settings['Initialized']==False)



def test_Instrument_error():
    Instr = Instrument()

    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    try:
        Instr.wedges=Ana
        assert False
    except AttributeError:
        assert True

    try:
        Instr.append("Wrong object type")
        assert False
    except AttributeError:
        assert True
    
    try:
        Instr.append(["List of",3.0,"wrong objects"])
        assert False
    except AttributeError:
        assert True

    try:
        Instr.settings = {'Name','New dictionary'}
        assert False
    except NotImplementedError:
        return True


def test_Instrument_warnings():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.wedges = wedge

    with warnings.catch_warnings(record=True) as w: # From https://docs.python.org/3.1/library/warnings.html
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        Instr.wedges = wedge
        # Verify some things
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert 'The list of wedges is not empty! Appending new wedges(s)' in str(w[0].message)


def test_Instrument_append():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.append([wedge,wedge])
    Instr.append(wedge)

    assert(len(Instr.wedges)==3)


def test_Instrument_plot():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    Instr.plot(ax)

def test_Instrument_Setting(): 
    Instr = Instrument()
    #Instr.updateSetting('New',True)
    Instr.settings['SettingVersion']=1.0
    assert(Instr.settings['SettingVersion']==1.0)
