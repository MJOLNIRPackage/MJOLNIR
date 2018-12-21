import sys
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import math,numpy as np
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector
from MJOLNIR import _tools
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Wedge(GeometryConcept.GeometryConcept):
    """Wedge object to keep track of analysers and detectors. To be used as a storage object and facilitate easy movement of multiple detectors and analysers as once."""
    @_tools.KwargChecker(include=[''])
    def __init__(self,position=(0.0,0.0,0.0),detectors=[],analysers=[],concept='ManyToMany',**kwargs):
        """
        Args:

            - position (float 3): Position of wedge (default (0,0,0))

        Kwargs:

            - detectors (list or single detector): Either a list or a single detector (default empty)

            - analysers (list or single analyser): Either a list or a single analyser (default empty)

            - concept (string "ManyToMany" or "OneToOne"): Setting to controle if there is a "one to one" correspondence between analysers and detectors or a "many to many" relationship.

        .. note::
            A wedge does not have a direction. The direction of analysers and detectors are to be set individually.

        """
        super(Wedge,self).__init__(position)
        self._analysers = []
        self._detectors = []

        self.append(analysers)
        self.append(detectors)
        self._settings = {}
        
        for key in kwargs:
            self.settings[key]=kwargs[key]
        self.settings['concept']=concept
        
        

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

    @property
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')

    def append(self,Object):
        """Append Object(s) to corresponding list.

        Args:

            - object (Detector(s)/Analyser(s)): Single detector/analyser of list of detectors/analysers
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

    @_tools.KwargChecker()
    def plot(self,ax,offset=(0,0,0)):
        """Recursive plotting routine."""
        for obj in self.analysers+self.detectors:
            obj.plot(ax,offset=np.array(self.position,dtype=float)+np.array(offset,dtype=float))

    def __str__(self):
        string = ''
        string+="{} with settings:\nPosition = {}\n".format(self.__class__,self.position)
        string+='Containing analysers:\n'
        for ana in self.analysers:
            string+='\t'+str(ana)+'\n'
        string+='Containing detectors:\n'
        for det in self.detectors:
            string+='\t'+str(det)+'\n'

        return string

    def calculateDetectorAnalyserPositions(self):
        """Find neutron position on analyser and detector. Assuming that the analyser is in the z=0 plane."""
        if(len(self.detectors)==0 or len(self.analysers)==0):
            raise ValueError('Wedge does not contain detectors and/or analysers.')

        detectorPixelPositions = []
        analyserPixelPositions = []
        if self.settings['concept']=='OneToOne':
            if len(self.detectors)!=len(self.analysers):
                raise RuntimeError('Concept set to OneToOne but number of detectors does not mach analysers ({}!?{}'.format(len(self.detectors),len(self.analysers)))
            detectorCounter = 0
            for det in self.detectors:
                PixelPos = det.getPixelPositions()+self.position
                if len(PixelPos)!=1:
                    raise ValueError("OneToOne concept chosen by detector split into multiple parts!")
                detectorPixelPositions.append(PixelPos)
                LDA = PixelPos[0]-self.analysers[detectorCounter].position # Detector - analyser vector
                LAS = self.analysers[detectorCounter].position+self.position
                vertical = np.array([0,0,1])

                perpVect = np.cross(vertical,LAS)

                deltaXD = np.dot(PixelPos,perpVect)

                LD = np.linalg.norm(LDA,axis=1)
                LA = np.linalg.norm(LAS)

                deltaXDprime = deltaXD/(LD/LA+1.0)

                analyserPixelPositions.append(np.outer(deltaXDprime,perpVect)+self.analysers[detectorCounter].position+self.position)
                detectorCounter+=1

        elif self.settings['concept']=='ManyToMany':
            for det in self.detectors:
                PixelPos = [x +self.position for x in det.getPixelPositions()]
                if len(PixelPos)!=len(self.analysers):
                    raise ValueError("ManyToMany concept chosen by detector split into number of parts not matching number of analysers!")
                detectorPixelPositions.append(np.concatenate(PixelPos))
                
                vertical = np.array([0,0,1])
                
                LAS = [self.analysers[i].position+self.position for i in range(len(PixelPos))]
                perpVect = [np.cross(vertical,LAS[i])/(np.linalg.norm(np.cross(vertical,LAS[i]))) for i in range(len(PixelPos))]
                
                deltaXD = [np.dot(PixelPos[i],perpVect[i]) for i in range(len(PixelPos))]
                LDA = [np.array(PixelPos[i])-np.array(deltaXD[i]).reshape(-1,1)*np.array(perpVect[i]).reshape(1,3)-self.analysers[i].position for i in range(len(PixelPos))]
                
                LD = [np.linalg.norm(LDA[i],axis=1) for i in range(len(PixelPos))]
                LA = [np.linalg.norm(LAS[i]) for i in range(len(PixelPos))]
                deltaXDprime = [deltaXD[i]/(LD[i]/LA[i]+1.0) for i in range(len(PixelPos))]


                analyserPixelPositions.append(np.concatenate([np.outer(deltaXDprime[i],perpVect[i])+self.analysers[i].position+self.position for i in range(len(PixelPos))]))

        else:
            raise ValueError("Wedge does not contain a Concept setting that is understood. Should be either 'OneToOne' or 'ManyToMany'")

        return detectorPixelPositions,analyserPixelPositions

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
    