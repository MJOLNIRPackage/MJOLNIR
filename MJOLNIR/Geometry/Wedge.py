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
