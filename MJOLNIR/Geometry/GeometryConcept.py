import sys
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle as pickle
from MJOLNIR import _tools

class GeometryConcept(object):
    """Abstract geometry concept. Used as base class for Wedge and Instrument."""

    @_tools.KwargChecker()
    def __init__(self,position=(0,0,0)):
        self.position = position
        """
        Kwargs:

            - Position (3vector): Position of object (default [0,0,0])

        Raises:

            - AttributeError

            - NotImplementedError

        >>> GenericConcept = GeometryConcept(position=(0.0,1.0,0.0))
        >>> print(GenericConcept.position)
        (0.0,1.0,0.0)
        """

    @property
    def position(self):
        return self._position

    @position.getter
    def position(self):
        return self._position

    @position.setter
    def position(self,position):
        
        position  = np.array(position,dtype=float)
        if position.ndim !=1 or len(position)!=3:
            raise AttributeError('Position is to be a 3 vector, got {}.'.format(position))
        self._position = position

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

    def __str__(self):
        returnString=('{} located at {}'.format(str(self.__class__).split('.')[-1][:-2],self.position))
        return returnString

    def plot(self,ax):
        """
        Args:

            - ax (matplotlib axis): 3D matplotlib axis into whicht plotting is performed

        .. warning::
            Method not incorporated, but acts as virtual method.

        """
        raise NotImplementedError

    def __eq__(self, other): 
        return np.logical_and(set(self.__dict__.keys()) == set(other.__dict__.keys()),self.__class__ == other.__class__)

    def save(self, filename):
        try:                                # Opening the given file with an error catch
            fileObject = open(filename, 'wb')
        except IOError as e:   # pragma: no cover
            print("Error in opening file:\n{}".format(e))
        else:
                pickle.dump(self, fileObject, -1)
                fileObject.close()  

    def load(self,filename):
        """Method to load an object from a pickled file."""
        try:                                # Opening the given file with an error catch
            fileObject = open(filename, 'rb')
        except IOError as e: # pragma: no cover
            print("Error in opening file:\n{}".format(e))
        else:
            tmp_dict = pickle.load(fileObject)
            
            fileObject.close()
            # TODO: Make checks that the object loaded is of correct format?
            self=tmp_dict
        

def test_Concept_init():
    Concept = GeometryConcept(position=(0.0,1.0,0.0))
    assert(np.all(Concept.position==np.array([0.0,1.0,0.0])))


def test_Concept_plot():
    Concept = GeometryConcept()
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    try:
        Concept.plot(ax)
        assert False
    except NotImplementedError:
        assert True
    print(str(Concept))
    



class GeometryObject(GeometryConcept):
    """Physical geometry object on which other physical MJOLNIR components are build. All of the components needed to create an instrument should
    inherit from this class in order enforce a uniform interface."""
    
    @_tools.KwargChecker()
    def __init__(self, position=(0.0,0.0,0.0), direction=(0,0,1)):
        """
        Kwargs:

            - Position (3vector): Position of object (default [0,0,0])

            - Direction (3vector): Direction along which the object points (default [0,0,1])

        Raises:
            
            - AttributeError

        >>> GenericObject = GeometryObject(position=(0.0,1.0,0.0),direction=(1.0,0,0))
        >>> print(GenericObject.position)
        (0.0,1.0,0.0)
        """
        super(GeometryObject,self).__init__(position)
        self.direction = direction
    


    @property
    def direction(self):
        return self._direction

    @direction.getter
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self,direction):
        direction  = np.array(direction,dtype=float)
        if direction.ndim !=1 or len(direction)!=3:
            raise AttributeError('Direction is to be a 3 vector')
        if np.abs(np.linalg.norm(direction))<1e-10:
            raise AttributeError('Length of direction is not allowed to be zero')
        
        direction/=np.linalg.norm(direction)
        self._direction = direction


    def __str__(self):
        return "Position = {}\tDirection = {}".format(self._position,self._direction)



# Test of GeometryObject
def test_Object_init():
    GenericObject = GeometryObject(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(np.all(GenericObject.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericObject.direction==(1.0,0.0,0.0)))

def test_position():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    GenericObject.position = (0.0,0.0,0.0)
    print(str(GenericObject))
    assert(np.all(GenericObject.position==(0.0,0.0,0.0)))

def test_Object_position_exception():
    GenericConcept = GeometryConcept(position=(0,1.0,0.0))
    try:
        GenericConcept.position=((0,0),(0,0))
        assert False
    except AttributeError:
        assert True

    try:
        GenericConcept.position=(0,0,0,0)
        assert False
    except AttributeError:
        assert True

def test_Object_direction():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    GenericObject.direction = (0.0,0.0,0.5)
    assert(np.all(GenericObject.direction==(0.0,0.0,1.0)))

def test_Object_direction_exception():
    GenericObject = GeometryObject(position=(0,1.0,0.0),direction=(1.0,0,0))
    try:
        GenericObject.direction=((0,0),(0,0))
        assert False
    except AttributeError:
        assert True

    try:
        GenericObject.direction=(0,0,0)
        assert False
    except AttributeError:
        assert True

    try:
        GenericObject.direction=(0,0,0,1)
        assert False
    except AttributeError:
        assert True