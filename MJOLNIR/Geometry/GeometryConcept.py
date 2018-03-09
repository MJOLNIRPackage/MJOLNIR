import math
import numpy as np

class GeometryConcept(object):
    """Abstract geometry concept. Used as base class for Wedge and Instrument."""
    def __init__(self,position=(0,0,0)):
        self.position = position
        """
        Kwargs:

            Position (3vector): Position of object (default [0,0,0])

        Raises:
            AttributeError

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
        
        position  = np.array(position)
        if position.ndim !=1 or len(position)!=3:
            raise AttributeError('Position is to be a 3 vector')
        self._position = position

def test_Concept_init():
    GenericObject = GeometryObject(position=(0.0,1.0,0.0))
    assert(np.all(GenericObject.position==np.array([0.0,1.0,0.0])))




class GeometryObject(GeometryConcept):
    """Physical geometry object on which other physical MJOLNIR components are build. All of the components needed to create an instrument should
    inherit from this class in order enforce a uniform interface."""
    

    def __init__(self, position=(0.0,0.0,0.0), direction=(0,0,1)):
        """
        Kwargs:

            Position (3vector): Position of object (default [0,0,0])

            Direction (3vector): Direction along which the object points (default [0,0,1])

        Raises:
            AttributeError

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

    # def setDirection(self,direction):
    #     self.mDirection = np.array(direction)
    
    # def setDirectionRadians(self,angles):
    #     xdirection = math.sin(angles[1])*math.cos(angles[0])
    #     ydirection = math.sin(angles[1])*math.sin(angles[0])
    #     zdirection = math.cos(angles[1])
    #     self.mDirection = (xdirection,ydirection,zdirection)

    # def setDirectionDegrees(self,anglesd):
    #     anglesConverted = (anglesd[1]*math.pi/180.0,anglesd[0]*math.pi/180.0,anglesd[1]*math.pi/180.0)
    #     self.setDirectionRadians(anglesConverted)

    # def getDirectionRadians(self):
    #     theta = math.atan2(self.mDirection[1],self.mDirection[0])
    #     phi = math.acos(self.mDirection[2])
    #     return [theta,phi]

    # def getDirectionDegrees(self):
    #     angles = getDirectionRadians()
    #     return [math.degrees(angles[0]),math.degrees(angles[1])]


    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.__dict__)

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