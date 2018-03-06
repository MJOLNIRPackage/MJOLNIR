"""
    GeometryObject
    ==============
    
    General geometry object on which other MJOLNIR components are build.
    Generates an empty object with attributes:

    * Position (3d vector)
    * Direction (3d vector)

    Further, the object has the following methods:

    * setPosition

        """


import math
import numpy as np


class GeometryObject(object):
    """
    GeometryObject
    ==============
    
    General geometry object on which other MJOLNIR components are build.
    Generates an empty object with attributes:

        * Position (3d vector)
        * Direction (3d vector)

    Further, the object has the following methods:
        * setPosition
        """

    

    def __init__(self, position=(0.0,0.0,0.0), direction=(0,0,1)):
        """
        Init
        ----
        
        Initialization of generic GeometryObject.
        Input is:
            * Position (3d vector)
            * Direction (3d vector)
        """

        self.set_position(position)
        self.set_direction(direction)
    
    @property
    def _position(self):
        return self._position

    def get_position(self):
        return self._position

    def set_position(self,position):
        position = np.atleast_3d(position)
        self._position = position

    @property
    def _direction(self):
        return self._direction


    def get_direction(self):
        return self._direction

    def set_direction(self,direction):
        direction = np.atleast_3d(direction)
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