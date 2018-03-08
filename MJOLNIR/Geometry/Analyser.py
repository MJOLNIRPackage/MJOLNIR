import math,numpy as np
#import CAMEA
#import CAMEA.component
from MJOLNIR.Geometry import GeometryObject

class Analyser(GeometryObject.GeometryObject):
    """Generic analyser object"""
    def __init__(self,position,direction):
        super(Analyser,self).__init__(position,direction)
        self.type = "Generic Analyser"

    @property
    def type(self):
        return self._type

    @type.getter
    def type(self):
        return self._type

    @type.setter
    def type(self,type):
        self._type = type
        
    def plot(self,ax,offset=(0.0,0.0,0.0)):
        """
        Args:

            ax (matplotlib.pyplot 3d axis): Axis object into which the analyser is plotted

        Kwargs:

            offset (3vector): Offset of analuser due to bank position (default [0,0,0])
        
        >>> GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
        >>> GenericAnalyser.plot(ax)
        """
        raise NotImplementedError

    def __str__(self):
        returnString=('{} located at {}'.format(self.type,self.position))
        return returnString

def test_init():
    GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(np.all(GenericAnalyser.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericAnalyser.direction==(1.0,0.0,0.0)))
    assert(GenericAnalyser.type=="Generic Analyser")

def test_Generic_plot():
    
    GenericAnalyser = Analyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    ax = []
    try:
        GenericAnalyser.plot(ax)
        assert False
    except NotImplementedError:
        assert True




















class FlatAnalyser(Analyser):
    """Flat Analyser"""
    def __init__(self,position,direction,width=0.05,height=0.1):
        """
        Args:

            Position (3vector): Position of object (default [0,0,0])

            Direction (3vector): Direction along which the object points (default [0,0,1])

        Kwargs:

            length (float): Length of detector tube in meters (default 0.25)

            pixels (int): Number of pixels (default 456)

            diameter (float): Diameter of tube in meters (default 0.02)

        Raises:
            AttributeError, NotImplimentedError
        

        """
        super(FlatAnalyser,self).__init__(position,direction)
        self.type = "Flat Analyser"
        self.width = width
        self.height = height

    @property
    def width(self):
        return self._width

    @width.getter
    def width(self):
        return self._width

    @width.setter
    def width(self,width):
        if(width<0):
            raise AttributeError('The width of the analyser cannot be negative.')
        self._width = width 


    @property
    def height(self):
        return self._height

    @height.getter
    def height(self):
        return self._height

    @height.setter
    def height(self,height):
        if(height<0):
            raise AttributeError('The height of the analyser must be grater than 0.')
        self._height = height



    def plot(self,ax,bankPos=np.array([0,0,0]),n=100):
        """
        Args:

            ax (matplotlib.pyplot 3d axis): Axis object into which the analyser is plotted

        Kwargs:

            offset (3vector): Offset of detector due to bank position (default [0,0,0])

            n (int): Number of points on the surface to be plotted (default 100)
        
        >>> Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0))
        >>> Analyser.plot(ax,offset=(0.0,0.0,0.0),n=100)
        """
        pos = self.position.copy()

        pos+=bankPos
        Width = self.width / 2.0
        Height =self.height / 2.0

        angle = np.arctan2(self.position[1],self.position[0]) # angle relative to x-axis
        if self.direction.shape!=():
            angleZ = self.direction[0]
        else:
            angleZ = self.direction # Energy-determining angle
        
        x=np.linspace(-Width, Width, n)
        y=np.linspace(-Height/2, Height/2, n)
        Xc, Yc=np.meshgrid(x, y)
        Zcrot = np.sin(angleZ)*Xc

        Xcrot = Xc*np.cos(angle)-Yc*np.sin(angle)
        Ycrot = Xc*np.sin(angle)+Yc*np.cos(angle)

        ## Draw parameters
        rstride = 20
        cstride = 20

        ax.plot_surface(Xcrot+pos[0], Ycrot+pos[1], Zcrot+pos[2], alpha=1, rstride=rstride, cstride=cstride)


def test_Analyser_init():
    Analyser = FlatAnalyser(position=(0.0,1.0,0.0),direction=(1.0,0,0),width=0.05,height=0.1)
    assert(np.all(Analyser.position==np.array([0.0,1.0,0.0])))
    assert(np.all(Analyser.direction==(1.0,0.0,0.0)))
    assert(Analyser.type=="Flat Analyser")
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






def plane(point,width,height,normal):
    nppoint = np.array(point)
    npnormal = np.array(normal)


    npnormal = npnormal/np.linalg.norm(npnormal)
    theta = np.arccos(np.dot(np.array([1,0,0]),npnormal))-math.pi/2.0
    #print theta
    inplaneHorizontal = np.array([math.cos(theta),math.sin(theta),0])*width/2.0

    inplaneVertical = np.cross(npnormal,inplaneHorizontal)/(np.linalg.norm(inplaneHorizontal))*height/2.0



    points = np.array([nppoint+inplaneHorizontal+inplaneVertical,nppoint+inplaneHorizontal-inplaneVertical,
                nppoint-inplaneHorizontal-inplaneVertical,nppoint-inplaneHorizontal+inplaneVertical,
                nppoint+inplaneHorizontal+inplaneVertical])

    x,y,z= points.T
    return x,y,z