from MJOLNIR.Geometry import GeometryObject
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Detector(GeometryObject.GeometryObject):
    """Generic detector being the base class of all detectors."""

    def __init__(self, position,direction):
        """
        Kwargs:
            Position (3vector): Position of object (default [0,0,0])

            Direction (3vector): Direction along which the object points (default [0,0,1])"""
        super(Detector,self).__init__(position,direction)
        self.type = "Generic Detector"

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

            ax (matplotlib.pyplot 3d axis): Axis object into which the detector is plotted

        Kwargs:

            offset (3vector): Offset of detector due to bank position (default [0,0,0])
            
        >>> GenericDetector = Detector(position=(0.0,1.0,0.0),direction=(1.0,0,0))
        >>> GenericDetector.plot(ax)
        Plots detector tube in provided axis object.

        """
        raise NotImplementedError

    def __str__(self):
        returnString=('{} located at {}'.format(self.type,self.position))
        return returnString

def test_init():
    GenericDetector = Detector(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    assert(np.all(GenericDetector.position==np.array([0.0,1.0,0.0])))
    assert(np.all(GenericDetector.direction==(1.0,0.0,0.0)))
    assert(GenericDetector.type=="Generic Detector")

def test_Generic_plot():
    GenericDetector = Detector(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    ax = []
    try:
        GenericDetector.plot(ax)
        assert False
    except NotImplementedError:
        assert True

class TubeDetector1D(Detector):
    """1D Tube detector used at PSI. The detector is assumed to be a perfect cylinder consisting of pixels."""
    def __init__(self, position, direction,length=0.25, pixels=456,diameter=0.02):
        """
        Args:

            Position (3vector): Position of object (default [0,0,0])

            Direction (3vector): Direction along which the object points (default [0,0,1])

        Kwargs:

            length (float): Length of detector tube in meters (default 0.25)

            pixels (int): Number of pixels (default 456)

            diameter (float): Diameter of tube in meters (default 0.02)

        Raises:
            AttributeError
        

        """

        super(TubeDetector1D,self).__init__(position,direction)
        self.type = "1D Tube Detector"
        self.pixels = pixels
        self.length = length
        self.diameter = diameter

    @property
    def pixels(self):
        return self._pixels

    @pixels.getter
    def pixels(self):
        return self._pixels

    @pixels.setter
    def pixels(self,pixels):
        if(pixels<1):
            raise AttributeError('The number of pixels need to be greater than 0.')
        self._pixels = int(pixels)   

    @property
    def length(self):
        return self._length

    @length.getter
    def length(self):
        return self._length

    @length.setter
    def length(self,length):
        if(length<0):
            raise AttributeError('The lenght of the detector tube must be grater than 0.')
        self._length = length

    @property
    def diameter(self):
        return self._diameter

    @diameter.getter
    def diameter(self):
        return self._diameter

    @diameter.setter
    def diameter(self,diameter):
        if(diameter<0):
            raise AttributeError('The diameter of the detector tube must be grater than 0')
        self._diameter = diameter



    def plot(self,ax,offset=(0.0,0.0,0.0),n=100):
        """
        Args:

            ax (matplotlib.pyplot 3d axis): Axis object into which the detector is plotted

        Kwargs:

            offset (3vector): Offset of detector due to bank position (default [0,0,0])

            n (int): Number of points on the surface to be plotted (default 100)
        
        >>> Detector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
        >>> Detector.plot(ax,offset=(0.0,0.0,0.0),n=100)
        Plots detector tube in provided axis object.

        """
        
        pos = self.position.copy()
        pos+=offset


        r = self.diameter/2

        x=np.linspace(-r, r, n)
        z=np.linspace(-self.length/2, self.length/2, n)
        Xc, Zc=np.meshgrid(x, z)
        Yc = np.sqrt(r**2-Xc**2)


        a = np.array([0.0,0.0,1.0]) # Original direction
        b = self.direction.copy()    # Correct direction

        b.shape = (1,3)


        # Rotation matrix found from http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

        v = np.cross(a.T,b)
        s = np.linalg.norm(v)
        vmat=np.array([0.0,-v[0,2],v[0,1], v[0,2],0.0,-v[0,0],  -v[0,1],v[0,0],0.0])
        vmat.shape = (3,3)
        c = np.dot(a.T,b[0,:])
        R = np.eye(3)+vmat+(1.0-c)/(pow(s,2.0))*np.dot(vmat,vmat)


        res = np.zeros((n,n,3))
        res2 = np.zeros((n,n,3))

        for i in range(n):
            for j in range(n):
                res[i,j,:]=np.dot(np.array([Xc[i,j],Yc[i,j],Zc[i,j]]), R)
                res2[i,j,:]=np.dot(np.array([Xc[i,j],-Yc[i,j],Zc[i,j]]), R)

        Xcrot=res[:,:,0]
        Ycrot=res[:,:,1]
        Zcrot=res[:,:,2]

        Xcrot2=res2[:,:,0]
        Ycrot2=res2[:,:,1]
        Zcrot2=res2[:,:,2]

        # Draw parameters
        rstride = np.round(n/5).astype(int)
        cstride = np.round(3*n/10).astype(int)

        ax.plot_surface(Xcrot+pos[0], Ycrot+pos[1], Zcrot+pos[2], alpha=1.0, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xcrot2+pos[0], Ycrot2+pos[1], Zcrot2+pos[2], alpha=1.0, rstride=rstride, cstride=cstride)





def test_TubeDetector_init():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0),pixels=20,length=0.3,diameter=0.025)
    assert(np.all(TubeDetector.position==np.array([0.0,1.0,0.0])))
    assert(np.all(TubeDetector.direction==(1.0,0.0,0.0)))
    assert(TubeDetector.type=="1D Tube Detector")
    assert(TubeDetector.pixels==20)
    assert(TubeDetector.length==0.3)
    assert(TubeDetector.diameter==0.025)

def test_TubeDetector_pixels():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        TubeDetector.pixels=0
        assert False
    except AttributeError:
        assert True

def test_TubeDetector_length():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        TubeDetector.length=-0.1
        assert False
    except AttributeError:
        assert True
    
def test_TubeDetector_diameter():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    try:
        TubeDetector.diameter=-0.1
        assert False
    except AttributeError:
        assert True

def test_TubeDetector1D_plot():
    TubeDetector = TubeDetector1D(position=(0.0,1.0,0.0),direction=(1.0,0,0))
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    TubeDetector.plot(ax)
    