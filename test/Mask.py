from xml.dom.minidom import Attr
import numpy as np
import warnings
# Compability of python 2 and 3 with metaclasses
# Python 2 and 3:
from six import with_metaclass
# or
from future.utils import with_metaclass
from MJOLNIR.Data.Mask import MaskingObject, lineMask, rectangleMask, circleMask, boxMask, indexMask, MultiMask,CurratAxeMask


def test_subclass_MaskingObject():
    # Generate a subclass of MaskingObect that is missing stuff
    
    try:
        class Missing__call__(MaskingObject):
            dimensionality = '2D'
            def __init__(self,*args,**kwargs):
                pass
            
            def plot(self,ax,transformation=None,**kwargs):
                pass
        assert False
    except TypeError as E:
        print(E)
        assert True
    
    try:
        class Missing__plot__(MaskingObject):
            dimensionality = '2D'
            def __init__(self,*args,**kwargs):
                pass
            
            def __call__(self,ax,**kwargs):
                pass
        assert False
    except TypeError as E:
        print(E)
        assert True
        
    try:
        class Wrong__plot__(MaskingObject):
            dimensionality = '2D'
            def __init__(self,*args,**kwargs):
                pass
            
            def __call__(self,ax,**kwargs):
                pass
            
            def plot(self,ax,**kwargs): # missing transformation argument
                pass
        assert False
    except TypeError as E:
        print(E)
        assert True
        
def test_BooleanAlgebra():
    class simpleMaskingObject(MaskingObject):
        dimensionality = '1D'
        def __init__(self,coordinates=None,maskInside=True):
            super(simpleMaskingObject,self).__init__(coordinates=coordinates,maskInside=maskInside)
            
        
        def call(self):
            return self.maskInside
        
    true =  simpleMaskingObject(maskInside=True)
    false = simpleMaskingObject(maskInside=False)
    
    # Adding gives a MultiMask
    multi = true+true
    assert(isinstance(multi,MultiMask))
    assert(multi()==True)
    
    negated = -multi
    assert(negated()==False)
    
    multi2 = true*true
    assert(multi2()==True)
    negated = -multi2
    assert(negated()==False)
    
    assert(negated() == false())
    
    negated = -false
    
    assert(negated()==True)
    
    assert((false-false)() == True)
    
    # Check 'division' i.e. oposite of *
    #   A   B   Op  Res
    #   1   1   /    0
    #   1   0   /    1
    #   0   1   /    0
    #   0   0   /    0
    assert((true/true  )() == False)
    assert((true/false )() == True)
    assert((false/true )() == False)
    assert((false/false)() == False)
    
    
    
    ## Following test check if expressions can be concatenated together
    # Generate truth table using + and *
    #
    #   A   B   Op  Res
    #   1   1   +    1
    #   1   1   *    1
    #   1   0   +    1
    #   0   1   +    1
    #   1   0   *    0
    #   0   1   *    0
    #   0   0   *    0
    #   0   0   +    0
    
    assert((true+true  )()==True)
    assert((true*true  )()==True)
    assert((true+false )()==True)
    assert((false+true )()==True)
    assert((true*false )()==False)
    assert((false*true )()==False)
    assert((false*false)()==False)
    assert((false+false)()==False)
    
    
    # Generate truth table using + and * for expression (A OP1 B) OP2 C
    #
    #  (A   OP1   B)  OP2   C   Res
    #  (1    *    0)   +    1    1
    #  (0    +    0)   *    1    0
    #  (1    *    1)   *    1    1
    
    val1 = (true*true)+true
    val2 = (false+false)*true
    val3 = (true*true)*true
    
    assert(val1()==True)
    assert(val2()==False)
    assert(val3()==True)
    
    # Generate truth table using + and * for expression (A OP1 B) OP2 (C OP3 D)
    #
    #  (A   OP1   B)  OP2  (C   OP1   D)  Res
    #  (1    *    0)   +   (1    +    1)   1
    #  (0    +    0)   *   (1    *    0)   0
    #  (0    +    0)   +   (1    *    0)   0
    #  (1    *    0)   *   (1    +    1)   0
    
    val1 = (true*true)+(true+true)
    val2 = (false+false)*(true*false)
    val3 = (false+false)+(true*false)
    val4 = (true*false)*(true+true)
    
    assert(val1()==True)
    assert(val2()==False)
    assert(val3()==False)
    assert(val4()==False)
    
def test_lineMask():
    X,Y,Z = np.array(np.meshgrid(np.linspace(0,1,50),np.linspace(0,1,50),np.linspace(0,1,50))).reshape(3,-1)
    
    
    class p(object):
        pass
    
    points = p()
    
    for name,coord in zip(['X','Y','Z'],[X,Y,Z]):
        setattr(points,name,coord)
        
    line = lineMask(0.5,1.0,coordinates='X')
    line2 = lineMask(0.5,1.0)
    line3 = lineMask(-0.001,0.50,maskInside=False)
    mask = line(points)
    mask2 = line2(X)
    mask3 = line3(X)
    print(np.sum(mask))
    print(np.sum(mask2))
    print(np.sum(mask3))
    assert(np.all(mask==mask2))
    assert(np.all(mask3==mask2))
    
    
def test_rectangleMask():
    X,Y,Z = np.array(np.meshgrid(np.linspace(0,3,50),np.linspace(0,3,50),np.linspace(0,3,50))).reshape(3,-1)
    
    class p(object):
        pass
    
    points = p()
    
    for name,coord in zip(['X','Y','Z'],[X,Y,Z]):
        setattr(points,name,coord)
        
    rec1 = rectangleMask(corner1=[0,2],corner2=[1,0],coordinates=['X','Y'])
    rec2 = rectangleMask(corner1=[0,2],corner2=[0,0],corner3=[1,0])
    
    mask1 = rec1(points)
    mask2 = rec2(X,Y)
    assert(np.all(mask1==mask2))
    assert(rec1.length == 1.0)
    assert(rec1.height == 2.0)
    
    # rotated mask
    X = np.ones(100)*0.5
    Y = np.linspace(-1,1,100)
    points = np.array([X,Y])
    
    # Points are inside  if y is less than 0.5 and bigger than -0.5
    rec2 = rectangleMask(corner1=[1,1],corner2=[0,0],corner3=[1,-1])
    assert(np.isclose(rec2.rotation,-0.75*np.pi))
    inside = np.logical_and(Y<0.5,Y>-0.5)
    mask = rec2(X,Y)
    assert(np.all(inside==mask))
    
    rec3= rectangleMask(corner1=[2,1],corner2=[0,0],corner3=[1,-2])
    assert(rec3.height == rec3.length)
    assert(rec3.length == 2.23606797749979)
    assert(np.isclose(rec3.rotation,-0.85241638*np.pi))
    
    
def test_circleMask():
    X,Y,Z = np.array(np.meshgrid(np.linspace(0,3,50),np.linspace(0,3,50),np.linspace(0,3,50))).reshape(3,-1)
    
    class p(object):
        pass
    
    points = p()
    
    for name,coord in zip(['X','Y','Z'],[X,Y,Z]):
        setattr(points,name,coord)
        
    try:
        circ0 = circleMask(center=[0.5,0.5])
        assert False
    except AttributeError:
        assert True

    circ1 = circleMask(center=[0.5,0.5],radiusPoint=[1,0.5],coordinates=['X','Y'])
    circ2 = circleMask(center=[0.5,0.5],radiusPoint=[1,0.5])
    circ3 = circleMask(center=[0.5,0.5],radius=0.5)

    assert(np.isclose(circ1.radius,0.5))
    assert(np.isclose(circ3.radius,0.5))
    mask = circ1(points)
    assert(np.sum(mask)==10450)
    assert(np.all(mask==circ2(X,Y)))
    assert(np.all(mask==circ3(X,Y)))
    assert(circ2(0.5,0.5))
    assert(circ2(0.9,0.9) == False)
    assert(circ2(np.cos(np.pi/4)*0.5+0.501,np.cos(np.pi/4)*0.5+0.501) == False)
    assert(circ2(np.cos(np.pi/4)*0.5+0.499,np.cos(np.pi/4)*0.5+0.499) == True)
    
def test_boxMask():
    X,Y,Z = np.array(np.meshgrid(np.linspace(0,2,100),np.linspace(0,2,100),np.linspace(0,2,100))).reshape(3,-1)
    
    class p(object):
        pass
    
    points = p()
    
    for name,coord in zip(['X','Y','Z'],[X,Y,Z]):
        setattr(points,name,coord)
        
    corners = np.array([[0,0,0],[1,0,0],[0,0.5,0],[0,0,0.5]])
    box = boxMask(*corners,coordinates=['X','Y','Z'])
    box2= boxMask([0,0,0],[1,0.5,0],[0,0,0.5])
    mask = box(points)
    mask2 = box2(X,Y,Z)
    assert(np.sum(mask)==int(100**3*(1/2*1/4**2)))
    assert(np.all(mask2==mask))
    
    box3 = boxMask([1,0,0],[0,0.5,0],[0,0,0.5],[0,0,0])
    assert(np.isclose(box3.length*box3.width*box3.height,0.25))
    assert(np.all(np.isclose(box3.center,np.array([ 0.34444444, -0.06111111,  0.13888889]))))
    
    
def test_indexMask():
    A = np.arange(120).reshape(4,5,6)
    imask = indexMask(1,4,axis=1) # mask all all but two edges of first axis
    mask = imask(A)
    assert(np.all(mask[:,0,:]==False))
    assert(np.all(mask[:,1,:]==True))
    assert(np.all(mask[:,2,:]==True))
    assert(np.all(mask[:,3,:]==True))
    assert(np.all(mask[:,4,:]==False))
    
    A = np.arange(120*3*7).reshape(4,5,6,3,7)
    imask = indexMask(1,2,axis=3) # mask all all but two edges of first axis
    mask = imask(A)
    assert(np.all(mask[:,:,:,0,:]==False))
    assert(np.all(mask[:,:,:,1,:]==True))
    assert(np.all(mask[:,:,:,2,:]==False))
    
def test_CurratAxeMask():
    # Real test is performed in the DataSet testing

    M = CurratAxeMask([[1,0,0]])
    X = np.linspace(0,1,11)
    
    try:
        M(X)
        assert False
    except AttributeError:
        assert True

def test_CurratAxeMask():
    # Real test is performed in the DataSet testing

    M = CurratAxeMask([[1,0,0]])
    X = np.linspace(0,1,11)
    
    try:
        M(X)
        assert False
    except AttributeError:
        assert True
