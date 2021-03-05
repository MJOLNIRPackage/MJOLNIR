from __future__ import division
from abc import ABCMeta
import numpy as np
import warnings
# Compability of python 2 and 3 with metaclasses
# Python 2 and 3:
from six import with_metaclass
# or
from future.utils import with_metaclass

def requiredArguments(func,requiredNames):
    """Return list of arguments not found in arguments of func"""
    varnames = func.__code__.co_varnames
    
    notInSignature = np.array([None if x in varnames else x for x in requiredNames])
    notInSignature = notInSignature[notInSignature != np.array(None)]
    return notInSignature

class MaskingObjectMeta(ABCMeta):
    """MetaClass to ensure call signatures of derived classes of MaskingObject"""

    def __new__(mcls,name,bases,namespace): # When a new class is created, that has this metaclass
        missingMethods = [] # Hold all methods missing
        wrongSignature = [] # Hold all methods having wrong signature
        
        if not name in ['MaskingObject','MultiMask']: # If MaskingObject, then just continue
            
            if not '__call__' in namespace: # Check if __call__ exists
                missingMethods.append('__call__')
                
            
            if namespace['dimensionality'] == '2D' or namespace['dimensionality'] == '3D':
                if 'plot' in namespace: # Check if plot has the right signature

                    callFunction = namespace['plot']
                    requiredNames = ['ax','transformation']
                    notInSignature = requiredArguments(callFunction,requiredNames)
                    if len(notInSignature)>0:
                        wrongSignature.append(['plot',requiredNames,notInSignature])
                else:
                    missingMethods.append('plot')
            
            # Write error report
            ErrorMessage = []
            if len(missingMethods)>0:
                ErrorMessage.append("Can't instantiate abstract class {} with abstract methods\n".format(name)+', '.join(missingMethods))
            
            for entry in wrongSignature:
                methodName = entry[0]
                requiredNames = entry[1]
                notInSignature = entry[2]
                ErrorMessage.append('The "{}" method on "{}" must have all of the following arguments in call signature [{}]. Following are missing: [{}]'.format(methodName,name,', '.join(requiredNames),', '.join(notInSignature)))
            if len(ErrorMessage)>0:
                raise TypeError('\n'.join(ErrorMessage))
        
        return ABCMeta.__new__(mcls,name,bases,namespace)
    
            
class MaskingObject(with_metaclass(MaskingObjectMeta)):
    """Base class for all masking objects"""
    #dimensionality = '2D'
    def __init__(self,coordinates=None,maskInside=True,name=None):

        self.maskInside = maskInside
        self._coordinates = []
        if coordinates is None:
            self.coordinates = []
        else:
            if isinstance(coordinates,str):
                coordinates = [coordinates]
            self.coordinates = coordinates
        self.bar = self.negate
        self.name = name

    def plot(self,ax): # pragma: no cover
        """plotting function to put the masked object onto the figure"""
        pass
    
    
    def __call__(self,X,Y=None,Z=None):
        pass
    
    def __add__(self,other):
        return MultiMask(masks=[self,other],operation=np.logical_or)
    
    def __mul__(self,other):
        return MultiMask(masks=[self,other],operation=np.logical_and)
    
    def __neg__(self):
        return self.negate()

    def __pos__(self):
        return self.clone()
    
    def __sub__(self,other):
        return self+(-other)
    
    def __truediv__(self,other):
        return self*(-other)
    
    def generateInputKwargs(self):
        """Generate dictionary with initial args/kwargs used for init"""
        kwargs = {}
        init = self.__class__.__init__
        initCode = init.__code__
        args = initCode.co_argcount
        for X in initCode.co_varnames[1:args]:
            kwargs[X] = self.__dict__[X]
        return kwargs
        
    def clone(self):
        """Return clone of self"""
        return self.__class__(**self.generateInputKwargs())
        
    def negate(self):
        
        newMask = self.clone()
        newMask.maskInside = not newMask.maskInside
        return newMask
        
class MultiMask(MaskingObject):
    def __init__(self,masks,operation=np.logical_and,negated=False):
        if not hasattr(masks,'__len__'):
            masks = [masks]
            
        self.masks = masks
        self.operation = operation
        self.negated = negated
        self.bar = self.negate
        coordinates = np.concatenate([np.array(m.coordinates) for m in self.masks]).flatten()
        self.coordinates = np.array(list(set(coordinates))) # unique coordinates

    
    def __call__(self,*args,**kwargs):
        masks = []
        for m in self.masks:
            masks.append(m(*args,**kwargs))
        
        results = self.operation(*masks)
        if self.negated:
            results = np.logical_not(results)
        return results
    
    def negate(self):
        newMask = self.clone()
        newMask.negated = not newMask.negated
        return newMask
        
        
    def plot(self,ax,transformation=None,*args,**kwargs):# pragma: no cover
        for m in self.masks:
            m.plot(ax,transformation,*args,**kwargs)
    
    def __add__(self,other):
        return MultiMask(masks=[self,other],operation=np.logical_or)
    
    def __mul__(self,other):
        return MultiMask(masks=[self,other],operation=np.logical_and)
    
    def __neg__(self):
        return self.negate()

    def __pos__(self):
        return self.clone()
    
    def __sub__(self,other):
        return self+(-other)
    
    def __truediv__(self,other):
        return self*(-other)
    
    def generateInputKwargs(self):
        """Generate dictionary with initial args/kwargs used for init"""
        kwargs = {}
        init = self.__class__.__init__
        initCode = init.__code__
        args = initCode.co_argcount
        for X in initCode.co_varnames[1:args]:
            kwargs[X] = self.__dict__[X]
        return kwargs
        
    def clone(self):
        """Return clone of self"""
        return self.__class__(**self.generateInputKwargs())



def RotationMatrix(theta,deg=True):
    """Generate 2D rotation matrix given angle theta"""
    if deg == True:
        theta = np.deg2rad(theta)
    M = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    return M

def RotationMatrix3D(theta,n=None,deg=True):
    """Generate 3D rotation matrix given angle and rotation axis"""
    if deg == True:
        theta = np.deg2rad(theta)
    if n is None:
        M = np.eye((3,3))
    else:
        n = np.array(n,dtype=float)
        n*=1.0/np.linalg.norm(n)
        cost = np.cos(theta)
        sint = np.sin(theta)
        M = np.zeros((3,3))
         # taken from http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
        M[0,0] = cost+n[0]**2*(1-cost)
        M[0,1] = n[0]*n[1]*(1-cost)-n[2]*sint
        M[0,2] = n[0]*n[2]*(1-cost)+n[1]*sint
        
        M[1,0] = n[0]*n[1]*(1-cost)+n[2]*sint
        M[1,1] = cost+n[1]**2*(1-cost)
        M[1,2] = n[1]*n[2]*(1-cost)-n[0]*sint
        
        M[2,0] = n[0]*n[2]*(1-cost)-n[1]*sint
        M[2,1] = n[1]*n[2]*(1-cost)+n[0]*sint
        M[2,2] = cost+n[2]**2*(1-cost)
        
    return M

class lineMask(MaskingObject):
    dimensionality = '1D'
    def __init__(self,start,end,maskInside=True,coordinates=None,name=None):
        """Generate a rectangle mask with side corner1 to corner 2 and use corner 3 to define height.
        
        args:
            start (1D point): Starting positon
            
            end (1D point): Ending corner
            
        kwargs:
            
            maskInside (bool): If true, points inside is masked otherwise outside (default True)
            
            coordinates (list of str): List containing names of attributes on point object. If None, use x,y (default None)
            
        Example:
            
            >>> P = np.array([[0.5,2.6])
            >>> b = lineMask(*P)
            >>> mask = b(points.T)
            >>> inside = points[mask].T
            
        """
        super(lineMask,self).__init__(coordinates=coordinates,maskInside=maskInside,name=name)
        
        self.start = start
        self.end = end
        
    
        
    
        
        
    def plot(self,ax,transformation=None,*args,**kwargs):# pragma: no cover
        warnings.warn('It is not possible to plot a 1D masks.')#raise NotImplementedError('It is not possible to plot a 1D masks.')
    
    def __call__(self,x,y=None,z=None):
        if len(self.coordinates)>0:
            points = np.array([getattr(x,coord) for coord in self.coordinates])
        else:
            if y is None:
                #if x.shape[0] != 1:
                #    raise AttributeError('Dimensionality of x is to be (1,N), received {}'.format(x.shape))
                points = x
            else:
                raise AttributeError('Dimensionality of input is to be (1,N), received both an x and y')
        
        mask = np.array([np.logical_and(points>self.start,points<=self.end)])
        
        if self.maskInside == False:
            mask = np.logical_not(mask)
        return mask

class rectangleMask(MaskingObject):
    dimensionality = '2D'
    def __init__(self,corner1,corner2,corner3=None,maskInside=True,coordinates=None,name=None):
        """Generate a rectangle mask with side corner1 to corner 2 and use corner 3 to define height.
        
        args:
            corner1 (2D point): First corner
            
            corner2 (2D point): Second corner
            
        kwargs:
            corner3 (2D point): Third corner used to define height of box. If corner3 is not set corner1 and corner2 are assumed to be opposite corners.
            
            maskInside (bool): If true, points inside is masked otherwise outside (default True)
            
            coordinates (list of str): List containing names of attributes on point object. If None, use x,y (default None)
            
        Example:
            >>> fig,ax = plt.subplots()
            >>> points = np.random.rand(3000,2)
            >>> ax.scatter(points[:,0],points[:,1],color='k')
            >>> P = np.array([[0.15, 0.50],[0.42, 0.16],[0.85, 0.2]])
            >>> b = rectangleMask(*P)
            >>> ax.set_aspect('equal')
            >>> b.plot(ax)
            >>> ax.scatter(*b.center,color='r')
            >>> ax.scatter(P[:,0],P[:,1],color='m')
            >>> mask = b(points.T)
            >>> inside = points[mask].T
            >>> plt.scatter(inside[0],inside[1],color='b')
            >>> plt.scatter(*np.array([b.corners]).T,color='g')
            
        """
        super(rectangleMask,self).__init__(coordinates=coordinates,maskInside=maskInside,name=name)
        
        self.corner1 = np.array(corner1,dtype=float)
        self.corner2 = np.array(corner2,dtype=float)
        
        if isinstance(corner3, np.ndarray):
            if corner3.shape == ():
                self.corner3 = None
            else:
                self.corner3 = corner3
        else:
            if corner3 is None:
                self.corner3 = None
            else:
                self.corner3 = np.array(corner3,dtype=float)
        
        #self.corners = [corner1,corner2]

        
        corner3,corner4 =self.calculateCorners(self.corner1,self.corner2,self.corner3)

        
        self.rotation = np.arctan2(self.edge[1],self.edge[0])
            
        self.rotationMatrix = RotationMatrix(-self.rotation,deg = False)
        
        # Calculate center position 
        self.center = np.mean(self.corners,axis=0).reshape(2,1)
    
    
        
    def calculateCorners(self,c1,c2,c3):
        if c3 is None:
            c3 = c2
            diff = c3-c1
            c2 = c1+np.array([0,diff[1]])
            
        self.edge = c2 - c1 # Defined as line between first and second corner
        self.height = np.linalg.norm(self.edge) # Find side length of rectangle
        self.edge*=1.0/self.height
        self.cross = np.array([-self.edge[1],self.edge[0]]) # Vector orthogonal to edge
        
        self.length = np.dot(self.cross,c3-c1)
        self.lengthSign = np.sign(self.length)
        self.length = np.abs(self.length)
        
        c3 = c2+self.lengthSign*self.cross*self.length/np.linalg.norm(self.cross)
        c4 = c3-c2+c1
        self.corners = [c1,c2,c3,c4]
        return c3,c4
        
        
    def plot(self,ax,transformation=None,*args,**kwargs):# pragma: no cover
        points = np.array([self.corners[0],self.corners[1],self.corners[2],self.corners[3],self.corners[0]]).T
        if not transformation is None:
            points = transformation(*points)
        ax.plot(*points,**kwargs)
    
    def __call__(self,x,y=None,z=None):
        if len(self.coordinates)>0:
            points = np.array([getattr(x,coord) for coord in self.coordinates])
        else:
            if y is None:
                if x.shape[0] != 2:
                    raise AttributeError('Dimensionality of x is to be (2,N), received {}'.format(x.shape))
                points = np.array(x)
            else:
                points = np.array([x,y])
        # calculate mask by rotating all points so rectangle is along coordinate axes
        # This makes the logical checks easy
        pointShape = points.shape[1:]
        points.shape = (2,-1)
        
        rotatedPoints = np.einsum('ij,jk->ik',self.rotationMatrix,np.array(points)-self.center)
        rotatedPoints.shape = np.concatenate([[2],pointShape])
        logic = np.array([np.logical_and(x>-0.5*edge,x<0.5*edge) for x,edge in zip(rotatedPoints,[self.height,self.length])])
        
        mask = np.all(logic,axis=0)
        if self.maskInside == False:
            mask = np.logical_not(mask)
        return mask
    

class boxMask(MaskingObject):
    dimensionality = '3D'
    def __init__(self,corner1,corner2,corner3,corner4=None,maskInside=True,coordinates=None,name=None):
        """Generate a box mask with side corner1, corner 2, corner 3 and use corner 4 to define height.
        
        args:
            corner1 (3D point): First corner
            
            corner2 (3D point): Second corner
            
            corner2 (3D point): Third corner
            
        kwargs:
            corner4 (3D point): Fourth corner used to define height of box. If corner4 is not set corner1, corner2, corner3 are assumed to be maximal positions.
            
            maskInside (bool): If true, points inside is masked otherwise outside (default True)
            
            coordinates (list of str): List containing names of attributes on point object. If None, use x,y (default None)
            
        
        """
        super(boxMask,self).__init__(coordinates=coordinates,maskInside=maskInside,name=name)
        self.corner1 = np.array(corner1,dtype=float)
        self.corner2 = np.array(corner2,dtype=float)
        self.corner3 = np.array(corner3,dtype=float)
        self.corner4 = np.array(corner4,dtype=float)
        
        

        if corner4 is None: # use provided corners and find bounding box
            X,Y,Z = [[np.min(x),np.max(x)] for x in np.array([corner1,corner2,corner3]).T]
            c1 = np.array([X[0],Y[0],Z[0]])
            c2 = np.array([X[1],Y[0],Z[0]])
            c3 = np.array([X[0],Y[1],Z[0]])
            c4 = np.array([X[1],Y[1],Z[0]])
            c5 = np.array([X[0],Y[0],Z[1]])
            c6 = np.array([X[1],Y[0],Z[1]])
            c7 = np.array([X[0],Y[1],Z[1]])
            c8 = np.array([X[1],Y[1],Z[1]])
            
            self.length = X[1]-X[0]
            self.width = Y[1]-Y[0]
            self.height = Z[1]-Z[0]
            
            self.planeNormal = np.array([0,0,1])#np.cross(edge1,edge2)
            #planeNormal *= 1.0/np.linalg.norm(planeNormal)
            
            #self.rotationMatrix = RotationMatrix3D(0.0,n=None,deg=False)
            
        else:
            c1 = np.array(corner1,dtype=float)
            c2 = np.array(corner2,dtype=float)
            
            # Start by finding common plane for c1, c2, and  c3
            edge1 = c2 - c1 # Defined as line between first and second corner
            edge2 = corner3 - c1 # Defined as line between first and second corner
            
            self.length = np.linalg.norm(edge1)
            
            edge1 *= 1.0/self.length
            edge2 *= 1.0/np.linalg.norm(edge2)
            
            planeNormal = np.cross(edge1,edge2)
            self.planeNormal = planeNormal/np.linalg.norm(planeNormal)
            
            cross = np.cross(edge1,self.planeNormal)
            
            
            width = np.dot(cross,corner3 - c1) # distance from c3 orthogonally on edge1
            widthSign = np.sign(width)
            self.width = np.abs(width)
            
            
            # calculate corner point for c3 and c4
            c4 = c2+widthSign*cross*self.width
            
            c3 = c4-c2+c1 # Find corresponding c4
            
            # find cube
            height = np.dot(self.planeNormal,corner4-c1)
            heightSign = np.sign(height)
            self.height = np.abs(height)
            
            c5 = c1+heightSign*self.planeNormal*self.height
            c6 = c5+c2-c1 
            c7 = c5+c3-c1 
            c8 = c5+c4-c1 
            
        self.corners = np.array([c1,c2,c3,c4,c5,c6,c7,c8])
        self.center = np.array([np.mean(X) for X in self.corners.T]).reshape(1,-1)
        
        # Create rotation matrix.

        # rotate planeNormal to lie along x,z plane (https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d)
        
        v = np.cross(self.planeNormal,np.array([0.0,1.0,0.0]))
        Vs = np.zeros((3,3))
        Vs[0,1] = -v[2]
        Vs[0,2] = +v[1]
        Vs[1,2] = -v[0]
        Vs[1,0] = -Vs[0,1]
        Vs[0,2] = -Vs[0,2]
        Vs[2,1] = -Vs[1,2]
        
        c = np.dot(self.planeNormal,np.array([0.0,1.0,0.0]))
        
        rot1 = np.eye(3)+Vs+np.dot(Vs,Vs)*1.0/(1+c)
        
        
        c1r1,c2r1,c3r1,c4r1,c5r1,c6r1,c7r1,c8r1 = [np.dot(rot1,(x-self.center).T) for x in [c1,c2,c3,c4,c5,c6,c7,c8]]
        
        # rotate around y into plane using edge 1 (c2-c1)
        edge = c2r1-c1r1
        theta = np.arctan2(edge[2],edge[0])
        
        rot2 = RotationMatrix3D(theta,n=[0,1,0],deg=False)
        self.rotationMatrix = np.einsum('ij,jk->ik',rot2,rot1)

        
        
    def plot(self,ax,transformation = None, *args,**kwargs):# pragma: no cover
        
        # Draw square for Z = 0
        points = np.array([self.corners[0],self.corners[1],self.corners[3],self.corners[2],self.corners[0]]).T
        if not transformation is None:
            points = transformation(*points)
        p = ax.plot3D(*points,**kwargs)
        if not 'color' in kwargs:
            color = p[0].get_color()
            kwargs['color'] = color
        
        # Draw square for Z = 1
        points = np.array([self.corners[4],self.corners[5],self.corners[7],self.corners[6],self.corners[4]]).T
        if not transformation is None:
            points = transformation(*points)
        ax.plot3D(*points,**kwargs)
        
        # Draw square for X = 0
        points = np.array([self.corners[0],self.corners[2],self.corners[6],self.corners[4],self.corners[0]]).T
        if not transformation is None:
            points = transformation(*points)
        ax.plot3D(*points,**kwargs)
        # Draw square for X = 1
        points = np.array([self.corners[1],self.corners[3],self.corners[7],self.corners[5],self.corners[1]]).T
        if not transformation is None:
            points = transformation(*points)
        ax.plot3D(*points,**kwargs)
    
    
    def __call__(self,x,y=None,z=None):
        if len(self.coordinates)>0:
            points = np.array([getattr(x,coord) for coord in self.coordinates])
        else:
            if y is None:
                if x.shape[0] != 3:
                    raise AttributeError('Dimensionality of x is to be (3,N), received {}'.format(x.shape))
                points = x
            else:
                if z is None:
                    raise AttributeError('Recieved only X and Y, but no Z')
                points = np.array([x,y,z])
        # calculate mask by rotating all points so rectangle is along coordinate axes
        # This makes the logical checks easy
        pointShape = points.shape[1:]
        
        points.shape = (3,-1)
        
        rotatedPoints = np.einsum('ij,jk->ik',self.rotationMatrix,np.array(points)-self.center.T)
        rotatedPoints.shape = np.concatenate([[3],pointShape]).astype(int)
        logic = np.array([np.logical_and(x>=-0.5*edge,x<=0.5*edge) for x,edge in zip(rotatedPoints,[self.length,self.height,self.width])])
        
        mask = np.all(logic,axis=0)
        if self.maskInside == False:
            mask = np.logical_not(mask)
        return mask


class circleMask(MaskingObject):
    
    dimensionality = '2D'
    
    def __init__(self,center,radiusPoint=None,radius=None,maskInside=True,coordinates=None,name=None):
        """Generate a circular mask with center at center and edge at radiusPoint.
        
        args:
            center (2d point): center
            
        kwargs:
            radiusPoint (2d point): Circumference goes through this point (default None)

            radius (float): radius if  circle/sphere. If set, it overrides radiusPoint (default none)

            maskInside (bool): If true, points inside is masked otherwise outside (default True)
            
            coordinates (list of str): List containing names of attributes on point object. If None, use x,y (default None)
            
        Example:
            >>> fig,ax = plt.subplots()
            >>> points = np.random.rand(3000,2)
            >>> ax.scatter(points[:,0],points[:,1],color='k')
            >>> P = np.array([[0.15, 0.50],[0.42, 0.16]])
            >>> b = circleMask(*P)
            >>> ax.set_aspect('equal')
            >>> b.plot(ax)
            >>> ax.scatter(*b.center,color='r')
            >>> ax.scatter(P[:,0],P[:,1],color='m')
            >>> mask = b(points.T)
            >>> inside = points[mask].T
            >>> plt.scatter(inside[0],inside[1],color='b')
            >>> plt.scatter(*np.array([b.corners]).T,color='g')
            
        """
        super(circleMask,self).__init__(coordinates=coordinates,maskInside=maskInside,name=name)
        
        self.center = np.array(center,dtype=float).reshape(-1,1)
        if self.center.shape[0] == 3: # 3D sphere
            self.dimensionality = '3D'

        if not radius is None:
            radiusVector= np.zeros_like(self.center)
            radiusVector[0] = radius
            radiusPoint = self.center + radiusVector
        else:
            if radiusPoint is None:
                raise AttributeError('Either radius or radiusPoint is to be set.')
        radiusPoint = np.array(radiusPoint,dtype=float).reshape(-1,1)
        self.radius = np.linalg.norm(radiusPoint-self.center)
        self.radiusPoint = radiusPoint

        
        
        
    def plot(self,ax,transformation=None,**kwargs):# pragma: no cover
        if self.dimensionality == '2D':
            theta = np.linspace(0,np.pi*2,200)
            
            points = np.array([np.cos(theta),np.sin(theta)])*self.radius+self.center
            if not transformation is None:
                points = transformation(*points)
            ax.plot(points[0],points[1],**kwargs)
        elif self.dimensionality == '3D':
            
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = np.cos(u)*np.sin(v)*self.radius+self.center[0]
            y = np.sin(u)*np.sin(v)*self.radius+self.center[1]
            z = np.cos(v)*self.radius+self.center[2]
            points = [x,y,z]
            if not transformation is None:
                points = transformation(*points)
            ax.plot_wireframe(*points, **kwargs)
    
    def __call__(self,x,y=None,z=None):
        if len(self.coordinates)>0:
            points = np.array([getattr(x,coord) for coord in self.coordinates])
        else:
            if self.dimensionality == '2D':
                if y is None:
                    if x.shape[0] != 2:
                        raise AttributeError('Dimensionality of x is to be (2,N), received {}'.format(x.shape))
                    points = x
                else:
                    points = np.array([x,y])
            elif self.dimensionality == '3D':
                if y is None or z is None:
                    if x.shape[0] != 3:
                        raise AttributeError('Dimensionality of x is to be (3,N), received {}'.format(x.shape))
                    points = x
                else:
                    points = np.array([x,y,z])
        
        pointShape = points.shape[1:]
        points.shape = (points.shape[0],-1)
        
        mask = np.linalg.norm(points-self.center,axis=0)<=self.radius
        mask.shape = pointShape
        
        if self.maskInside == False:
            mask = np.logical_not(mask)
        return mask
    
class indexMask(MaskingObject):
    dimensionality = '1D'
    
    def __init__(self,start,end,axis=0,maskInside=True,name=None):
        """Mask indices along a given axis.
        
        args:
            start (int): starting point included in mask
        
            end (int): end point not included in mask
            
        kwargs:
            axis (int): Axis along which masking is to be performed (default 0)
            
            maskInside (bool): If true, points inside is masked otherwise outside (default True)
        
        """
        super(indexMask,self).__init__(coordinates=None,maskInside=maskInside,name=name)
        self.axis = axis
        self.start = start
        self.end = end
        
    def plot(self,ax,transformation=None,*args,**kwargs):# pragma: no cover
        raise NotImplementedError('It is not possible to plot a 1D masks.')
    
    def __call__(self,x,*args):
        points = x#np.concatenate([[x],args])
        
        if not len(points.shape)>self.axis:
            raise AttributeError('Masking axis is {}, but shape of x is {}.'.format(self.axis,points))
        
        arange = np.arange(points.shape[self.axis]).flatten()
       
        transpose = np.arange(len(points.shape))
        transpose[self.axis] = 0
        transpose[0] = self.axis
        
        otherAxes = int(points.size/points.shape[self.axis])
        inside = np.array([np.logical_and(arange>=self.start,arange<self.end)]).flatten()
                
        # create mask with masking axis along 0
        newShape = np.array(points.shape)[transpose]
        mask = np.repeat(inside,otherAxes).reshape(newShape).transpose(transpose)
        
        
        if self.maskInside == False:
            mask = np.logical_not(mask)
        return mask


        
