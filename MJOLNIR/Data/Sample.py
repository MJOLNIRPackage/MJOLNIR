import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import scipy
import matplotlib.pyplot as plt
import numpy as np
import h5py as hdf
import warnings
from MJOLNIR import _tools
import datetime
import math
from MJOLNIR import TasUBlibDEG as TasUBlib
from MJOLNIR._tools import Marray
import MJOLNIR.Data.DataFile


def cosd(x):
    return np.cos(np.deg2rad(x))

def sind(x):
    return np.sin(np.deg2rad(x))

def camelCase(string,split='_'):
    """Convert string to camel case from <split> seperated"""

    if not split in string:
        return string
    splitString = string.split(split)
    first = splitString[0]
    others = [x.title() for x in splitString[1:]]
    
    combi = [first]+others
    return ''.join([str(x) for x in combi])

class Sample(object):
    """Sample object to store all information of the sample from the experiment"""
    @_tools.KwargChecker()
    def __init__(self,a=1.0,b=1.0,c=1.0,alpha=90,beta=90,gamma=90,sample=None,name='Unknown',projectionVector1=None, projectionVector2 = None):
        if isinstance(sample,hdf._hl.group.Group):
            self.name = str(np.array(sample.get('name'))[0].decode())
            self.orientationMatrix = np.array(sample.get('orientation_matrix'))*2*np.pi

            self.planeNormal = np.array(sample.get('plane_normal'))
            
            self.polarAngle = np.array(sample.get('polar_angle'))
            self.rotationAngle = np.array(sample.get('rotation_angle'))
            self.unitCell = np.array(sample.get('unit_cell'))
            self.plane_vector1 = np.array(sample.get('plane_vector_1'))
            self.plane_vector2 = np.array(sample.get('plane_vector_2'))
            crossProduct = np.cross(self.plane_vector1[:3],self.plane_vector2[:3])
            if not np.all(np.isclose(crossProduct,[0,0,0])):
                self.planeNormal = crossProduct
            self.A3Off = np.array([0.0])#
            if not np.isclose(np.linalg.norm(self.plane_vector1[:3].astype(float)),0.0) or not np.isclose(np.linalg.norm(self.plane_vector2[:3].astype(float)),0.0): # If vectors are not zero
                self.projectionVector1,self.projectionVector2 = calcProjectionVectors(self.plane_vector1.astype(float),self.plane_vector2.astype(float))
            else:
                self.projectionVector1,self.projectionVector2 = [np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0])]
            self.initialize()
            self.calculateProjections()

            attributes = ['azimuthal_angle','x','y','sgu','sgu_zero','sgl','sgl_zero']
            values = [camelCase(x) for x in attributes]
            for att, val in zip(attributes,values):
                setattr(self,val,np.array(sample.get(att)))
            
            
        elif np.all([a is not None,b is not None, c is not None]):
            self.unitCell = np.array([a,b,c,alpha,beta,gamma])
           
            self.polarAngle = np.array(None)
            self.rotationAngle = np.array(0)
            self.name=name
            if projectionVector1 is None or projectionVector2 is None:
                projectionVector1,projectionVector2 = [np.array([1.0,0.0,0.0,0.0,0.0,0.0,0.0,5.0,5.0]),np.array([0.0,1.0,0.0,0.0,0.0,0.0,0.0,5.0,5.0])]

            r1 = projectionVector1
            r2 = projectionVector2
            self.plane_vector1 = r1
            self.plane_vector2 = r2

            self.planeNormal = np.cross(self.plane_vector1[:3],self.plane_vector2[:3])
            
            
            cell = TasUBlib.calcCell(self.unitCell)

            self.orientationMatrix = TasUBlib.calcTasUBFromTwoReflections(cell, r1, r2)
            #self.orientationMatrix = TasUBlib.calcTasUBFromTwoReflections(self.cell,self.plane_vector1,self.plane_vector2)
            self.projectionVector1,self.projectionVector2 = calcProjectionVectors(self.plane_vector1.astype(float),self.plane_vector2.astype(float))#,self.planeNormal.astype(float))
            self.initialize()
            self.calculateProjections()
        else:
            print(sample)
            print(a,b,c,alpha,beta,gamma)
            raise AttributeError('Sample not understood')

            
    @property
    def unitCell(self):
        return self._unitCelll

    @unitCell.getter
    def unitCell(self):
        return np.array([self.a,self.b,self.c,self.alpha,self.beta,self.gamma])#self._unitCell

    @unitCell.setter
    def unitCell(self,unitCell):
        self._unitCell = unitCell
        self.a = unitCell[0]
        self.b = unitCell[1]
        self.c = unitCell[2]
        self.alpha = unitCell[3]
        self.beta  = unitCell[4]
        self.gamma = unitCell[5]
        self.updateCell()
        
    @property
    def a(self):
        return self._a

    @a.getter
    def a(self):
        return self._a

    @a.setter
    def a(self,a):
        if a>0:
            self._a = a
        else:
            raise AttributeError('Negative or null given for lattice parameter a')

    @property
    def b(self):
        return self._b

    @b.getter
    def b(self):
        return self._b

    @b.setter
    def b(self,b):
        if b>0:
            self._b = b
        else:
            raise AttributeError('Negative or null given for lattice parameter b')

    @property
    def c(self):
        return self._c

    @c.getter
    def c(self):
        return self._c

    @c.setter
    def c(self,c):
        if c>0:
            self._c = c
        else:
            raise AttributeError('Negative or null given for lattice parameter c')


    @property
    def alpha(self):
        return self._alpha

    @alpha.getter
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self,alpha):
        if alpha>0 and alpha<180:
            self._alpha = alpha
        else:
            raise AttributeError('Negative,null or above 180 degrees given for lattice parameter alpha')

    @property
    def beta(self):
        return self._beta

    @beta.getter
    def beta(self):
        return self._beta

    @beta.setter
    def beta(self,beta):
        if beta>0 and beta<180:
            self._beta = beta
        else:
            raise AttributeError('Negative,null or above 180 degrees given for lattice parameter beta')

    @property
    def gamma(self):
        return self._gamma

    @gamma.getter
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self,gamma):
        if gamma>0 and gamma<180:
            self._gamma = gamma
        else:
            raise AttributeError('Negative,null or above 180 degrees given for lattice parameter gamma')


    def __eq__(self,other):
        if not isinstance(other,type(self)):
            return False
        return np.all(self.unitCell==other.unitCell)#,\
        #np.all(self.orientationMatrix==other.orientationMatrix)])

    def initialize(self):
        """Initialize the Sample object. Automatically called during __init__method."""
        # From http://gisaxs.com/index.php/Unit_cell
        self.updateCell()
        self.B = TasUBlib.calculateBMatrix(self.cell)
        #self.reciprocalMatrix = np.array([self.reciprocalVectorA,self.reciprocalVectorB,self.reciprocalVectorC]).T

    def calculateProjections(self):
        """Calculate projections and generate projection angles."""
        checks = np.array(['unitCell','orientationMatrix','projectionVector1','projectionVector2']) # 
        boolcheck = np.logical_not(np.array([hasattr(self,x) for x in checks]))
        if np.any(boolcheck):
            raise AttributeError('Sample object is missing: {}.'.format(', '.join(str(x) for x in checks[boolcheck])))
        
        if self.projectionVector1[np.argmax(np.abs(self.projectionVector1))]<0:
            self.projectionVector1*=-1
        if self.projectionVector2[np.argmax(np.abs(self.projectionVector2))]<0:
            self.projectionVector2*=-1

        V1 = self.projectionVector1.copy()
        #V1/=np.linalg.norm(V1)
        V2 = self.projectionVector2.copy()
        #V2/=np.linalg.norm(V2)
        pV1Q = np.dot(self.B,V1)
        pV2Q = np.dot(self.B,V2)
        self.projectionAngle = _tools.vectorAngle(pV1Q,pV2Q)

        if np.isclose(0.0,self.projectionAngle):
            raise AttributeError("The provided orientations are equal.")

        
        UB= self.orientationMatrix
        self.UB = UB

        self.orientationMatrixINV = np.linalg.inv(UB)
        

        p23 = np.array([[1,0,0],[0,1,0]]) # To extract Qx, Qy only
        PM = np.array([V1,V2]).T # Projection matrix
        self.PM = PM
        self.convert = np.dot(p23,np.einsum('ij,jk->ik',UB,PM)) # Convert from projX,projY to Qx, Qy
        self.convertHKL = np.dot(p23,UB) # Convert from HKL to Qx, Qy

        # Calculate 'misalignment' of the projection vector 1
        try:
            self.theta = -TasUBlib.calcTasMisalignment(UB,self.planeNormal,V1)
        except AttributeError:
            self.theta = 0
        
        self.RotMat = _tools.Rot(self.theta) # Create 2x2 rotation matrix
        
        self.convertinv = np.linalg.inv(self.convert) # Convert from Qx, Qy to projX, projY

        self.convertHKLINV = _tools.invert(self.convertHKL) # Convert from Qx, Qy to HKL

        # Calcualte RotationMatrix for UB as block diagonal matrix. Only Qx,Qy part rotates as'
        # calculated in RotMat
        self.RotMat3D = np.eye(3)
        self.RotMat3D[:2,:2] = self.RotMat
        #self.orientationMatrixINV = np.linalg.inv(np.dot(self.RotMat3D,UB))
        self.orientationMatrixINV = np.linalg.inv(self.UB)
        

    def tr(self,p0,p1):
        """Convert from projX, projY coordinate to Qx,QY coordinate."""
        p0, p1 = np.asarray(p0), np.asarray(p1)
        
        P = np.array([p0,p1])
        
        Pos = np.einsum('ij,j...->i...',self.convertinv,P)
        return Pos[0],Pos[1]
        
    def inv_tr(self, x,y):
        """Convert from Qx,QY  coordinate to projX, projY coordinate."""
        x, y = np.asarray(x), np.asarray(y)
        P = np.array([x,y])
        Pos = np.einsum('ij,j...->i...',self.convert,P)
        return Pos[0],Pos[1]   


    def format_coord(self,x,y):
        """Format coordinates from QxQy in rotated frame into HKL."""
        x, y = np.asarray(x), np.asarray(y)
        rlu = self.calculateQxQyToHKL(x,y)#np.dot(self.orientationMatrixINV,np.array([x,y,0]))
        return "h = {0:.3f}, k = {1:.3f}, l = {2:.3f}".format(rlu[0],rlu[1],rlu[2])

    def calculateQxQyToHKL(self,x,y):
        """convert from Qx,Qy to HKL."""
        pos = np.array([x,y,np.zeros_like(x)])
        return np.einsum('ij,j...->i...',self.orientationMatrixINV,pos)

    def calculateHKLToQxQy(self,H,K,L):
        """convert HKL to Qx,Qy."""
        pos = np.array([H,K,L])
        return np.einsum('ij,j...->i...',self.orientationMatrix,pos)[:2]

    def calculateHKLtoProjection(self,H,K,L):
        """convert from projections to HKL."""
        HKL = np.array([H,K,L])
        #points = np.einsum('i...,ij...->i...',HKL,self.PM)
        points = np.einsum('ij,j...->i...',_tools.invert(self.PM),HKL)
        return points


    def __str__(self):
        returnStr = 'Sample ' + self.name + '\n'
        #if not self.temperature is None: returnStr+= 'Temperatur: '+str(self.temperature)+'\n'
        #if not self.magneticField is None: returnStr+= 'Magnetic Field: '+str(self.magneticField)+'\n'
        #if not self.electricField is None: returnStr+= 'Electric Field: '+str(self.electricField)+'\n'
        returnStr+= 'Unit cell: \n' + str(self.unitCell) + '\n'
        returnStr+= 'Orientation matrix: \n' + str(self.orientationMatrix) +'\n'

        return returnStr

    @_tools.KwargChecker()
    def CurratAxe(self,Ei,Ef,Bragg,spurionType='Monochromator',HKL=False,Projection=False):
        """Function to calculate Currat-Axe position in QxQy coordinate system.
    
        Args:
            
            - Ei (float/list): Incoming energy in meV.
            
            - Ef (float/list): Outgoing energy in meV.
            
            - Bragg (list): Bragg peak in HKL or list of.
            
        Kwargs:
            
            - spurionType (str): Either "Monochromator" or "Analyser" for origin of "wrong" energy (default "Monochromator").
            
            - HKL (bool): Whether or not to recalculate to HKL instead of Qx, Qy, Qz (default False).
            
            - Projection (Bool): Whether or not to recalculate to Projection vectors instead of Qx, Qy, Qz (default False).
            
        Returns:
            
            - Position (list): List of size (len(Bragg),len(Ei),len(Ef),3), where last axis is Qx, Qy, Qz
            
        """
        A3Off = self.theta
        
        Ei = np.asarray(Ei).flatten() # Shape (m)
        Ef = np.asarray(Ef).flatten() # Shape (n)
        Bragg = np.asarray(Bragg).reshape(-1,3) # shape (l,3)
        
        QlLocal = []
        if spurionType.lower() == 'monochromator':
            for B in Bragg:
                Ql = []
                Angles = np.array([TasUBlib.calcTasQAngles(self.orientationMatrix,self.planeNormal,1.0,A3Off,np.array([B[0],B[1],B[2],e,e]))[:2] for e in Ef])
                for ei in Ei:
                    Ql.append(np.array([TasUBlib.calcTasQH(self.orientationMatrixINV,angle,ei,e,0) for angle,e in zip(Angles,Ef)])[:,1])
                QlLocal.append(Ql)
        elif spurionType.lower() == 'analyser':
            for B in Bragg:
                Ql = []
                for ei in Ei:
                    Angles2 = np.array(TasUBlib.calcTasQAngles(self.orientationMatrix,self.planeNormal,1.0,A3Off,np.array([B[0],B[1],B[2],ei,ei]))[:2])
                    Ql.append(np.array([TasUBlib.calcTasQH(self.orientationMatrixINV,Angles2,ei,e,A3Off*0.0) for e in Ef])[:,1]) # Extract Qx,Qy
                QlLocal.append(Ql)
        else:
            raise AttributeError('Provided spurionType not understood. Expected "Monochromator" or "Analyser" but recieved "{}".'.format(spurionType))

        returnVal = np.array(QlLocal) # Shape (l,m,n,3)
        
        if HKL == True or Projection == True: # need to calculate HKL for Projection calculation
            returnValShape = np.array(returnVal.shape)
            returnVal = self.calculateQxQyToHKL(returnVal[:,:,:,0].flatten(),returnVal[:,:,:,1].flatten())
            if Projection == True:
                
                toProjection = self.calculateHKLtoProjection(returnVal[0],returnVal[1],returnVal[2])
                returnVal = np.array(self.inv_tr(toProjection[0],toProjection[1]))
                returnValShape[-1]=2 # reshape Qx,Qy,Qz dimension to P1,P2 (3 -> 2)
            
            returnVal.shape = returnValShape # Shape (l,m,n,3) or (l,m,n,2)
        return returnVal

    def updateCell(self,unitCell=None):
        """Update cell parameters with current unit cell values.

        Kwargs:
            
            - unitCell (list): List of a,b,c,alpha,beta,gamma. If None, use self.unitCell (default None)

        """
        if unitCell is None:
            unitCell = self.unitCell
        else:
            self.unitCell = unitCell

        self.realVectorA = np.array([self.a,0,0])
        self.realVectorB = self.b*np.array([cosd(self.gamma),sind(self.gamma),0.0])#np.dot(np.array([self.b,0,0]),rotationMatrix(0,0,self.gamma))
        self.realVectorC = self.c*np.array([cosd(self.beta),(cosd(self.alpha)-cosd(self.beta)*cosd(self.gamma))/sind(self.gamma),
        np.sqrt(1-cosd(self.beta)**2-((cosd(self.alpha)-cosd(self.beta)*cosd(self.gamma))/sind(self.gamma))**2)])#np.dot(np.array([self.c,0,0]),rotationMatrix(0,self.beta,0))
        
        self.volume = np.abs(np.dot(self.realVectorA,np.cross(self.realVectorB,self.realVectorC)))
        self.reciprocalVectorA = 2*np.pi*np.cross(self.realVectorB,self.realVectorC)/self.volume
        self.reciprocalVectorB = 2*np.pi*np.cross(self.realVectorC,self.realVectorA)/self.volume
        self.reciprocalVectorC = 2*np.pi*np.cross(self.realVectorA,self.realVectorB)/self.volume
        

        bv1,bv2,bv3 = self.reciprocalVectorA,self.reciprocalVectorB,self.reciprocalVectorC
        a1,a2,a3,alpha1,alpha2,alpha3= self.unitCell
        
        b1,b2,b3 = [np.linalg.norm(x) for x in [bv1,bv2,bv3]]
        beta1 = np.rad2deg(_tools.vectorAngle(bv2,bv3))
        beta2 = np.rad2deg(_tools.vectorAngle(bv3,bv1))
        beta3 = np.rad2deg(_tools.vectorAngle(bv1,bv2))
        self.cell = [a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]


    def updateSampleParameters(self,unitCell):
        """Update the sample parameters and change UB matrix as well.

        Args: 

            - unitCell (list): List of cell parameters. If None use self.unitCell

        """
        
        UB = self.UB
        B = self.B
        U = np.dot(UB,np.linalg.inv(B))
        
        self.updateCell(unitCell)
        newB = TasUBlib.calculateBMatrix(self.cell)
        newUB = np.dot(U,newB)
        
        self.UB = newUB
        self.orientationMatrix = newUB
        self.calculateProjections()




def calcProjectionVectors(R1,R2,norm=None):
    r1 = R1[:3]
    r2 = R2[:3]
    if not norm is None:
        NV = norm
    else:
        NV = np.cross(r1,r2)
    NV/= np.linalg.norm(NV)
    
    Zeros = np.isclose(NV,0.0)
    if np.sum(Zeros)==3:
        raise AttributeError('The two plane vectors are equivalen, {}, {}!'.format(r1,r2))
    
    
    if np.sum(Zeros) == 2 or np.sum(Zeros)==1: # Easy case where the two vectors are to be along the x, y, or z directions
        if Zeros[0] == True:
            V1 = np.array([1.0,0.0,0.0])
            V2 = np.cross(NV,V1)
        elif Zeros[1]:
            V1 = np.array([0.0,1.0,0.0])
            V2 = np.cross(NV,V1)
        else:
            V1 = np.array([0.0,0.0,1.0])
            V2 = np.cross(NV,V1)
    else: # The tricky case of all vectors having non-zero components.
        V1 = r1
        V2 = r2
            
    V1 = _tools.LengthOrder(V1)
    V2 = _tools.LengthOrder(V2)
    #for V in [V1,V2]: # Flip sign if needed
    #    maxArg = np.argmax(np.abs(V))
    #    V*=np.sign(V[maxArg])
    return V1,V2





