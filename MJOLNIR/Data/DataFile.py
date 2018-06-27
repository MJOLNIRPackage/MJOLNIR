import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import scipy
import matplotlib.pyplot as plt
import numpy as np
import h5py as hdf
import warnings

class DataFile(object):
    """Object to load and keep track of HdF files and their conversions"""
    def __init__(self,fileLocation):
        if fileLocation.split('.')[-1]=='nxs':
            self.type='nxs'
            with hdf.File(fileLocation) as f:
                sample=f.get('/entry/sample')
                self.sample = Sample(sample=f.get('/entry/sample'))
                self.I=np.array(f.get('entry/data/intensity'))
                self.qx=np.array(f.get('entry/data/qx'))
                self.qy=np.array(f.get('entry/data/qy'))
                self.energy=np.array(f.get('entry/data/en'))
                self.Norm=np.array(f.get('entry/data/normalization'))
                self.Monitor=np.array(f.get('entry/data/monitor'))
                instr = getInstrument(f)
                self.Ei = np.array(instr.get('monochromator/energy'))
                self.A3 = np.array(f.get('entry/sample/rotation_angle'))
                self.A4 = np.array(instr.get('detector/polar_angle'))
                self.A3Off = np.array(f.get('entry/zeros/A3'))
                if not f.get('entry/zeros/A4') is None:
                        self.A4Off = np.array(f.get('entry/zeros/A4'))
                else:
                    self.A4Off = 0.0
                self.binning = np.array(f.get('entry/reduction/MJOLNIR_algorithm_convert/binning'))[0]
                self.instrumentCalibration = np.array(f.get('entry/calibration/{}_pixels'.format(str(self.binning))))
                self.temperature = np.array(sample.get('temperature'))
                self.magneticField = np.array(sample.get('magnetic_field'))
                self.electricField = np.array(sample.get('electric_field'))

        elif fileLocation.split('.')[-1]=='h5':
            self.type='h5'
            with hdf.File(fileLocation) as f:
                sample=f.get('/entry/sample')
                self.sample = Sample(sample=f.get('/entry/sample'))
                self.Norm=np.array(f.get('entry/data/normalization'))
                self.Monitor=np.array(f.get('entry/data/monitor'))
                instr = getInstrument(f)
                self.Ei = np.array(instr.get('monochromator/energy'))
                self.A3 = np.array(f.get('entry/sample/rotation_angle'))
                self.A4 = np.array(instr.get('detector/polar_angle'))
                self.A3Off = np.array(f.get('entry/zeros/A3'))
                self.A4Off = np.array(f.get('entry/zeros/A4'))
                self.temperature = np.array(sample.get('temperature'))
                self.magneticField = np.array(sample.get('magnetic_field'))
                self.electricField = np.array(sample.get('electric_field'))
        else:
            raise AttributeError('File is not of type nxs or h5.')
        self.name = fileLocation.split('/')[-1]
        self.fileLocation = fileLocation
        self.sample.calculateProjections()


    @property
    def A3Off(self):
        return self._A3Off

    @A3Off.getter
    def A3Off(self):
        return self._A3Off

    @A3Off.setter
    def A3Off(self,A3Off):
        if A3Off is None:
            self._A3Off = 0.0
        else:
            self._A3Off = A3Off

    @property
    def A4Off(self):
        return self._A4Off

    @A4Off.getter
    def A4Off(self):
        return self._A4Off

    @A4Off.setter
    def A4Off(self,A4Off):
        if A4Off is None:
            self._A4Off = 0.0
        else:
            self._A4Off = A4Off

    def __eq__(self,other):
        return(self.fileLocation==other.fileLocation)
    
    def __str__(self):
        returnStr = 'Data file {} from the MJOLNIR software package of type {}\n'.format(self.name,self.type)
        returnStr+= 'Ei: '+ str(self.Ei) + '\nA3: ' + ','.join([str(x) for x in self.A3])
        returnStr+= '\nA4: ' + ','.join([str(x) for x in self.A4]) + '\nSample: '+str(self.sample)
        return returnStr

    def __add__(self,other):
        pass

class Sample(object):
    """Sample object to store all infortion of the sample from the experiment"""
    def __init__(self,a=None,b=None,c=None,alpha=90,beta=90,gamma=90,sample=None,name='Unknown'):
        if isinstance(sample,hdf._hl.group.Group):
            self.name = str(np.array(sample.get('name'))[0])
            self.orientationMatrix = np.array(sample.get('orientation_matrix'))
            self.planeNormal = np.array(sample.get('plane_normal'))
            self.polarAngle = np.array(sample.get('polar_angle'))
            self.rotationAngle = np.array(sample.get('rotation_angle'))
            self.unitCell = np.array(sample.get('unit_cell'))
            
        elif np.all([a is not None,b is not None, c is not None]):
            self.unitCell = np.array([a,b,c,alpha,beta,gamma])
            self.orientationMatrix = np.array([[1,0,0],[0,1,0]])
            self.planeNormal = np.array([0,0,1])
            self.polarAngle = np.array(0)
            self.rotationAngle = np.array(0)
            self.name=name
        else:
            print(sample)
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
        return np.all([self.name==other.name,np.all(self.unitCell==other.unitCell),\
        np.all(self.orientationMatrix==other.orientationMatrix)])

    def calculateProjections(self):
        """Calculate projections and generate projection angles."""
        try:
            self.unitCell
            self.orientationMatrix
        except:
            raise AttributeError('No unit cell is set for sample object.')
        self.realVectorA = np.array([self.a,0,0])
        self.realVectorB = np.dot(-np.array([self.b,0,0]),rotationMatrix(0,0,self.gamma))
        self.realVectorC = np.dot(np.array([self.c,0,0]),rotationMatrix(0,self.beta,0))
        
        self.volume = np.abs(np.dot(self.realVectorA,np.cross(self.realVectorB,self.realVectorC)))
        self.reciprocalVectorA = 2*np.pi*np.cross(self.realVectorB,self.realVectorC)/self.volume
        self.reciprocalVectorB = 2*np.pi*np.cross(self.realVectorC,self.realVectorA)/self.volume
        self.reciprocalVectorC = 2*np.pi*np.cross(self.realVectorA,self.realVectorB)/self.volume
        ## Ensure that aStar is along the x-axis
        RotMatrix = rotate2X(self.reciprocalVectorA)
        #angle = vectorAngle(self.reciprocalVectorA,np.array([1,0,0])) # TODO: make general!!!
        self.reciprocalVectorA=np.dot(self.reciprocalVectorA,RotMatrix)
        self.reciprocalVectorB=np.dot(self.reciprocalVectorB,RotMatrix)
        self.reciprocalVectorC=np.dot(self.reciprocalVectorC,RotMatrix)

        reciprocalMatrix = np.array([self.reciprocalVectorA,self.reciprocalVectorB,self.reciprocalVectorC])
        self.projectionVector1 = np.array(np.dot(self.orientationMatrix[0],reciprocalMatrix))
        self.projectionVector2 = np.array(np.dot(self.orientationMatrix[1],reciprocalMatrix))
        self.projectionAngle = vectorAngle(self.projectionVector1,self.projectionVector2)

        if np.isclose(0,self.projectionAngle):
            raise AttributeError("The provided orientations are equal.")

        self.projectionMatrix = np.array([\
        [np.linalg.norm(self.projectionVector1),np.cos(self.projectionAngle)*np.linalg.norm(self.projectionVector2)]\
        ,[0,np.sin(self.projectionAngle)*np.linalg.norm(self.projectionVector2)]])

    def tr(self,h,k):
        """Convert from curved coordinate to rectlinear coordinate."""
        h, k = np.asarray(h), np.asarray(k)
        Pos = np.dot(self.projectionMatrix,[h,k])
        return Pos[0],Pos[1]
        
    def inv_tr(self, x,y):
        """Convert from rectlinear coordinate to curved coordinate."""
        x, y = np.asarray(x), np.asarray(y)
        Pos = np.dot(np.linalg.inv(self.projectionMatrix),[x,y])
        return Pos[0],Pos[1]   


    def format_coord(self,x,y):
        pos = self.inv_tr(x,y)
        rlu = self.orientationMatrix[0]*pos[0]+self.orientationMatrix[1]*pos[1]
        return "h = {0:.3f}, k = {1:.3f}, l = {2:.3f}".format(rlu[0],rlu[1],rlu[2])

    def __str__(self):
        returnStr = 'Sample ' + self.name + '\n'
        #if not self.temperature is None: returnStr+= 'Temperatur: '+str(self.temperature)+'\n'
        #if not self.magneticField is None: returnStr+= 'Magnetic Field: '+str(self.magneticField)+'\n'
        #if not self.electricField is None: returnStr+= 'Electric Field: '+str(self.electricField)+'\n'
        returnStr+= 'Unit cell: \n' + str(self.unitCell) + '\n'
        returnStr+= 'Orientation matrix: \n' + str(self.orientationMatrix) +'\n'

        return returnStr


def rotationMatrix(alpha,beta,gamma,format='deg'):
    if format=='deg':
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    Rx = np.array([[1,0,0],[0,np.cos(alpha),-np.sin(alpha)],[0,np.sin(alpha),np.cos(alpha)]])
    Ry = np.array([[np.cos(beta),0,np.sin(beta)],[0,1,0],[-np.sin(beta),0,np.cos(beta)]])
    Rz = np.array([[np.cos(gamma),-np.sin(gamma),0],[np.sin(gamma),np.cos(gamma),0],[0,0,1]])
    return np.dot(Rz,np.dot(Ry,Rx))

def vectorAngle(V1,V2):
    return np.arccos(np.dot(V1,V2.T)/(np.linalg.norm(V1)*np.linalg.norm(V2)))

def getNX_class(x,y,attribute):
    try:
        variableType = y.attrs['NX_class']
    except:
        variableType = ''
    if variableType==attribute:
        return x
    
def getInstrument(file):
    location = file.visititems(lambda x,y: getNX_class(x,y,b'NXinstrument'))
    return file.get(location)

def extractData(files):
    if not isinstance(files,list):
        files = [files]
    I = []
    posx = []
    posy = []
    energy = []
    Norm = []
    Monitor = []
    a3 = []
    a4 = []
    Ei = []

    for datafile in files:
        I.append(datafile.I)
        posx.append(datafile.qx)
        posy.append(datafile.qy)
        energy.append(datafile.energy)
        Norm.append(datafile.Norm)
        Monitor.append(datafile.Monitor)
        if np.array(datafile.A3Off).shape is ():
            datafile.A3Off = 0.0
        a3.append(datafile.A3-datafile.A3Off)

        if np.array(datafile.A4Off).shape is ():
            datafile.A4Off = 0.0
        a4.append(datafile.A4-datafile.A4Off)
        Ei.append(datafile.Ei)
        
    
    I = np.concatenate(I)
    qx = np.concatenate(posx)
    qy = np.concatenate(posy)
    energy = np.concatenate(energy)
    Norm = np.concatenate(Norm)
    Monitor = np.concatenate(Monitor)
    a3 = np.concatenate(a3)
    a4 = np.concatenate(a4)
    Ei = np.concatenate(Ei)

    return I,qx,qy,energy,Norm,Monitor,a3,a4,files[-1].instrumentCalibration,Ei # FIXME: Might be a problem if combining dataset with different calibrations!

def rotMatrix(v,theta): # https://en.wikipedia.org/wiki/Rotation_matrix
    v/=np.linalg.norm(v)
    m11 = np.cos(theta)+v[0]**2*(1-np.cos(theta))
    m12 = v[0]*v[1]*(1-np.cos(theta))-v[2]*np.sin(theta)
    m13 = v[0]*v[2]*(1-np.cos(theta))+v[1]*np.sin(theta)
    m21 = v[0]*v[1]*(1-np.cos(theta))+v[2]*np.sin(theta)
    m22 = np.cos(theta)+v[1]**2*(1-np.cos(theta))
    m23 = v[1]*v[2]*(1-np.cos(theta))-v[0]*np.sin(theta)
    m31 = v[0]*v[2]*(1-np.cos(theta))-v[1]*np.sin(theta)
    m32 = v[1]*v[2]*(1-np.cos(theta))+v[0]*np.sin(theta)
    m33 = np.cos(theta)+v[2]**2*(1-np.cos(theta))
    return np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])


def rotate2X(v):
    if np.isclose(v[2]/np.linalg.norm(v),1): # v is along z
        return rotMatrix([0,1,0],np.pi/2)
    # Find axis perp to v and proj v into x-y plane -> rotate 2 plane and then to x
    vRotInPlane = np.array([-v[1],v[0],0])
    vPlan = np.array([v[0],v[1],0])
    theta = np.arccos(np.dot(v,vPlan)/(np.linalg.norm(v)*np.linalg.norm(vPlan)))
    R = rotMatrix(vRotInPlane,theta)
    v2 = np.dot(R,v)
    theta2 = np.arccos(np.dot(v2,np.array([1,0,0]))/np.linalg.norm(v2))
    R2 = rotMatrix(np.array([0,0,1.0]),-theta2)
    
    Rotation = np.dot(R2,R)
    return Rotation
# --------------------------- TESTS -------------------------

def test_sample_exceptions():
    try: # No parameters given
        s1 = Sample()
        assert False
    except:
        assert True

    try: # negative parameters given
        s1 = Sample(a=-1,b=1,c=1)
        assert False
    except:
        assert True

    try: # negative parameters given
        s1 = Sample(a=1,b=-1,c=1)
        assert False
    except:
        assert True
    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=-1)
        assert False
    except:
        assert True

    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=1,alpha=200)
        assert False
    except:
        assert True
    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=1,beta=-10)
        assert False
    except:
        assert True
    try: # negative parameters given
        s1 = Sample(a=1,b=1,c=1,gamma=-10)
        assert False
    except:
        assert True

def test_vectorAngle():
    v1 = np.array([1,0,0])
    v2 = np.array([0,1,0])
    v3 = np.array([1,1,1])
    theta1 = vectorAngle(v1,v2)
    theta2 = vectorAngle(v1,v3)
    theta3 = vectorAngle(v2,v3)
    print(theta2)
    assert(np.isclose(theta1,np.pi/2))
    assert(np.isclose(theta2,theta3))
    assert(np.isclose(theta2,0.955316618125))
    
def test_Rotation_matrix():
    M1 = rotationMatrix(0,0,0)

    assert(np.all(np.isclose(M1,np.identity(3))))

    M2 = rotationMatrix(90,0,0)
    M3 = rotationMatrix(np.pi/2,np.pi/3,np.pi/4,format='rad')
    M4 = rotationMatrix(90,60,45)
    
    assert(np.all(np.isclose(M2,np.array([[1,0,0],[0,0,-1],[0,1,0]]))))
    assert(np.all(np.isclose(M3,M4)))
    
    M4_check = np.array([[3.53553391e-01,6.12372436e-01,7.07106781e-01],[3.53553391e-01,6.12372436e-01,-7.07106781e-01],[-8.66025404e-01,5.00000000e-01,3.06161700e-17]])
    assert(np.all(np.isclose(M4,M4_check)))

def test_equality():
    s1 = Sample(1,2,3,90,90,120)
    s2 = Sample(1,2,3,90,90,120)
    assert(s1==s2)

def test_calculateProjections():

    s1 = Sample(np.pi*2,np.pi*2,np.pi*2,90,90,120)

    s1.orientationMatrix = np.array([[1,0,0],[0,1,0]])
    s1.calculateProjections()

    theta = 2*np.pi/3 # 120 degrees between projection vectors
    AStar = np.linalg.norm(s1.reciprocalVectorA)
    BStar = np.linalg.norm(s1.reciprocalVectorB)
    CStar = np.linalg.norm(s1.reciprocalVectorC)
    ## Projection is given by [[a*,cos(theta)b*],[0,sin(tetha)b*]]
    assert(np.all(np.isclose(s1.projectionMatrix,np.array([[AStar*1,BStar*np.cos(theta)],[0,BStar*np.sin(theta)]]))))
    assert(np.all(np.isclose(s1.tr(1,1),np.array([AStar+np.cos(theta)*BStar,np.sin(theta)*BStar]))))

    point = s1.tr(3.5,7.2)
    reverse = s1.inv_tr(point[0],point[1])
    assert(np.all(np.isclose(reverse,np.array([3.5,7.2]))))

    string = s1.format_coord(point[0],point[1])
    assert(string=='h = 3.500, k = 7.200, l = 0.000')


def test_DataFile():
    try:
        DF = DataFile('/nope.txt')
        assert False
    except:
        assert True
    files = ['TestData/cameasim2018n000001.h5',
             'TestData/cameasim2018n000001.nxs']
    DF1 = DataFile(files[0])
    DF2 = DataFile(files[1])
            
    assert(DF1.sample == DF2.sample)

def test_DataFile_rotations():
    vectors = [np.array([0,0,3.0]),np.array([1.0,0.0,0.0]),np.array([0.0,1.0,0.0]),np.random.rand(3)]
    rotations = [rotate2X(v) for v in vectors]
    rotVector = [np.dot(rotations[i],vectors[i]) for i in range(len(vectors))]
    for i in range(len(rotVector)):
        assert(np.isclose(np.linalg.norm(vectors[i]),np.linalg.norm(rotVector[i])))
        assert(np.isclose(rotVector[i][0],np.linalg.norm(rotVector[i])))
