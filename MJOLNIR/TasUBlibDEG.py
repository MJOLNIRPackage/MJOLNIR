import numpy as np
from MJOLNIR import _tools

def cosd(x):
    return np.cos(np.deg2rad(x))

def sind(x):
    return np.sin(np.deg2rad(x))

def tand(x):
    return np.tan(np.deg2rad(x))

def arccosd(x):
    return np.rad2deg(np.arccos(x))

def arcsind(x):
    return np.rad2deg(np.arcsin(x))

def arctan2d(y,x):
    return np.rad2deg(np.arctan2(y,x))

factorsqrtEK = 0.694692

def uFromAngles(om,sgu,sgl):
    return np.array([cosd(om)*cosd(sgl),-sind(om)*cosd(sgu)+cosd(om)*sind(sgl)*sind(sgu),
                      sind(om)*sind(sgu)+cosd(om)*sind(sgl)*cosd(sgu)]).T
    
    
def calcUBFromAngles(B, om, sgu, sgl):
    N = np.array([[1.0,0,0],[0,cosd(sgu),-sind(sgu)],[0,sind(sgu),cosd(sgu)]])
    M = np.array([[cosd(sgl),0,sind(sgl)],[0,1,0],[-sind(sgl),0,cosd(sgl)]])    
    OM= np.array([[cosd(om),-sind(om),0],[sind(om),cosd(om),0],[0,0,1]])
    
    UB = np.dot(OM,M)
    UB = np.dot(UB,N)
    UB = np.dot(UB,B)

    return UB  

def calculateBMatrix(cell):
    [a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3] = cell
    
    B = np.array([[b1,b2*cosd(beta3),b3*cosd(beta2)],
              [0,b2*sind(beta3),-b3*sind(beta2)*cosd(alpha1)],
              [0,0,2*np.pi/a3]]) # 2*np.pi/a3
    return B

def matFromTwoVectors(v1,v2):
    a1 = v1/np.linalg.norm(v1)
    a2 = v2/np.linalg.norm(v2)
    a3 = np.cross(a1,a2)
    a3/=np.linalg.norm(a3)
    a2 = np.cross(a1,a3)
    return np.array([a1,a2,a3]).T


def calcTheta(ki,kf,stt):
    return arctan2d(np.abs(ki) - np.abs(kf) * cosd(stt), #rtan
              np.abs(kf) * sind(stt));
    
def calcTasUVectorFromAngles(r):
    A3 = r[3]
    A4 = r[4]
    ss = np.sign(A4) # Scattering sense
    stt = np.abs(A4) # Sample 2theta
    Ei = r[7]
    Ef = r[8]
    ki = np.sqrt(Ei)*factorsqrtEK
    kf = np.sqrt(Ef)*factorsqrtEK
    sgu = r[5]
    sgl = r[6]
    
    theta = calcTheta(ki, kf, stt)
    try:
        om = A3-ss*theta
    except:
        om = A3.reshape(-1,1,1)-ss*theta
    return uFromAngles(om,sgu,sgl)
  
def calcTasUBFromTwoReflections(cell, r1, r2):
    B = calculateBMatrix(cell)
    R1 = np.array(r1[:3])
    R2 = np.array(r2[:3])
    h1 = np.dot(B,R1)
    h2 = np.dot(B,R2)
    
    HT = matFromTwoVectors(h1, h2)
    u1 = calcTasUVectorFromAngles(r1)
    u2 = calcTasUVectorFromAngles(r2)
    
    UT = matFromTwoVectors(u1, u2)
    HTT = HT.T
    
    U = np.dot(UT,HTT)
    UB = np.dot(U,B)
    return UB
    
def tasReflectionToQC(qe,UB):
    qe = qe[:3]
    return np.dot(UB,qe)

def buildTVMatrix(U1V,U2V):
    U1V /=np.linalg.norm(U1V)
    U2V /=np.linalg.norm(U2V)
    T3V = np.cross(U1V.T, U2V.T).T
    
    return np.array([U1V.flatten(),U2V.flatten(),T3V.flatten()])

def buildRMatrix(UB, planeNormal, qe):
    U1V = tasReflectionToQC(qe,UB)
    U1V/=np.linalg.norm(U1V)
    U2V = np.cross(np.dot(UB,planeNormal).T, U1V.T).T
    if np.linalg.norm(U2V)<0.001:
        raise AttributeError('Calculate length of U2V too small ({})'.format(U2V))
    
    TV = buildTVMatrix(U1V, U2V)
    TVinv = np.linalg.inv(TV)
    return TVinv

def calcTasQAngles(UB,planeNormal,ss,A3Off,qe):
    R = buildRMatrix(UB,planeNormal,qe)
    cossgl = np.sqrt(R[0,0]*R[0,0]+R[1,0]*R[1,0])
    sgl = ss*arctan2d(-R[2,0],cossgl)
    om = arctan2d(R[1,0]/cossgl, R[0,0]/cossgl)
    sgu = arctan2d(R[2,1]/cossgl, R[2,2]/cossgl)
    
    QC = tasReflectionToQC(qe, UB)
    q = np.linalg.norm(QC);
    
    Ei = qe[3]
    Ef = qe[4]
    
    ki = np.sqrt(Ei)*factorsqrtEK
    kf = np.sqrt(Ef)*factorsqrtEK
    
    cos2t =(ki**2 + kf**2 - q**2) / (2. * np.abs(ki) * np.abs(kf))
    
    A4 = arccosd(cos2t)
    theta = calcTheta(ki, kf, A4)
    A3 = om + ss*theta + A3Off
    A3 = np.mod(A3 + ss*180.0,360.0) - ss*180.0
    
    return -A3,-A4,sgu,sgl
    

def calcTasMisalignment(UB,planeNormal,qe):
    R = buildRMatrix(UB,planeNormal, qe[:3])
    om = arctan2d(R[1,0], R[0,0])
    return om


def calcTasQH(UBINV,angles,Ei,Ef,A3Off=0):

  A3,A4 = angles
  r = [0,0,0,A3+A3Off,A4,0.0,0.0,Ei,Ef]
  
  ki = np.sqrt(Ei)*factorsqrtEK
  kf = np.sqrt(Ef)*factorsqrtEK
  
  QV = calcTasUVectorFromAngles(r);
  
  q = np.sqrt(ki**2 +kf**2-
           2. *ki *kf * cosd(A4));
  
  if len(QV.shape)>2: # If multidimensional QV
    q = np.expand_dims(np.swapaxes(q,0,2),3)

  QV *=q
  
  Q = np.einsum('ij,...j->...i',UBINV,QV)
  
  if len(QV.shape)>2:
    QxQy = np.swapaxes(QV,0,2)
    return Q,QxQy[:,:,:,0],QxQy[:,:,:,1]

  return Q,QV

def addAuxReflection(cell,r1,r2,ss=1.0):
    h,k,l = r2[:3]#np.array([0,0,1])
    
    Ei,Ef = r1[-2:]
    
    ki,kf = np.sqrt([Ei,Ef])*factorsqrtEK
    
    B = calculateBMatrix(cell)
    
    #qe = np.array([h,k,l,Ei,Ef])
    
    
    stt1 = r1[4]
    ### makeAuxReflection
    theta = calcTheta(ki,kf,stt1)
    om = r1[3]-ss*theta
    

    om+= tasAngleBetweenReflections(B,r1[:3],[h,k,l])
    
    
    QC = np.dot(B,[h,k,l])
    q = np.linalg.norm(QC)
    
    cos2t = (ki * ki + kf * kf - q * q) /\
                (2. * np.abs(ki) * np.abs(kf))
    r2stt = np.rad2deg(np.arccos(cos2t))
    
    r2A3 = om+ss*theta
    r2A3 = np.mod(r2A3+ss*180.0,360.0)-ss*180.0
    
    r2 = np.array([h,k,l,r2A3,r2stt,0.0,0.0,Ei,Ef])

    return r2

def calTwoTheta(B,qe,ss):

  h,k,l,Ei,Ef = qe
  ki,kf = np.sqrt([Ei,Ef])*factorsqrtEK
  QC = np.dot(B,qe[:3])
  
  q = np.linalg.norm(QC);
 
  cos2t = (ki * ki + kf * kf - q * q) /\
              (2. * np.abs(ki) * np.abs(kf))
      
  return ss * np.rad2deg(np.arccos(cos2t))
 
def tasAngleBetweenReflections(B,r1,r2):
  ## Angle between the two reflections
  H1 = r1[:3]
  H2 = r2[:3]
  chi1 = np.dot(B,H1)
  chi2 = np.dot(B,H2)
  angle = np.rad2deg(np.arccos(np.dot(chi1,chi2)/(np.linalg.norm(chi1)*np.linalg.norm(chi2))))
  return angle


def calcCell(cell):
    """Calculate reciprocal lattice vectors from real space cell parameters"""

    a1,a2,a3,alpha1,alpha2,alpha3 = cell
    realVectorA = np.array([a1,0,0])
    realVectorB = a2*np.array([cosd(alpha3),sind(alpha3),0.0])#np.dot(np.array([b,0,0]),rotationMatrix(0,0,gamma))
    realVectorC = a3*np.array([cosd(alpha2),(cosd(alpha1)-cosd(alpha2)*cosd(alpha3))/sind(alpha3),
    np.sqrt(1-cosd(alpha2)**2-((cosd(alpha1)-cosd(alpha2)*cosd(alpha3))/sind(alpha3))**2)])#np.dot(np.array([c,0,0]),rotationMatrix(0,beta,0))

    volume = np.abs(np.dot(realVectorA,np.cross(realVectorB,realVectorC)))
    reciprocalVectorA = 2*np.pi*np.cross(realVectorB,realVectorC)/volume
    reciprocalVectorB = 2*np.pi*np.cross(realVectorC,realVectorA)/volume
    reciprocalVectorC = 2*np.pi*np.cross(realVectorA,realVectorB)/volume


    bv1,bv2,bv3 = reciprocalVectorA,reciprocalVectorB,reciprocalVectorC
    

    b1,b2,b3 = [np.linalg.norm(x) for x in [bv1,bv2,bv3]]
    beta1 = np.rad2deg(_tools.vectorAngle(bv2,bv3))
    beta2 = np.rad2deg(_tools.vectorAngle(bv3,bv1))
    beta3 = np.rad2deg(_tools.vectorAngle(bv1,bv2))
    cell = [a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]
    return cell