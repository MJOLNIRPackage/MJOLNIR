import numpy as np

factorsqrtEK = 0.694692

def uFromAngles(om,sgu,sgl):
    return np.array([np.cos(om)*np.cos(sgl),-np.sin(om)*np.cos(sgu)+np.cos(om)*np.sin(sgl)*np.sin(sgu),
                      np.sin(om)*np.sin(sgu)+np.cos(om)*np.sin(sgl)*np.cos(sgu)]).T
    
    
def calcUBFromAngles(B, om, sgu, sgl):
    N = np.array([[1.0,0,0],[0,np.cos(sgu),-np.sin(sgu)],[0,np.sin(sgu),np.cos(sgu)]])
    M = np.array([[np.cos(sgl),0,np.sin(sgl)],[0,1,0],[-np.sin(sgl),0,np.cos(sgl)]])    
    OM= np.array([[np.cos(om),-np.sin(om),0],[np.sin(om),np.cos(om),0],[0,0,1]])
    
    UB = np.dot(OM,M)
    UB = np.dot(UB,N)
    UB = np.dot(UB,B)

    return UB  

def calculateBMatrix(cell):
    [a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3] = cell
    
    B = np.array([[b1,b2*np.cos(beta3),b3*np.cos(beta2)],
              [0,b2*np.sin(beta3),-b3*np.sin(beta2)*np.cos(alpha1)],
              [0,0,1/a3]])
    return B

def matFromTwoVectors(v1,v2):
    a1 = v1/np.linalg.norm(v1)
    a2 = v2/np.linalg.norm(v2)
    a3 = np.cross(a1,a2)
    a3/=np.linalg.norm(a3)
    a2 = np.cross(a1,a3)
    return np.array([a1,a2,a3])
#
#def rtan(y, x):
#
#    val = 0.0
#    if ((x == 0.) and (y == 0.)):
#        return .0
#      
#    if (x == 0.):
#        if (y < 0.):
#            return -np.pi / 2.;
#        else:
#          return np.pi / 2. 
#    if(np.abs(y) < np.abs(x)):
#        val = np.arctan(np.abs(y / x))
#        if (x < 0.):
#            val = np.pi - val;
#        
#        if (y < 0.):
#            val = -val
#            
#        return val
#    else:
#        val = np.pi / 2. - np.arctan(np.abs(x / y))
#        if (x < 0.):
#            val = np.pi - val;
#        if (y < 0.):
#            val = -val
#        return val;


def calcTheta(ki,kf,stt):
    return np.arctan2(np.abs(ki) - np.abs(kf) * np.cos(stt), #rtan
              np.abs(kf) * np.sin(stt));
    
def calcTasUVectorFromAngles(r):
    A3 = np.deg2rad(r[3])
    A4 = np.deg2rad(r[4])
    ss = np.sign(A4) # Scattering sense
    stt = np.abs(A4) # Sample 2theta
    Ei = r[7]
    Ef = r[8]
    ki = np.sqrt(Ei)*factorsqrtEK
    kf = np.sqrt(Ef)*factorsqrtEK
    sgu = np.deg2rad(r[5])
    sgl = np.deg2rad(r[6])
    
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
    return np.dot(UB,qe    )

def buildTVMatrix(U1V,U2V):
    U1V /=np.linalg.norm(U1V)
    U2V /=np.linalg.norm(U2V)
    T3V = np.cross(U1V.T, U2V.T).T
    
    return np.array([U1V,U2V,T3V])

def buildRMatrix(UB, planeNormal, qe):
    U1V = tasReflectionToQC(qe,UB)
    U1V/=np.linalg.norm(U1V)
    U2V = np.cross(planeNormal.T, U1V.T).T
    if np.linalg.norm(U2V)<0.001:
        raise AttributeError('Calculate length of U2V too small ({})'.format(U2V))
    
    TV = buildTVMatrix(U1V, U2V)
    TVinv = np.linalg.inv(TV)
    return TVinv

def calcTasQAngles(UB,planeNormal,ss,A3Off,qe):
    R = buildRMatrix(UB,planeNormal,qe)
    cossgl = np.sqrt(R[0,0]*R[0,0]+R[1,0]*R[1,0])
    sgl = ss*np.arctan2(-R[2,0],cossgl)
    om = np.arctan2(R[1,0]/cossgl, R[0,0]/cossgl)
    sgu = np.arctan2(R[2,1]/cossgl, R[2,2]/cossgl)
    
    QC = tasReflectionToQC(qe, UB)
    q = np.linalg.norm(QC)*2.0*np.pi;
    
    Ei = qe[3]
    Ef = qe[4]
    
    ki = np.sqrt(Ei)*factorsqrtEK
    kf = np.sqrt(Ef)*factorsqrtEK
    
    cos2t =(ki**2 + kf**2 - q**2) / (2. * np.abs(ki) * np.abs(kf))
    
    A4 = np.arccos(cos2t)
    theta = calcTheta(ki, kf, A4)
    A3 = om + ss*theta + A3Off
    A3 = np.mod(A3 + ss*np.pi,2*np.pi) - ss*np.pi
    
    return A3,A4,sgu,sgl
    

def calcTasMisalignment(UB,planeNormal,qe):
    R = buildRMatrix(UB,planeNormal, qe[:3])
    om = np.arctan2(R[1,0], R[0,0])
    return om


def calcTasQH(UBINV,angles,Ei,Ef):

  A3,A4 = angles
  r = [0,0,0,A3,A4,0.0,0.0,Ei,Ef]
  
  ki = np.sqrt(Ei)*factorsqrtEK
  kf = np.sqrt(Ef)*factorsqrtEK
  
  QV = calcTasUVectorFromAngles(r);
  
  q = np.sqrt(ki**2 +kf**2-
           2. *ki *kf * np.cos(np.deg2rad(A4)));
  q /= 2. * np.pi
  
  q = np.expand_dims(np.swapaxes(q,0,2),3)
  QV *=q
  
  Q = np.einsum('ij,klmj->klmi',UBINV,QV)
  QxQy = np.swapaxes(QV*np.pi*2,0,2)
  return Q,QxQy[:,:,:,0],QxQy[:,:,:,1]