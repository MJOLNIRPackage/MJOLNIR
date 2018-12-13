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
  
  if len(QV.shape)>2:
    q = np.expand_dims(np.swapaxes(q,0,2),3)
  QV *=q

  Q = np.einsum('ij,...j->...i',UBINV,QV)
  if len(QV.shape)>2:
    QxQy = np.swapaxes(QV*np.pi*2,0,2)
    return Q,QxQy[:,:,:,0],QxQy[:,:,:,1]
  else:
    QxQy = QV*np.pi*2
    return Q,QxQy
  
  



def test_TasUBDeg(): # Test that the two libraries are equivalent in calculating the UB matrices
    from . import TasUBlibDEG
        #[a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]
    
    r1 = [ 1.0000000e+00,  0.0000000e+00,  0.0000000e+00,  8.4377825e-06,
       -4.4549999e+01, -1.4062500e-05,  3.3378601e-06,  5.0000076e+00,
        4.9678469e+00]
    r2 = [ 0.0000000e+00, -1.0000000e+00,  0.0000000e+00, -6.0000084e+01,
       -4.4549999e+01, -1.4062500e-05,  3.3378601e-06,  5.0000076e+00,
        4.9678469e+00]

    cellDeg = [6.11,   6.11,  11.35, 1.187430040454027, 1.1874300210500532, 0.5535845899562842, 90.  ,  90.  , 120., 90., 90., 60.]
    cell = [6.11,   6.11,  11.35, 1.187430040454027/(2*np.pi), 1.1874300210500532/(2*np.pi), 0.5535845899562842/(2*np.pi), 90.  ,  90.  , 120., 90., 90., 60.]

    cell[6:] = np.deg2rad(cell[6:])
    UB = calcTasUBFromTwoReflections(cell, r1, r2)
    UBDeg = TasUBlibDEG.calcTasUBFromTwoReflections(cellDeg, r1, r2)

    UBFile = np.array([[ 1.7459887e-01,  2.4665387e-02,  2.1624945e-08],
       [-7.2323568e-02, -1.8736884e-01, -5.1334577e-09],
       [ 4.7066340e-08,  1.6969278e-08, -8.8105723e-02]])
    
    test1 = np.all(np.isclose(UB,UBFile,atol=1e-6))
    test2 = np.all(np.isclose(UBDeg,2*np.pi*UBFile,atol=1e-6))
    if not test1 or not test2:
        print(cell)
        print('--------------')
        print(UB)
        print('--------------')
        print(UBFile)
        print('--------------')
        print(UBDeg/(2*np.pi))
        assert(test1)
        assert(test2)


    # All UB matrices are equal up to 2 pi for Deg
    UBInv = np.linalg.inv(UB)
    UBDegInv = np.linalg.inv(UBDeg)
    A3,A4 = [-23,-55]
    angles = np.array([A3,A4])
    Ei = 5.05
    Ef = 3.2
    ## Values from SIX
    qh =  1.174406
    qk =  -0.127853
    ql =  0.000000
    qm =  1.325155

    QH = calcTasQH(UBInv,angles,Ei,Ef)
    QHDeg = TasUBlibDEG.calcTasQH(UBDegInv,np.array([A3,A4]),Ei,Ef)
    
    assert(np.all(np.isclose([QH[0]],[qh,qk,ql],atol=1e-4)))
    assert(np.all(np.isclose([QHDeg[0]],[qh,qk,ql],atol=1e-4)))

    assert(np.all(np.isclose(QH,QHDeg,atol=1e-8)))