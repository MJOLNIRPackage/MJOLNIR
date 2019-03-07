import numpy as np

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
              [0,0,b3]]) # 2*np.pi/a3
    return B

def matFromTwoVectors(v1,v2):
    a1 = v1/np.linalg.norm(v1)
    a2 = v2/np.linalg.norm(v2)
    a3 = np.cross(a1,a2)
    a3/=np.linalg.norm(a3)
    a2 = np.cross(a1,a3)
    return np.array([a1,a2,a3])


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

def test_TasUBDeg(): # Test that the two libraries are equivalent in calculating the UB matrices
                
    
    r1 = [ 1.0000000e+00,  0.0000000e+00,  0.0000000e+00,  8.4377825e-06,
       -4.4549999e+01, -1.4062500e-05,  3.3378601e-06,  5.0000076e+00,
        4.9678469e+00]
    r2 = [ 0.0000000e+00, -1.0000000e+00,  0.0000000e+00, -6.0000084e+01,
       -4.4549999e+01, -1.4062500e-05,  3.3378601e-06,  5.0000076e+00,
        4.9678469e+00]
        # H, K, L, A3, A4, sgu, sgl, Ei, Ef

    cellDeg = [6.11,   6.11,  11.35, 1.187430040454027, 1.1874300210500532, 0.5535845899562842, 90.  ,  90.  , 120., 90., 90., 60.]
    #[a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]


    UB = calcTasUBFromTwoReflections(cellDeg, r1, r2)

    UBFile = np.array([[ 1.7459887e-01,  2.4665387e-02,  2.1624945e-08],
       [-7.2323568e-02, -1.8736884e-01, -5.1334577e-09],
       [ 4.7066340e-08,  1.6969278e-08, -8.8105723e-02]])
    
    test = np.all(np.isclose(UB,2*np.pi*UBFile,atol=1e-6))
    
    print(cellDeg)
    print('--------------')
    print(UB)
    print('--------------')
    print(UBFile)
    print('--------------')
    assert(test)


    UBInv = np.linalg.inv(UB)
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
        
    assert(np.all(np.isclose([QH[0]],[qh,qk,ql],atol=1e-4)))

def test_TasUBDeg_CreateUB():

    UBFile = np.array([[ 1.7459887e-01,  2.4665387e-02,  2.1624945e-08],
       [-7.2323568e-02, -1.8736884e-01, -5.1334577e-09],
       [ 4.7066340e-08,  1.6969278e-08, -8.8105723e-02]])*2*np.pi

    cell = [6.11,   6.11,  11.35, 1.187430040454027, 1.1874300210500532, 0.5535845899562842, 90.  ,  90.  , 120., 90., 90., 120.]
    B = calculateBMatrix(cell)
    OM = 360-22.500 # Offset known 
    sgu = 0.0
    sgl = 0.0
    UB = calcUBFromAngles(B,OM,sgu,sgl)
    
    assert(np.all(np.isclose(UBFile,UB,atol=4)))

def test_TasUBDEG_CalculateAngles(): # TODO: Redo these calculations as one needs to do the trick with a3 to calculate anything correctly...
    cell = [6.11,   6.11,  11.35, 1.187430040454027, 1.1874300210500532, 0.5535845899562842, 90.  ,  90.  , 120., 90., 90., 60.]
    B = calculateBMatrix(cell)
    OM = 0.00 # Offset known 
    sgu = 0.0
    sgl = 0.0
    UB = calcUBFromAngles(B,OM,sgu,sgl)

    planeNormal = np.array([0,0,1.0])
    qe = np.array([1.0,0.0,0.0,5.0,5])
    ss = 1 # Scattering sense
    A3Off = 0.0
    
    UBINV = np.linalg.inv(UB)
    QE = np.array([[1.0,0.0,0.0,5.0,5],
                   [0.0,0.5,0.0,8.0,9.2],
                   [-1.1,-0.1,0.0,5.4,4.4]])
    for qe in QE:

        a3,a4,sgu,sgl = calcTasQAngles(UB,planeNormal,ss,A3Off,qe)

        print('------------------\nA3:{}\nA4{}'.format(a3,a4))
        hkl = calcTasQH(UBINV,[a3,a4],qe[3],qe[4])[0]
        print('{}'.format(hkl))
        assert(np.all(np.isclose([sgu,sgl],0.0))) # Sgu and sgl are 0 by definition
        assert(np.all(np.isclose(hkl,qe[:3])))
    
