from MJOLNIR.TasUBlibDEG import calcTasUBFromTwoReflections, calcTasQH, calculateBMatrix, calcUBFromAngles, calcTasQAngles, calcCell
import numpy as np

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
    

def test_TasUBDeg_calculateCell():
    cellParams = np.array([5.103,5.103,13.755,90.0,90.0,120.0])

    cell = calcCell(cellParams)
    result = np.array([5.103, 5.103, 13.755, 1.4217514122941153, 1.4217514122941153, 0.4567928249494428, 90.0, 90.0, 120.0, 90., 90., 60.])
    assert(np.all(np.isclose(cell,result)))

    cellParams = np.array([5.103,6.103,13.755,80.0,93.0,120.0])

    cell = calcCell(cellParams)
    result = np.array([5.103, 6.103, 13.755, 1.4229152753265009, 1.2064634997894217, 0.4642192953135897, 80.0, 93.0, 120.0, 99.81858754108276, 92.31754754830274, 60.05495280628465])
    assert(np.all(np.isclose(cell,result)))