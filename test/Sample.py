from MJOLNIR import TasUBlibDEG
from MJOLNIR.Data.Sample import Sample
import MJOLNIR.Data.DataFile
from MJOLNIR import _tools
import os
import numpy as np
import MJOLNIR.TasUBlibDEG as TasUBlib

dataPath = 'Data'



def test_sample_exceptions():
    try: # No parameters given
        s1 = Sample(a=None)
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

def test_Sample_conversions():
    df = MJOLNIR.Data.DataFile.DataFile(os.path.join(dataPath,'camea2018n000178.hdf'))
    sample = df.sample

    # Check that tr and inv_tr results in the same
    p0,p1 = 2.25,-1.23
    Qx,Qy = sample.tr(p0,p1)
    P0,P1 = sample.inv_tr(Qx,Qy)

    assert(np.all(np.isclose([p0,p1],[P0,P1])))

    # Check that format_coord prints something
    string = sample.format_coord(Qx,Qy) # format is h = ??, k = ??, l = ??
    hkl = [float(x.split('=')[-1]) for x in string.split(',')]
    print(string)
    print(hkl)

    QxQyFromSample = sample.calculateHKLToQxQy(*hkl)
    print(QxQyFromSample)
    assert(np.all(np.isclose(QxQyFromSample,[Qx,Qy],atol=1e-3)))




def test_equality():
    s1 = Sample(1,2,3,90,90,120)
    s2 = Sample(1,2,3,90,90,120)
    assert(s1==s2)

def test_calculateProjections(): # TODO: Update test

    s1 = Sample(np.pi*2,np.pi*2,np.pi*2,90,90,60)

    s1.calculateProjections()

    theta = 2*np.pi/3 # 120 degrees between projection vectors
    AStar = np.linalg.norm(s1.reciprocalVectorA)
    BStar = np.linalg.norm(s1.reciprocalVectorB)
    CStar = np.linalg.norm(s1.reciprocalVectorC)
    ## Projection is given by [[a*,cos(theta)b*],[0,sin(tetha)b*]]
    #assert(np.all(np.isclose(s1.projectionMatrix,np.array([[AStar*1,BStar*np.cos(theta)],[0,BStar*np.sin(theta)]]))))
    #assert(np.all(np.isclose(s1.tr(1,1),np.array([AStar+np.cos(theta)*BStar,np.sin(theta)*BStar]))))

    #point = s1.tr(3.5,7.2)
    #reverse = s1.inv_tr(point[0],point[1])
    #assert(np.all(np.isclose(reverse,np.array([3.5,7.2]))))

    #string = s1.format_coord(point[0],point[1])
    #assert(string=='h = 3.500, k = 7.200, l = 0.000')

def test_DataFile_Sample_UB():
    df = MJOLNIR.Data.DataFile.DataFile(os.path.join(dataPath,'camea2018n000136.hdf'))
    s = df.sample
    b1,b2,b3 = [np.linalg.norm(x) for x in [s.reciprocalVectorA,s.reciprocalVectorB,s.reciprocalVectorC]]
    a1,a2,a3 = [np.linalg.norm(x) for x in [s.realVectorA,s.realVectorB,s.realVectorC]]
    beta1,beta2,beta3 = _tools.vectorAngle(s.reciprocalVectorB,s.reciprocalVectorC),_tools.vectorAngle(s.reciprocalVectorC,s.reciprocalVectorA),_tools.vectorAngle(s.reciprocalVectorA,s.reciprocalVectorB)
    alpha1,alpha2,alpha3 = _tools.vectorAngle(s.realVectorB,s.realVectorC),_tools.vectorAngle(s.realVectorC,s.realVectorA),_tools.vectorAngle(s.realVectorA,s.realVectorB)

    cell = s.cell#[a1,a2,a3,b1,b2,b3,alpha1,alpha2,alpha3,beta1,beta2,beta3]
    r1 = s.plane_vector1
    r2 = s.plane_vector2

    UBCalc = TasUBlib.calcTasUBFromTwoReflections(cell,r1,r2)
    comparison = np.isclose(UBCalc,s.orientationMatrix,atol=1e-5) # Assert that they are equal to 1e-5 (Numerical error)
    print(np.round(UBCalc,5))
    print(np.round(s.orientationMatrix,5))
    assert(np.all(comparison))


def test_DataFile_Sample_Projection():
    df = MJOLNIR.Data.DataFile.DataFile(os.path.join(dataPath,'camea2018n000136.hdf')) # In A-B plane
    print(df.sample.projectionVector1,df.sample.projectionVector2)
    assert(np.all(np.isclose(df.sample.projectionVector1,np.array([1.0,0.0,0.0]))))
    assert(np.all(np.isclose(df.sample.projectionVector2,np.array([0.0,1.0,0.0]))))

    df2 = MJOLNIR.Data.DataFile.DataFile(os.path.join(dataPath,'camea2018n000178.hdf')) # in 1 1 0 and 0 0 1 plane
    assert(np.all(np.isclose(df2.sample.projectionVector1,np.array([0.0,0.0,1.0]))))
    assert(np.all(np.isclose(df2.sample.projectionVector2,np.array([1.0,1.0,0.0]))))



def test_Sample_CurratAxe():
    df = MJOLNIR.Data.DataFile.DataFile(os.path.join(dataPath,'camea2018n000178.hdf'))
    sample = df.sample

    Ei = [5,7,9]
    Ef = np.array([[4,4,4],[5,5,5]])
    Bragg = [[1,0,0],[0,1,1]]


    pos = sample.CurratAxe(Ei,Ef,Bragg)

    assert(pos.shape == (len(Bragg),np.array(Ei).size,Ef.size,3))

    assert(np.all(np.isclose(pos[:,:,:,2],0.0))) # All Qz are to be zero


    POS = np.array([[[[ 1.59337451e-01,  6.96316473e-01,  0.00000000e+00],
            [ 1.59337451e-01,  6.96316473e-01,  0.00000000e+00],
            [ 1.59337451e-01,  6.96316473e-01,  0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00]],

            [[ 4.35859084e-01,  7.63659414e-01,  0.00000000e+00],
            [ 4.35859084e-01,  7.63659414e-01,  0.00000000e+00],
            [ 4.35859084e-01,  7.63659414e-01,  0.00000000e+00],
            [ 2.78156829e-01,  7.17745441e-01,  0.00000000e+00],
            [ 2.78156829e-01,  7.17745441e-01,  0.00000000e+00],
            [ 2.78156829e-01,  7.17745441e-01,  0.00000000e+00]],

            [[ 6.74964311e-01,  8.21890116e-01,  0.00000000e+00],
            [ 6.74964311e-01,  8.21890116e-01,  0.00000000e+00],
            [ 6.74964311e-01,  8.21890116e-01,  0.00000000e+00],
            [ 5.18676001e-01,  7.69828565e-01,  0.00000000e+00],
            [ 5.18676001e-01,  7.69828565e-01,  0.00000000e+00],
            [ 5.18676001e-01,  7.69828565e-01,  0.00000000e+00]]],


        [[[ 1.11373538e+00,  4.08688259e-01,  0.00000000e+00],
            [ 1.11373538e+00,  4.08688259e-01,  0.00000000e+00],
            [ 1.11373538e+00,  4.08688259e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00]],

            [[ 1.33491786e+00,  2.29586062e-01,  0.00000000e+00],
            [ 1.33491786e+00,  2.29586062e-01,  0.00000000e+00],
            [ 1.33491786e+00,  2.29586062e-01,  0.00000000e+00],
            [ 1.19906946e+00,  3.22887295e-01,  0.00000000e+00],
            [ 1.19906946e+00,  3.22887295e-01,  0.00000000e+00],
            [ 1.19906946e+00,  3.22887295e-01,  0.00000000e+00]],

            [[ 1.52617192e+00,  7.47183562e-02,  0.00000000e+00],
            [ 1.52617192e+00,  7.47183562e-02,  0.00000000e+00],
            [ 1.52617192e+00,  7.47183562e-02,  0.00000000e+00],
            [ 1.38306144e+00,  1.59458179e-01,  0.00000000e+00],
            [ 1.38306144e+00,  1.59458179e-01,  0.00000000e+00],
            [ 1.38306144e+00,  1.59458179e-01,  0.00000000e+00]]]])

    assert(np.all(np.isclose(pos,POS)))

    POSAnalyser = np.array([[[[ 1.60279692e-01,  6.22804387e-01,  0.00000000e+00],
            [ 1.60279692e-01,  6.22804387e-01,  0.00000000e+00],
            [ 1.60279692e-01,  6.22804387e-01,  0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00],
            [-1.20179144e-08,  6.57512083e-01, -0.00000000e+00]],

            [[ 4.41363759e-01,  5.77272258e-01,  0.00000000e+00],
            [ 4.41363759e-01,  5.77272258e-01,  0.00000000e+00],
            [ 4.41363759e-01,  5.77272258e-01,  0.00000000e+00],
            [ 2.80013947e-01,  6.06605614e-01,  0.00000000e+00],
            [ 2.80013947e-01,  6.06605614e-01,  0.00000000e+00],
            [ 2.80013947e-01,  6.06605614e-01,  0.00000000e+00]],

            [[ 6.85994178e-01,  5.47926749e-01,  0.00000000e+00],
            [ 6.85994178e-01,  5.47926749e-01,  0.00000000e+00],
            [ 6.85994178e-01,  5.47926749e-01,  0.00000000e+00],
            [ 5.24052917e-01,  5.73796337e-01,  0.00000000e+00],
            [ 5.24052917e-01,  5.73796337e-01,  0.00000000e+00],
            [ 5.24052917e-01,  5.73796337e-01,  0.00000000e+00]]],


        [[[ 1.00477106e+00,  3.48941284e-01,  0.00000000e+00],
            [ 1.00477106e+00,  3.48941284e-01,  0.00000000e+00],
            [ 1.00477106e+00,  3.48941284e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00],
            [ 9.86285477e-01,  5.11890634e-01,  0.00000000e+00]],

            [[ 1.06290681e+00,  6.98843240e-02,  0.00000000e+00],
            [ 1.06290681e+00,  6.98843240e-02,  0.00000000e+00],
            [ 1.06290681e+00,  6.98843240e-02,  0.00000000e+00],
            [ 1.03489627e+00,  2.31469030e-01,  0.00000000e+00],
            [ 1.03489627e+00,  2.31469030e-01,  0.00000000e+00],
            [ 1.03489627e+00,  2.31469030e-01,  0.00000000e+00]],

            [[ 1.13033943e+00, -1.67701473e-01,  0.00000000e+00],
            [ 1.13033943e+00, -1.67701473e-01,  0.00000000e+00],
            [ 1.13033943e+00, -1.67701473e-01,  0.00000000e+00],
            [ 1.09633291e+00, -7.27153902e-03,  0.00000000e+00],
            [ 1.09633291e+00, -7.27153902e-03,  0.00000000e+00],
            [ 1.09633291e+00, -7.27153902e-03,  0.00000000e+00]]]])


    posAna = sample.CurratAxe(Ei,Ef,Bragg,spurionType='ANAlYser')

    assert(np.all(np.isclose(posAna[:,:,:,2],0.0))) # All Qz are to be zero
    assert(np.all(np.isclose(posAna,POSAnalyser)))


    posMonoProjection = sample.CurratAxe(Ei,Ef,Bragg,spurionType='monoChrOmaTOR',Projection=True)
    posMonoHKL = sample.CurratAxe(Ei,Ef,Bragg,spurionType='monoChrOmaTOR',HKL=True)
    assert(posMonoProjection.shape == (len(Bragg),np.array(Ei).size,Ef.size,2))
    assert(posMonoHKL.shape == (len(Bragg),np.array(Ei).size,Ef.size,3))

    try:
        sample.CurratAxe(Ei,Ef,Bragg,spurionType='WRONG!!')
        assert False
    except AttributeError:
        assert True
