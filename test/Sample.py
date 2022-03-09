from MJOLNIR import TasUBlibDEG
from MJOLNIR.Data.Sample import Sample
import MJOLNIR.Data.DataFile
from MJOLNIR import _tools
import os
import numpy as np
import MJOLNIR.TasUBlibDEG as TasUBlib

dataPath = 'samlpedata'



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


def test_parameters():
    s1 = Sample(1,2,3,90,90,120)
    pars = np.array([getattr(s1,x) for x in ['a','b','c','alpha','beta','gamma']])

    assert(np.all(np.isclose(pars,np.array([1,2,3,90,90,120]))))

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


    POS = np.array([[[[ 5.99262072e-01, -3.88754450e-01,  0.00000000e+00],
         [ 5.99262072e-01, -3.88754450e-01,  0.00000000e+00],
         [ 5.99262072e-01, -3.88754450e-01,  0.00000000e+00],
         [ 6.17531540e-01, -2.25780728e-01,  0.00000000e+00],
         [ 6.17531540e-01, -2.25780728e-01,  0.00000000e+00],
         [ 6.17531540e-01, -2.25780728e-01,  0.00000000e+00]],

        [[ 5.67556390e-01, -6.71586618e-01,  0.00000000e+00],
         [ 5.67556390e-01, -6.71586618e-01,  0.00000000e+00],
         [ 5.67556390e-01, -6.71586618e-01,  0.00000000e+00],
         [ 5.78587070e-01, -5.07707339e-01,  0.00000000e+00],
         [ 5.78587070e-01, -5.07707339e-01,  0.00000000e+00],
         [ 5.78587070e-01, -5.07707339e-01,  0.00000000e+00]],

        [[ 5.40140833e-01, -9.16148495e-01,  0.00000000e+00],
         [ 5.40140833e-01, -9.16148495e-01,  0.00000000e+00],
         [ 5.40140833e-01, -9.16148495e-01,  0.00000000e+00],
         [ 5.44912214e-01, -7.51486190e-01,  0.00000000e+00],
         [ 5.44912214e-01, -7.51486190e-01,  0.00000000e+00],
         [ 5.44912214e-01, -7.51486190e-01,  0.00000000e+00]]],


       [[[ 1.39600635e-03, -1.18635182e+00,  0.00000000e+00],
         [ 1.39600635e-03, -1.18635182e+00,  0.00000000e+00],
         [ 1.39600635e-03, -1.18635182e+00,  0.00000000e+00],
         [ 1.42087648e-01, -1.10208991e+00,  0.00000000e+00],
         [ 1.42087648e-01, -1.10208991e+00,  0.00000000e+00],
         [ 1.42087648e-01, -1.10208991e+00,  0.00000000e+00]],

        [[-2.42766804e-01, -1.33258385e+00,  0.00000000e+00],
         [-2.42766804e-01, -1.33258385e+00,  0.00000000e+00],
         [-2.42766804e-01, -1.33258385e+00,  0.00000000e+00],
         [-1.08490328e-01, -1.23703420e+00,  0.00000000e+00],
         [-1.08490328e-01, -1.23703420e+00,  0.00000000e+00],
         [-1.08490328e-01, -1.23703420e+00,  0.00000000e+00]],

        [[-4.53891710e-01, -1.45902909e+00,  0.00000000e+00],
         [-4.53891710e-01, -1.45902909e+00,  0.00000000e+00],
         [-4.53891710e-01, -1.45902909e+00,  0.00000000e+00],
         [-3.25162358e-01, -1.35371906e+00,  0.00000000e+00],
         [-3.25162358e-01, -1.35371906e+00,  0.00000000e+00],
         [-3.25162358e-01, -1.35371906e+00,  0.00000000e+00]]]])

    assert(np.all(np.isclose(pos,POS)))

    POSAnalyser = np.array([[[[ 0.52989639, -0.36439634,  0.        ],
         [ 0.52989639, -0.36439634,  0.        ],
         [ 0.52989639, -0.36439634,  0.        ],
         [ 0.61753154, -0.22578073,  0.        ],
         [ 0.61753154, -0.22578073,  0.        ],
         [ 0.61753154, -0.22578073,  0.        ]],

        [[ 0.39061242, -0.61275375,  0.        ],
         [ 0.39061242, -0.61275375,  0.        ],
         [ 0.39061242, -0.61275375,  0.        ],
         [ 0.47356748, -0.47128762,  0.        ],
         [ 0.47356748, -0.47128762,  0.        ],
         [ 0.47356748, -0.47128762,  0.        ]],

        [[ 0.27904852, -0.83243237,  0.        ],
         [ 0.27904852, -0.83243237,  0.        ],
         [ 0.27904852, -0.83243237,  0.        ],
         [ 0.35895352, -0.68922135,  0.        ],
         [ 0.35895352, -0.68922135,  0.        ],
         [ 0.35895352, -0.68922135,  0.        ]]],


       [[[-0.01730114, -1.06349686,  0.        ],
         [-0.01730114, -1.06349686,  0.        ],
         [-0.01730114, -1.06349686,  0.        ],
         [ 0.14208765, -1.10208991,  0.        ],
         [ 0.14208765, -1.10208991,  0.        ],
         [ 0.14208765, -1.10208991,  0.        ]],

        [[-0.29935284, -1.02227325,  0.        ],
         [-0.29935284, -1.02227325,  0.        ],
         [-0.29935284, -1.02227325,  0.        ],
         [-0.13797499, -1.05145191,  0.        ],
         [-0.13797499, -1.05145191,  0.        ],
         [-0.13797499, -1.05145191,  0.        ]],

        [[-0.54564749, -1.00402183,  0.        ],
         [-0.54564749, -1.00402183,  0.        ],
         [-0.54564749, -1.00402183,  0.        ],
         [-0.38329525, -1.02717256,  0.        ],
         [-0.38329525, -1.02717256,  0.        ],
         [-0.38329525, -1.02717256,  0.        ]]]])


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
