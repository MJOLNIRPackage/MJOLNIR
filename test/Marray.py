from MJOLNIR.Marray import Marray
import numpy as np




def test_Marray_initialize():
    Empty = Marray()
    assert(len(Empty)==0)

def test_Marray_1D():
    shortList = Marray([0,1,2,3,4])
    assert(shortList.multidimensional == False)

    shortList[0] = 10
    shortList[1] = -10
    assert(np.max(shortList) == 10)
    assert(np.min(shortList) == -10)

    shortList = 10*shortList
    assert(np.max(shortList) == 100)
    assert(np.min(shortList) == -100)

    shortList = 10+shortList-5
    assert(np.max(shortList) == 105)
    assert(np.min(shortList) == -95)
    try:
        shortList.mask = [0,0,0,0,0,0,0,0,0,0,0]
        assert False # Wrong shape for list
    except AttributeError:
        assert True

    shortList = shortList+shortList

    mask = np.zeros(len(shortList),dtype=bool)
    mask[-2:] = True
    shortList.mask = mask
    assert(np.all(shortList.mask == mask))
    assert(len(shortList.extractData())==np.sum(1-mask))
    assert(np.all(shortList.extractData() == shortList.compress()))

    shortList.mask = False
    assert(len(shortList)==5)
    assert(len(shortList.flatten())==5)

    shortList.shape = (-1)
    shortList.shape = -1
    assert(shortList.shape == (len(shortList),))

def test_Marray_multiD():
    A = Marray([np.random.rand(i,2,3) for i in [1,2,3]])
    assert(len(A) == 3)
    assert(np.all([np.all(np.isclose(AData,AD)) for AData,AD in zip(A._data,A)]))
    assert(A.multidimensional == True)

    _temp = A[-1].copy()
    del A[-1]
    assert(len(A) == 2)
    A.append(_temp)
    assert(len(A) == 3)

    string = str(A)

    A[2][0,1,1] = 20
    A[2][0,0,0] = -20

    assert(np.min(A) == -20)
    assert(np.max(A) == 20)

    assert(len(np.max(A,axis=0))==3)
    assert(len(np.min(A,axis=0))==3)
    A.mask = -1

    mask = [np.zeros(d.shape,dtype=bool) for d in A]
    mask[0][0,1,:] = True
    mask[2][1,0,1] = True
    A.mask = mask
    assert(len(A.extractData())==np.sum([np.array(M).size for M in mask])-np.sum([np.sum(M) for M in mask]))
    assert(np.all(A.extractData() == A.compress()))

    A.shape = (3,-1)
    assert(np.all(A.shape == np.array([(6,), (12,), (18,)])))
    A.shape = (3,[(1, 2, 3), (2, 2, 3), (3, 2, 3)])
    assert(np.all(A.shape == np.array([(1, 2, 3), (2, 2, 3), (3, 2, 3)])))

    0-A-A+2*A-0

    B = Marray([np.random.rand(i,2,3) for i in [1,2,3]])
    B.mask = False

    A+B
    B+A