"""
MJOLNIR Module
^^^^^^^^^^^^^^
"""
import numpy as np

def initial(a,b,c):
    return np.linspace(a,b,c)


def test_initial():
    assert( np.all(initial(0,10,11) == [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]))