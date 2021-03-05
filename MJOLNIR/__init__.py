"""
MJOLNIR Module
^^^^^^^^^^^^^^

"""
#import Geometry, Statistic, Data

import sys,os
sys.path.append('.')

__version__ = '1.1.16'
__author__ = 'Jakob Lass'

__multiFLEXXNormalization__ = os.path.join(os.path.dirname(__file__),'CalibrationMultiFLEXX.csv') 
__flatConeNormalization__ = os.path.join(os.path.dirname(__file__),'CalibrationFlatCone.csv')
