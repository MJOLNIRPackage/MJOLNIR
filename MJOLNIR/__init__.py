"""
MJOLNIR Module
^^^^^^^^^^^^^^

"""
#import Geometry, Statistic, Data

import sys,os
sys.path.append('.')

__version__ = '1.3.5'
__author__ = 'Jakob Lass'

__multiFLEXXNormalization__ = os.path.join(os.path.dirname(__file__),'CalibrationMultiFLEXX.csv') 
__flatConeNormalization__ = os.path.join(os.path.dirname(__file__),'CalibrationFlatCone.csv')
__CAMEANormalizationBinning1__ = os.path.join(os.path.dirname(__file__), 'Normalization_1.calib')
__CAMEANormalizationBinning2__ = os.path.join(os.path.dirname(__file__), 'Normalization_2.calib')
__CAMEANormalizationBinning3__ = os.path.join(os.path.dirname(__file__), 'Normalization_3.calib')
__CAMEANormalizationBinning4__ = os.path.join(os.path.dirname(__file__), 'Normalization_4.calib')
__CAMEANormalizationBinning5__ = os.path.join(os.path.dirname(__file__), 'Normalization_5.calib')
__CAMEANormalizationBinning6__ = os.path.join(os.path.dirname(__file__), 'Normalization_6.calib')
__CAMEANormalizationBinning7__ = os.path.join(os.path.dirname(__file__), 'Normalization_7.calib')
__CAMEANormalizationBinning8__ = os.path.join(os.path.dirname(__file__), 'Normalization_8.calib')
__bambusNormalization__ = os.path.join(os.path.dirname(__file__),'CalibrationBambus.csv')
