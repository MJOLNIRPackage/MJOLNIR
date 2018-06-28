import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument

Instr = Instrument.Instrument(fileName='../TestData/CAMEA_Full.xml')
Instr.initialize()

VanNormFile = '../TestData/VanNormalization.h5'
A4NormFile = '../TestData/A4Normalization.h5'

Instr.generateCalibration(Vanadiumdatafile=VanNormFile ,A4datafile=A4NormFile,savelocation='../TestData/',plot=True,tables=[8]) 
