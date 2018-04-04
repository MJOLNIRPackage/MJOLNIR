import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument
from MJOLNIR.Data import DataSet

Instr = Instrument.Instrument(filename='../TestData/CAMEA_Full.xml')
Instr.initialize()

NF = '../TestData/VanNormalization.h5'

dataset = DataSet.DataSet(instrument=Instr,normalizationfiles=NF)
dataset.EnergyCalibration(NF,'../TestData/',plot=True,tables=['PrismaticHighDefinition',3]) 