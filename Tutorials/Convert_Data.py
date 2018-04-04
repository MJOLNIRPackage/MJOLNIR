import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument
from MJOLNIR.Data import DataSet

Instr = Instrument.Instrument(filename='../TestData/CAMEA_Full_2.xml')
Instr.initialize()

NF = '../TestData/VanNormalization.h5'
dataset = DataSet.DataSet(instrument=Instr,normalizationfiles=NF)

normalizationfile = '../TestData/EnergyNormalization_8.calib'

DataFiles = []
for i in range(5):
	DataFiles.append('../../CAMEA_Test_Files/cameasim2018n00000'+str(i+1)+'.h5')

dataset.ConvertDatafile(datafiles=DataFiles,normalizationfile=normalizationfile)
