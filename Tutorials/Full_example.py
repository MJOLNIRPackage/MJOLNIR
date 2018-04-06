import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument
from MJOLNIR.Data import DataSet,Viewer3D
import numpy as np
import matplotlib.pyplot as plt

Instr = Instrument.Instrument(filename='../TestData/CAMEA_Full.xml')
Instr.initialize()

NF = '../TestData/VanNormalization.h5'
DataFile='../TestData/cameasim2018n000005.h5'

dataset = DataSet.DataSet(instrument=Instr,normalizationfiles=NF,datafiles=DataFile)
dataset.EnergyCalibration(tables=[8],savelocation='./')
dataset.ConvertDatafile(savelocation='./')

Data,bins = dataset.binData3D(0.05,0.05,0.2,)


BinnedData = Data

Intensity = np.divide(BinnedData[0]*BinnedData[3],BinnedData[1]*BinnedData[2])

viewer = Viewer3D.Viewer3D(Intensity,bins)
plt.show()