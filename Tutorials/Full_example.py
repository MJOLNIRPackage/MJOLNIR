import sys
sys.path.append('..')

#from MJOLNIR.Geometry import Instrument
from MJOLNIR.Data import DataSet,Viewer3D
import numpy as np
import matplotlib.pyplot as plt

#Instr = Instrument.Instrument(filename='/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/CAMEA_Full.xml')
#Instr.initialize()

#NF = '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/VanNormalization.h5'
#DataFile='../TestData/cameasim2018n000011.h5'

DataFile = ['/home/lass/Dropbox/PhD/Software/DataSimulation/BeFilterTestIn10.h5','/home/lass/Dropbox/PhD/Software/DataSimulation/BeFilterTestOut10.h5','/home/lass/Dropbox/PhD/Software/DataSimulation/T0Phonon10meV.h5']

dataset = DataSet.DataSet(datafiles=DataFile[2])
dataset.ConvertDatafile(savelocation='/home/lass/Dropbox/PhD/Software/DataSimulation/')

Data,bins = dataset.binData3D(0.08,0.08,0.25)


BinnedData = Data

Intensity = np.divide(BinnedData[0]*BinnedData[3],BinnedData[1]*BinnedData[2])

viewer = Viewer3D.Viewer3D(Intensity,bins)
plt.show()