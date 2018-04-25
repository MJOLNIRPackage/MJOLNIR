import sys
sys.path.append('..')

from MJOLNIR.Data import DataSet,Viewer3D
import warnings
import matplotlib.pyplot as plt
import numpy as np
DataFile = '/home/lass/Dropbox/PhD/Software/DataSimulation/BeFilterTestOut10.h5'#['../TestData/cameasim2018n000011.h5']

dataset = DataSet.DataSet(dataFiles=DataFile)
dataset.convertDataFile(saveLocation='/home/lass/Dropbox/PhD/Software/DataSimulation/')

Data,bins = dataset.binData3D(0.08,0.08,0.25)

warnings.simplefilter('ignore')
Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
warnings.simplefilter('once')
viewer = Viewer3D.Viewer3D(Intensity,bins)
plt.plot()

