import sys
sys.path.append('..')

from MJOLNIR.Data import Viewer3D,DataSet
import numpy as np
import matplotlib.pyplot as plt
import h5py as hdf

ConvertedDataFile=['../TestData/cameasim2018n000005.nxs','../TestData/cameasim2018n000006.nxs']

I = []
qx = []
qy = []
energy = []
Norm = []
Monitor = []

for data in ConvertedDataFile:
    
    file = hdf.File(data,'r')

    I.append(np.array(file.get('entry/data/data')))
    qx.append(np.array(file.get('entry/data/qx')))
    qy.append(np.array(file.get('entry/data/qy')))
    energy.append(np.array(file.get('entry/data/en')))
    Norm.append(np.array(file.get('entry/data/normalization')))
    Monitor.append(np.array(file.get('entry/data/monitor')))
    file.close()

I = np.concatenate(I)
qx = np.concatenate(qx)
qy = np.concatenate(qy)
energy = np.concatenate(energy)
Norm = np.concatenate(Norm)
Monitor = np.concatenate(Monitor)

pos=[qx,qy,energy]

Data,bins = DataSet.binData3D(0.075,0.075,0.1,pos,I,norm=Norm,mon=Monitor)

Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])

viewer = Viewer3D.Viewer3D(Intensity,bins)
plt.show()