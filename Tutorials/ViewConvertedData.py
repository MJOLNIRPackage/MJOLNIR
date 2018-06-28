import sys
sys.path.append('..')

from MJOLNIR.Data import Viewer3D,DataSet
import numpy as np
import matplotlib.pyplot as plt
import h5py as hdf

ConvertedDataFile=['/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000001.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000004.nxs',
                   '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000005.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000006.nxs',
                   '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000007.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000008.nxs',
                   '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000009.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000010.nxs',
                   '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000011.nxs']

I = []
qx = []
qy = []
energy = []
Norm = []
Monitor = []

for data in ConvertedDataFile:
    
    file = hdf.File(data,'r')

    I.append(np.array(file.get('entry/data/intensity')))
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

Data,bins = DataSet.binData3D(0.05,0.05,0.2,pos,I,norm=Norm,mon=Monitor)

Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])

a=4.95
astar = 2*np.pi/a

viewer = Viewer3D.Viewer3D(Intensity,bins)
#for i in np.arange(0,2.1,2):
#    for j in np.arange(0,2.1,2):
#        viewer.ax.scatter(i*astar,j*astar)
viewer.caxis=(0.5e-6,5e-5)
plt.show()