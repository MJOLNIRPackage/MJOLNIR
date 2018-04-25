import sys
sys.path.append('..')
from MJOLNIR.Data import DataSet
import numpy as np
import matplotlib.pyplot as plt
import h5py as hdf

# Convert raw data to NXSqom

Datapath='path_to_data/'
data=['066274','066275','066276','066277','066278','066279','066280']

DataFile=[]
hfend='.h5'
for i in range(0,len(data)):
    DataFile.append(Datapath+data[i]+hfend)
    
dataset = DataSet.DataSet(dataFiles=DataFile)
dataset.convertDataFile(binning=1)

ConvertedDataFile=[]
nxsend='.nxs'
for i in range(0, len(data)):
    ConvertedDataFile.append(Datapath+data[i]+nxsend)

# Extract intensities and positions from files
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

# Bin data in polar coordinates

r = np.linalg.norm([qx,qy],axis=0)
theta = np.arctan2(qy,qx)

[I_bin,Monitor_bin,Normalization_bin,NormCount_bin],[r_bin,theta_bin,energy_bin] = \
DataSet.binData3D(0.04,np.deg2rad(2.0),0.5,[r,theta,energy],data=I,norm=Norm,mon=Monitor)

Qx = np.cos(theta_bin)*r_bin
Qy = np.sin(theta_bin)*r_bin

Int = np.divide(I_bin*NormCount_bin,Monitor_bin*Normalization_bin)

# Plot energy slice of data

Eslice=1

VMIN=1e-10
VMAX=1e-7

fig=plt.figure(figsize=(8,8))
plt.pcolormesh(Qx[:,:,Eslice].T,Qy[:,:,Eslice].T,Int[:,:,Eslice].T,vmin=VMIN,vmax=VMAX)
ax = fig.add_subplot(111)


plt.ylabel('$q_y$, [1/AA]')
plt.xlabel('$q_x$, [1/AA]')
plt.title('$\hbar \omega =$' + str(Eslice) + ' meV')
plt.axis([-1, 2.5, -2.7, 2.7])
ax.set_aspect('equal', 'datalim')

plt.show()
