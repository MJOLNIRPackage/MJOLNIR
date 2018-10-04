import sys
sys.path.append('..')
from MJOLNIR.Data import DataSet
import numpy as np
import matplotlib.pyplot as plt
import h5py as hdf

# Convert raw data to NXSqom

ConvertedDataFile=['/home/lass/Dropbox/PhD/Software/DataSimulation/T0Phonon10meV.nxs']
DS = DataSet.DataSet(convertedFiles=ConvertedDataFile)

# Extract all the data
I,qx,qy,energy,Norm,Monitor = DS.I,DS.qx,DS.qy,DS.energy,DS.Norm,DS.Monitor

# Reshape it 
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
DataSet.binData3D(0.01,np.deg2rad(1.0),0.5,[r.flatten(),theta.flatten(),energy.flatten()],data=I,norm=Norm,mon=Monitor)
Qx = np.cos(theta_bin)*r_bin
Qy = np.sin(theta_bin)*r_bin


Int = np.divide(I_bin*NormCount_bin,Monitor_bin*Normalization_bin)

# Plot energy slice of data
Eslice=2

VMIN=1e-10
VMAX=1e-6

fig=plt.figure(figsize=(8,8))
pc = plt.pcolormesh(Qx[:,:,Eslice].T,Qy[:,:,Eslice].T,Int[:,:,Eslice].T,vmin=VMIN,vmax=VMAX,zorder=10)
ax = fig.add_subplot(111)


plt.ylabel('$Q_y$ [1/AA]')
plt.xlabel('$Q_x$ [1/AA]')
plt.title('$\hbar \omega =$ {:.02f}'.format(np.mean(energy_bin[:,:,Eslice])) + ' meV')
plt.axis([-2.7, 2.7, -2.7, 2.7])
ax.set_aspect('equal', 'datalim')
plt.grid('on')
plt.show()
plt.colorbar(pc)
plt.show()