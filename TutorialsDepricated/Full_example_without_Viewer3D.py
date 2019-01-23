import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')
from MJOLNIR.Data import DataSet
def test_Full_example_without_Viewer3D(show=False):
    import numpy as np
    import matplotlib.pyplot as plt

    # Convert raw data to NXSqom

    DataFile=['Data/camea2018n000137.hdf']
    DS = DataSet.DataSet(dataFiles=DataFile)
    DS.convertDataFile(saveFile=False)

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

    VMIN=0
    VMAX=20

    fig=plt.figure(figsize=(8,8))
    pc = plt.pcolormesh(Qx[:,:,Eslice].T,Qy[:,:,Eslice].T,Int[:,:,Eslice].T,vmin=VMIN,vmax=VMAX,zorder=10)
    ax = fig.add_subplot(111)


    plt.ylabel('$Q_y$ [1/AA]')
    plt.xlabel('$Q_x$ [1/AA]')
    plt.title('$\hbar \omega =$ {:.02f}'.format(np.mean(energy_bin[:,:,Eslice])) + ' meV')
    plt.axis([-1.8, 1.8, -1.8, 1.8])
    ax.set_aspect('equal', 'datalim')
    plt.grid(True)
    
    plt.colorbar(pc)
    if show:
        plt.show()

if __name__=='__main__':
    test_Full_example_without_Viewer3D(True)