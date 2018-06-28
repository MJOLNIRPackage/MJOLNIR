from MJOLNIR.Data import DataSet
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mayavi import mlab
from matplotlib.colors import ListedColormap
cmap = plt.cm.jet
#my_cmap = cmap(np.arange(cmap.N))
#my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
#my_cmap = ListedColormap(my_cmap)

my_cmap = cmap

file = '/home/lass/Dropbox/PhD/Software/MJOLNIR/TestData/cameasim2018n000011.nxs'
DataObj = DataSet.DataSet(convertedFiles=file)

I = DataObj.I
qx = DataObj.qx
qy = DataObj.qy
energy = DataObj.energy
Norm = DataObj.Norm
Monitor = DataObj.Monitor


r = np.linalg.norm([qx,qy],axis=0)
theta = np.arctan2(qy,qx)

[I_bin,Monitor_bin,Normalization_bin,NormCount_bin],[r_bin,theta_bin,energy_bin] = \
DataSet.binData3D(0.01,np.deg2rad(0.5),0.5,[r.flatten(),theta.flatten(),energy.flatten()],data=I,norm=Norm,mon=Monitor)
Qx = np.cos(theta_bin)*r_bin
Qy = np.sin(theta_bin)*r_bin


Int = np.divide(I_bin*NormCount_bin,Monitor_bin*Normalization_bin)

# Plot energy slice of data
Eslice=1

if False:
    fig=plt.figure(figsize=(8,8))
    ax = plt.gca(projection='3d')
    for Eslice in range(len(energy_bin[0,0])-1):
    
        VMIN=1e-10
        VMAX=1e-6
        
        pc = ax.plot_surface(Qx[:,:,Eslice].T,Qy[:,:,Eslice].T,np.ones_like(Qy[:,:,Eslice]).T*energy_bin[0,0,Eslice],rstride=1,cstride=1,facecolors=my_cmap((Int[:,:,Eslice]/np.nanmax(Int[:,:,:])).T))#,vmin=VMIN,vmax=VMAX,zorder=10,cmap=my_cmap)
    #pc = ax.pcolormesh(Qx[:,:,Eslice].T,Qy[:,:,Eslice].T,Int[:,:,Eslice].T,vmin=VMIN,vmax=VMAX,zorder=10,cmap=my_cmap)
    ax.set_ylabel('$Q_y$ [1/AA]')
    ax.set_xlabel('$Q_x$ [1/AA]')
    plt.title('$\hbar \omega =$ {:.02f}'.format(np.mean(energy_bin[:,:,Eslice])) + ' meV')

else:
    for Eslice in range(len(energy_bin[0,0])-1):
        surf = mlab.mesh(Qx[:-1,:-1,Eslice].T,Qy[:-1,:-1,Eslice].T,energy_bin[0,0,Eslice]+Int[:,:,Eslice].T)#,warp_scale="auto")#np.ones_like(Qy[:,:,Eslice]).T*energy_bin[0,0,Eslice],color=(my_cmap((Int[:,:,Eslice].flatten()/np.nanmax(Int[:,:,:])).T)[:,:3]))
    
    xx = [np.min(Qx),np.max(Qx)]
    x0 = [xx[0],xx[0]]
    yy = [np.min(Qy),np.max(Qy)]
    y0 = [yy[0],yy[0]]
    zz = [np.min(energy_bin),np.max(energy_bin)]
    z0 = [zz[0],zz[0]]
    
    lensoffset=2.0
    mlab.plot3d(xx,y0,z0,line_width=0.01,tube_radius=0.01)# xaxis
    mlab.plot3d(x0,yy,z0,line_width=0.01,tube_radius=0.01)# yaxis
    mlab.plot3d(x0,y0,zz,line_width=0.01,tube_radius=0.01)# zaxis
    
    mlab.text3d(0.5*(np.sum(xx)), yy[0], zz[0], 'Qx [1/A]',scale=0.1)
    mlab.text3d(xx[0],0.5*(np.sum(yy)), zz[0], 'Qy [1/A]',scale=0.1)
    mlab.text3d(xx[0],yy[0],0.5*(np.sum(zz)), 'E [meV]',scale=0.1)
    
    #    mlab.plot3d(zx,zy+lensoffset,zz,line_width=0.01,tube_radius=0.01)
    #    mlab.plot3d(xx,xy+lensoffset,xz,line_width=0.01,tube_radius=0.01)
    
    #mlab.xlabel('Qx [1/A]')
    #mlab.ylabel('Qy [1/A]')
    #mlab.zlabel('E [meV]')
    mlab.show()


#plt.axis([-2.7, 2.7, -2.7, 2.7])
#ax.set_aspect('equal', 'datalim')
#plt.grid('on')
#plt.colorbar(pc)
