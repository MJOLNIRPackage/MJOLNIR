import sys
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

import matplotlib.pyplot as plt
import numpy as np
from MJOLNIR.Data import DataSet,Viewer3D
File1 = '/home/lass/Dropbox/PhD/Software/DataSimulation/NewNormalization/T0Phonon10meV.nxs'
File2 = '/home/lass/Dropbox/PhD/Software/DataSimulation/NewNormalization/T0Phonon10meV90A4InterlaceA3.nxs'
File3 = '/home/lass/Dropbox/PhD/Software/DataSimulation/NewNormalization/T0Phonon10meV93_5A4.nxs'
File4 = '/home/lass/Dropbox/PhD/Software/DataSimulation/NewNormalization/T0Phonon10meV93_5A4InterlaceA3.nxs'

DS = DataSet.DataSet(convertedFiles=[File1,File2,File3,File4])#,dataFiles =[File2,File3])

files = DS.convertedFiles

if True:
    planes2 = list(np.arange(64).reshape(8,8)) # Plot all planes binned with 8 pixels together
    
    ax = [DS.createRLUAxes() for _ in range(len(planes2))] # Create custom axes for plotting
    
    ax2 = DS.plotQPatches(files,planes=planes2,ax=ax,A4Extend=2,A3Extend=5)
    
    counter = 0
    for ax in ax2: # loop through axes to increase size and save
        fig = ax.get_figure()
        fig.set_size_inches(10.5, 10.5, forward=True)
        fig.tight_layout()
        fig.savefig('QPatches/{:03d}.png'.format(counter),format='png')
        counter+=1
    
    plt.close('all')
#file = hdf.File(fileName,'r')

if False:
    
    I = DS.I#.flatten()#np.array(file.get('entry/data/intensity'))
    posx = DS.qx#.flatten()#np.array(file.get('entry/data/qx'))
    posy = DS.qy#.flatten()#np.array(file.get('entry/data/qy'))
    energy = DS.energy#.flatten()#np.array(file.get('entry/data/en'))
    Norm = DS.Norm#.flatten()# np.array(file.get('entry/data/normalization'))
    Monitor = DS.Monitor#.flatten()# np.array(file.get('entry/data/monitor'))
    #title = np.string_(file.get('entry/title').value)
    #file.close()
    
    pos = [posx,posy,energy]
    
    Data,bins = DataSet.binData3D(0.02,0.02,0.2,pos,I,norm=Norm,mon=Monitor)
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    
    Viewer = Viewer3D.Viewer3D(Intensity,bins,axis=2)
    #Viewer.ax.set_title(str(title)[2:-1])
    #plt.show()