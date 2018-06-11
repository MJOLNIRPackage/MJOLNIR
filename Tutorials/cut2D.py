import sys
sys.path.append('../')

from MJOLNIR.Data import DataSet
import numpy as np
import matplotlib.pyplot as plt
file = '../TestData/cameasim2018n000011.nxs'
#file = ['/home/lass/Dropbox/PhD/Software/DataSimulation/T0Phonon10meV.nxs','/home/lass/Dropbox/PhD/Software/DataSimulation/T0Phonon8meV.nxs']
DataObj = DataSet.DataSet(convertedFiles=file)

energy = DataObj.energy

EnergyBins = DataSet.binEdges(energy,tolerance=0.125)
q1 = np.array([1.0,0])
q2 = np.array([0,1.0])
width = 0.1 # 1/A
minPixel = 0.01

ax,DataList,qBnLit,centerPos,binDIstance = DataObj.plotCutQE(q1,q2,width,minPixel,EnergyBins)
plt.colorbar(ax.pmeshs[0])

ax.set_clim(0,1e-6)
plt.tight_layout()

## Cut and plot 1D
ax,DataList,Bins,binCenter,binDistance = DataObj.plotCut1D(q1,q2,width,minPixel,Emin = 5.2, Emax = 5.7,plotCoverage=True)