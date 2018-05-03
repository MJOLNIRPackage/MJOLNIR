import sys
sys.path.append('../')

from MJOLNIR.Data import DataSet
import numpy as np
#file = '../TestData/cameasim2018n000011.nxs'
file = ['/home/lass/Dropbox/PhD/Software/DataSimulation/T0Phonon10meV.nxs','/home/lass/Dropbox/PhD/Software/DataSimulation/T0Phonon8meV.nxs']
DataObj = DataSet.DataSet(convertedFiles=file)

energy = DataObj.energy

EnergyBins = DataSet.binEdges(energy,tolerance=0.125)
q1 = np.array([0.0,0])
q2 = np.array([0,2.0])
width = 0.1 # 1/A
minPixel = 0.01

ax,DataList,qBnLit,centerPos,binDIstance    = DataObj.plotCutQE(q1,q2,width,minPixel,EnergyBins)

ax.set_clim(0,1e-6)
