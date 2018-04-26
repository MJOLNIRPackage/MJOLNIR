import sys
sys.path.append('../')

from MJOLNIR.Data import DataSet
import numpy as np
file = '../TestData/cameasim2018n000011.nxs'

DataObj = DataSet.DataSet(convertedFiles=file)

energy = DataObj.energy

EnergyBins = DataSet.binEdges(energy,tolerance=0.125)
q1 = np.array([1,0])
q2 = np.array([0,1])
width = 0.1 # 1/A
minPixel = 0.01

ax,DataList,qBnLit,centerPos,binDIstance    = DataObj.plotCutQE(q1,q2,width,minPixel,EnergyBins)

ax.set_clim(0,1e-6)
