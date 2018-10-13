import sys
sys.path.append('../')

from MJOLNIR.Data import DataSet
import numpy as np
import matplotlib.pyplot as plt
file = '../TestData/1024/Magnon_ComponentA3Scan.nxs'

Data = DataSet.DataSet(convertedFiles=file)
EMin = np.min(Data.energy)+1.5
EMax = EMin+0.05

ax = Data.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.025,yBinTolerance=0.025,
                     enlargen=False,log=False,ax=None,RLUPlot=True,vmin=0,vmax=10)
plt.colorbar(ax.pmeshs[0])
ax.set_clim(0,10)

ax2 = Data.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.025,yBinTolerance=0.025,
                     enlargen=False,log=False,ax=None,RLUPlot=True,vmin=0,vmax=10)

plt.colorbar(ax2.pmeshs[0])
plt.show()

fig = plt.figure(1)
fig.savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/Tutorials/PlotQPlanePolar.png')
fig = plt.figure(2)
fig.savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/Tutorials/PlotQPlaneXY.png')