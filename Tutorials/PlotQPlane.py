import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')


from MJOLNIR.Data import DataSet
def test_Plot_Q_Plane(save=False):
    import numpy as np
    import matplotlib.pyplot as plt
    file = 'TestData/1024/Magnon_ComponentA3Scan.h5'
    Data = DataSet.DataSet(dataFiles=file)
    Data.convertDataFile()
    EMin = np.min(Data.energy)+1.5
    EMax = EMin+0.05

    ax = Data.plotQPlane(EMin,EMax,binning='polar',xBinTolerance=0.025,yBinTolerance=0.025,
                        enlargen=False,log=False,ax=None,RLUPlot=True,vmin=0,vmax=10)
    plt.colorbar(ax.pmeshs[0])
    ax.set_clim(0,10)

    ax2 = Data.plotQPlane(EMin,EMax,binning='xy',xBinTolerance=0.025,yBinTolerance=0.025,
                        enlargen=False,log=False,ax=None,RLUPlot=True,vmin=0,vmax=10)

    plt.colorbar(ax2.pmeshs[0])

    if save:
        fig1 = plt.figure(1)
        fig2 = plt.figure(2)
        fig1.savefig('Tutorials/PlotQPlanePolar.png')
        fig2.savefig('Tutorials/PlotQPlaneXY.png')
        plt.show()
    else:
        if os.path.exists('TestData/1024/Magnon_ComponentA3Scan.nxs'):
            os.remove('TestData/1024/Magnon_ComponentA3Scan.nxs')


if __name__ == '__main__':
    test_Plot_Q_Plane(True)