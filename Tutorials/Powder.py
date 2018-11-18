import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')
from MJOLNIR.Data import DataFile
DataFile.assertFile('Data/Magnon_ComponentA3Scan.nxs')
from MJOLNIR.Data import DataSet
from MJOLNIR import _tools
def test_Powder(show=False):
    import matplotlib.pyplot as plt
    file = 'Data/Magnon_ComponentA3Scan.hdf'

    DataObj = DataSet.DataSet(dataFiles=file)
    DataObj.convertDataFile()
    I = DataObj.I
    qx = DataObj.qx
    qy = DataObj.qy
    energy = DataObj.energy
    Norm = DataObj.Norm
    Monitor = DataObj.Monitor

    EBinEdges = _tools.binEdges(energy,tolerance=0.125)

    ax,Data,qbins = DataObj.plotCutPowder(EBinEdges,qMinBin=0.05)
    plt.colorbar(ax.pmeshs[0])

    ax2,Data2,qbins2 = DataSet.plotCutPowder([qx,qy,energy],I,Norm,Monitor,EBinEdges,qMinBin=0.05)
    plt.colorbar(ax2.pmeshs[0])

    Data3,qbins3 = DataObj.cutPowder(EBinEdges)

    ax2.set_clim(0,0.01)
    if show:
        plt.show()
    else:
        if os.path.exists('TestData/1024/Magnon_ComponentA3Scan.nxs'):
            os.remove('TestData/1024/Magnon_ComponentA3Scan.nxs')

if __name__=='__main__':
    test_Powder(True)