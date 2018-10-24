import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataSet,Viewer3D
def test_Binning_data(view = False):
    import numpy as np
    import h5py as hdf
    import matplotlib.pyplot as plt
    fileName = 'TestData/1024/Magnon_ComponentA3Scan.h5'
    ds = DataSet.DataSet(dataFiles=fileName)
    ds.convertDataFile()

    I = ds.convertedFiles[0].I
    qx = ds.convertedFiles[0].qx
    qy = ds.convertedFiles[0].qy
    energy = ds.convertedFiles[0].energy
    Norm = ds.convertedFiles[0].Norm
    Monitor = ds.convertedFiles[0].Monitor
    title = 'Magnon ComponentA3Scan'

    pos = [qx,qy,energy]

    Data,bins = DataSet.binData3D(0.02,0.02,0.1,pos,I,norm=Norm,mon=Monitor)
    import warnings
    warnings.simplefilter("ignore")
    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])
    warnings.simplefilter('once')

    Viewer = Viewer3D.Viewer3D(Intensity,bins,axis=2)

    Viewer.caxis=(0,40)

    Viewer.ax.set_title(str(title)[2:-1])
    if view:
        plt.show()
    else:
        if os.path.exists('TestData/1024/Magnon_ComponentA3Scan.nxs'):
            os.remove('TestData/1024/Magnon_ComponentA3Scan.nxs')

if __name__=='__main__':
    test_Binning_data(True)

