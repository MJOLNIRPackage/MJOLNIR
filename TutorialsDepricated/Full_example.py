import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataSet,Viewer3D
def test_Full_example(save=False,show=False):
    import warnings
    import matplotlib.pyplot as plt
    import numpy as np
    DataFile = 'Data/camea2018n000017.hdf'

    dataset = DataSet.DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='Data/',saveFile=save)
    viewer = dataset.View3D(0.02,0.02,0.1,rlu=False)
    
    viewer.caxis=(0,40)
    if show:
        plt.show()

if __name__=='__main__':
    test_Full_example(False,True)

