import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataSet,Viewer3D
def test_Full_example(save=False,show=False):
    import warnings
    import matplotlib.pyplot as plt
    import numpy as np
    DataFile = 'TestData/1024/Magnon_ComponentA3Scan.h5'

    dataset = DataSet.DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='../TestData/',saveFile=save)
    viewer = dataset.View3D(0.02,0.02,0.1)
    
    viewer.caxis=(0,40)
    if show:
        plt.show()

if __name__=='__main__':
    test_Full_example(False,True)

