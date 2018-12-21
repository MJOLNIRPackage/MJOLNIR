if False:
    import sys
    sys.path.append('..')

    from MJOLNIR.Data import Viewer3D,DataSet
    import numpy as np
    import matplotlib.pyplot as plt
    import h5py as hdf

    ConvertedDataFile=['/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000001.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000004.nxs',
                    '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000005.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000006.nxs',
                    '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000007.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000008.nxs',
                    '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000009.nxs','/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000010.nxs',
                    '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000011.nxs']
    DS = DataSet.DataSet(convertedFiles=ConvertedDataFile)

    # Extract all the data
    I,qx,qy,energy,Norm,Monitor = DS.I,DS.qx,DS.qy,DS.energy,DS.Norm,DS.Monitor

    pos=[qx,qy,energy]

    Data,bins = DataSet.binData3D(0.05,0.05,0.2,pos,I,norm=Norm,mon=Monitor)

    Intensity = np.divide(Data[0]*Data[3],Data[1]*Data[2])

    a=4.95
    astar = 2*np.pi/a

    viewer = Viewer3D.Viewer3D(Intensity,bins)
    viewer.caxis=(0.5e-6,5e-5)
    plt.show()
