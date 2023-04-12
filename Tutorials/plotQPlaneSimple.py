import sys 
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from Tutorial_Class import Tutorial


def Tester():
    from MJOLNIR.Data import DataSet
    from MJOLNIR import _tools # Usefull tools useful across MJOLNIR 
    import numpy as np
    
    numbers = '483-489,494-500' # String of data numbers
    fileList = _tools.fileListGenerator(numbers,'/home/lass/Dropbox/PhD/CAMEAData/',2018) # Create file list from 2018 in specified folder
    
    # Create the data set
    ds = DataSet.DataSet(fileList)
    ds.convertDataFile(saveFile=False)
        
    # Choose energy limits for binning
    EMin = 3.5
    EMax = 4.0
    # Generate a figure making use of binning in polar coordinates
    Data,ax = ds.plotQPlane(EMin=EMin, EMax=EMax,xBinTolerance=0.03,yBinTolerance=0.03,binning='polar',vmin=2e-7,vmax=1e-5)
    
    fig = ax.get_figure() # Extract figure from returned axis
    fig.colorbar(ax.pmeshs[0]) # Create colorbar from plot
    ax.set_xlim(-2.6,0.68)
    ax.set_ylim(-1.76,2.58)
    fig.set_size_inches(4.3,4)
    fig.savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Quick/ConstantEnergy3DPolar.png',format='png',dpi=300)
    
    # Generate a figure making use of binning in regular orthonormal coordinates
    Data2,ax2 =  ds.plotQPlane(EMin=EMin, EMax=EMax,xBinTolerance=0.03,yBinTolerance=0.03,binning='xy',vmin=2e-7,vmax=1e-5)

    fig2 = ax2.get_figure()
    fig2.colorbar(ax2.pmeshs[0])
    fig2.set_size_inches(4.3,4)
    fig2.savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Quick/ConstantEnergy3DXY.png',format='png',dpi=300)
    
    
title = 'Plotting a Q plane for constant energy'

introText = 'After one has gotten an overview of the data measured through the Viewer3D (explained in `<QuickView3D.html>`_ '\
+'the next step is to plot only as single plane of constant energy as function of H, K, and L or :math:`Q_x` and :math:`Q_y` '\
+'depending on the boolean state of the "rlu" key word argument. Two different binning methods are currently provided: Polar '\
+'and XY. What is done is that the points measured are binned either. For explanation of units and size of bins, see below'

outroText = 'The above code creates the two figures shown below. The difference between the two is that the former performs the '\
'creation of bins in polar coordinates. This means that the provided bin tolerances are to be understood as the angular and radial bin '\
+'sizes respectively for x and y. If the rlu is set to true the units are in projections along the in plane scattering vectors. Otherwise '\
+'units are in 1/AA.'\
+'\n\n'\
+'|pic1| |pic2|\n\n'\
+'.. |pic1| image:: ConstantEnergy3DPolar.png\n  :width: 45%\n\n'\
+'.. |pic2| image:: ConstantEnergy3DXY.png\n  :width: 45%\n\n'\
+'\n\nFor further examples and the usage of 3D axis in this method, see `Plotting Q planes for constant energy <../Advanced/ConstantEnergy.html>`_.'\
+'\n\nBinnings explained\n'+18*'-'+'\n\n'\
+'The bin sizes depends on the other parameters provided to the method. The tables below seeks to show all of the possibilities:\n\n'\
+'+-------------------------+---------------+------+----------------------+---------------------------+\n'\
+'| Binning with rlu==False | Name          | Unit | Limits               | Explanation               |\n'\
+'+-------------------------+---------------+------+----------------------+---------------------------+\n'\
+'| XY                      | xBinTolerance | 1/AA | (0, :math:`\infty` ) | Binning along :math:`Q_x` |\n'\
+'+-------------------------+---------------+------+----------------------+---------------------------+\n'\
+'| XY                      | yBinTolerance | 1/AA | (0, :math:`\infty` ) | Binning along :math:`Q_y` |\n'\
+'+-------------------------+---------------+------+----------------------+---------------------------+\n'\
+'| Polar                   | xBinTolerance | rad  | (0, :math:`2\pi` ]   | Angular binning           |\n'\
+'+-------------------------+---------------+------+----------------------+---------------------------+\n'\
+'| Polar                   | yBinTolerance | 1/AA | (0, :math:`\infty` ) | Radial binning            |\n'\
+'+-------------------------+---------------+------+----------------------+---------------------------+\n'\
+'\n\nWith the rlu==True shown below, it is assumed that data is plotted in the RLU axis container.\n\n'\
+"+------------------------+---------------+------+----------------------+---------------------------------+\n"\
+"| Binning with rlu==True | Name          | Unit | Limits               | Explanation                     |\n"\
+"+------------------------+---------------+------+----------------------+---------------------------------+\n"\
+"| XY                     | xBinTolerance | N/A  | (0, :math:`\infty` ) | Binning along first projection  |\n"\
+"+------------------------+---------------+------+----------------------+---------------------------------+\n"\
+"| XY                     | yBinTolerance | N/A  | (0, :math:`\infty` ) | Binning along second projection |\n"\
+"+------------------------+---------------+------+----------------------+---------------------------------+\n"\
+"| Polar                  | xBinTolerance | rad  | (0, :math:`2\pi` ]   | Angular binning                 |\n"\
+"+------------------------+---------------+------+----------------------+---------------------------------+\n"\
+"| Polar                  | yBinTolerance | N/A  | (0, :math:`\infty` ) | Radial binning                  |\n"\
+"+------------------------+---------------+------+----------------------+---------------------------------+\n"\
+'\nFor further explanation of the RLU axis see `Reciprocal Lattice Unit Axis <../Tools/RLUAxis.html>`_.'

introText = title+'\n'+'^'*len(title)+'\n'+introText


    
Example = Tutorial('ConstantEnergy',introText,outroText,Tester,fileLocation = '/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Quick')

def test_ConstantEnergy():
    Example.test()

#if __name__ == '__main__':
Example.generateTutorial()
