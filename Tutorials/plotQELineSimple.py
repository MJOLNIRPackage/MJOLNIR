import sys 
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from Tutorial_Class import Tutorial


def Tester():
    from MJOLNIR.Data import DataSet
    from MJOLNIR import _tools # Usefull tools useful across MJOLNIR 
    import numpy as np
    import matplotlib.pyplot as plt
    
    numbers = '483-489,494-500' # String of data numbers
    fileList = _tools.fileListGenerator(numbers,'/home/lass/Dropbox/PhD/CAMEAData/',2018) # Create file list from 2018 in specified folder
    
    ds = DataSet.DataSet(fileList)
    ds.convertDataFile(saveFile=False)
    
    # Define the positions to be cut through
    Q1 = np.array([0,0,0])
    Q2 = np.array([0,0,1])
    Q3 = np.array([-1,0,1])
    Q4 = np.array([-1,0,0])
    Q5 = np.array([0,0,1])
    # Collect them into one array
    QPoints = np.array([Q1,Q2,Q3,Q4,Q5])
    
    # Define orthogonal width and minimum pixel size along Q-cut
    width = 0.05 # 1/AA
    minPixel = 0.01 # 1/AA
    
    # Define energy bins
    Energies = np.concatenate(ds.energy,axis=0)
    EnergyBins = np.linspace(np.min(Energies),np.max(Energies),31)
    
    fig = plt.figure(figsize=(14,6))
    ax = fig.gca()
    ax,DataLists,Bins,BinCenters,Offsets = \
    ds.plotCutQELine(QPoints=QPoints, width=width, minPixel=minPixel, \
                     ax=ax, EnergyBins=EnergyBins)
    
    # Change the colorbar of the plot
    ax.set_clim(0,2e-5)
    
    fig.savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Quick/plotCutQELineMnF2.png',format='png',dpi=300)
    
title = 'Cutting through data in Q and E'

introText = 'After the initial viewing of the data acquired through `Viewer3D <Viewer3D>`_ the next step is to perform cuts. '\
+'Normally, what is wanted is the intensity measured as a function of :math:`\\vec{Q}` and energy in a 2D plot along a specified '\
+'direction through space. Often these are specified using the RLU units and questions about the number neutrons hitting a '\
+'point or if one actually has covered reciprocal space in such a fasion that an article picture can be generated. All of this can '\
+'be check through the use of the plotCutQELine method as shown below. This method creates a dedicated axis element (and figure if '\
+'needed) containing the intensity as function of position in :math:`\\vec{Q}` and energy with a given width, certain height in energy '\
+'and a minimum pixel size along the :math:`\\vec{Q}` direction. As indicated by the name, one can provide a list of :math:`\\vec{Q}` points '\
+'to be visited, with the only restriction that they should contain data and that they are connected. That is, providing three positions '\
+'results in a plot with two segments connecting point 1 and 2 as well as 2 and 3.'

outroText = '\nWith the above code the 2D cut through the MnF2 data set has been done. Is this performed in an interactive '\
+'python version using the mouse allows for further study of the data. Hovering over a given pixel reveals the position as '\
+'well as the intensity. Further, clicking (with the pointer) prints the position, intensity, and specific informations about '\
+'how many neutron counts (cts), total normalization (Norm), total monitor (Mon), and pixels binned together (NormCount) in that bin. '\
+'Many further tweaks are possible through the use of kwargs as explained in the `Advanced plotCutQELine <../Advanced/plotCutQELine.html>`_ tutorial.'\
+'\n\n'\
+'.. figure:: plotCutQELineMnF2.png\n  :width: 95%\n  :align: center\n\n'\
+'\n\nIn the code above, the function fileListGenerator is further explained in `Tools Tutorials <../Scripting.html>`_.'


introText = title+'\n'+'^'*len(title)+'\n'+introText


    
Example = Tutorial('plotQELineSimple',introText,outroText,Tester,fileLocation = '/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Quick')

def test_plotQELineSimple():
    Example.test()

#if __name__ == '__main__':
Example.generateTutorial()
