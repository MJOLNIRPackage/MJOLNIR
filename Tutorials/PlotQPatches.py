import sys
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

import matplotlib.pyplot as plt
import numpy as np
from MJOLNIR.Data import DataSet
<<<<<<< HEAD
File1 = '../TestData/NewNormalization/T0Phonon10meV.nxs'
File2 = '../TestData/NewNormalization/T0Phonon10meV93_5A4.nxs'
=======
File1 = '../TestData/T0Phonon10meV.nxs'
File2 = '../TestData/T0Phonon10meV93_5A4.nxs'
>>>>>>> 7b6969e06dab2cbbb038ab41ef7d6923c34d3ee2

DS = DataSet.DataSet(convertedFiles=[File1,File2])

files = DS.convertedFiles


<<<<<<< HEAD
planes2 = list(np.arange(64).reshape(8,8)) # Plot all planes binned with 8 pixels together
ax = [DS.createRLUAxes() for _ in range(len(planes2))] # Create custom axes for plotting

ax2 = DS.plotQPatches([files[0]],planes=planes2,ax=ax,binningDecimals=2,A4Extend=3,A3Extend=5)
=======
planes2 = list(np.arange(64).reshape(8,8))[1:] # Plot all planes binned with 8 pixels together
ax = [DS.createRLUAxes() for _ in range(len(planes2))] # Create custom axes for plotting

ax2 = DS.plotQPatches([files[0]],planes=planes2,ax=ax,binningDecimals=2,A4Extend=2,A3Extend=3)
>>>>>>> 7b6969e06dab2cbbb038ab41ef7d6923c34d3ee2

counter = 0
for ax in ax2: # loop through axes to increase size and save
    fig = ax.get_figure()
    fig.set_size_inches(10.5, 10.5, forward=True)
    fig.tight_layout()
    fig.savefig('QPatches/{:03d}.png'.format(counter),format='png')
    counter+=1