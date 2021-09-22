"""Settings for the interactive plotting used in the MJOLNIR package. Collected here for a common reference and easier adaptability and future-proofing."""

import numpy as np

## Viewer 3D

Viewer3DSettings = {'upwards':         ['+','up'  ,'right'],
                    'upwardsScroll':   ['up'],
                    'downwards':       ['-','down','left'],
                    'downwardsScroll': ['down'],
                    'home':     ['home'],
                    'end':      ['end'],
                    'QxE':      ['0'],
                    'QyE':      ['1'],
                    'QxQy':     ['2'],
                    }

# 1D rectangular cuts

cut1DSettings = {'new':  ['n'],
                 'move': ['m'],
                 'cut':  ['c']
                }

# Concatenate all possible settings for 1D rectangular cuts
cut1DSettingsAll = list(np.concatenate(list(cut1DSettings.values())))

# kwargs for plotting of the rectangles
cut1DKkwargs = {'linewidth':1, 'edgecolor':'r', 'facecolor':'none','zorder':22}

# Global variable to give back generated 1DCuts to the user in scripting mode
global cut1DHolder 
cut1DHolder = []   

# Colors for selected and deselected cuts
selectedColor = (1.0,0.0,0.0,1)
deselectedColor = (0.0,1.0,0.0,1)

from enum import Enum
states = {}
# States are:
#    - inactive: Not possible to draw shape
#    - initial: No point has been selected
#    - direction: Start point selected, define direction
#    - width: Start and end points selected, define width 
for i,t in enumerate(['inactive','initial','direction','width']):
    states[t.upper().replace(' ','')]=i
States= Enum('States', states)