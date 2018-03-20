import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument,Detector,Analyser,Wedge

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

Instr = Instrument.Instrument()

Det = Detector.TubeDetector1D(position=(1,0,1),direction=(1,0,0))
Ana = Analyser.FlatAnalyser(position=(1,0,0),direction=(1,0,1))

wedge = Wedge.Wedge(position=(0,0,0),detectors=Det,analysers=Ana)

Instr.append(wedge)

fig = plt.figure()
ax = fig.gca(projection='3d')

Instr.plot(ax)

ax.set_xlim(-0.1,1.1)
ax.set_ylim(-0.1,1.1)
ax.set_zlim(-0.1,1.1)

plt.savefig('../docs/_templates/Build_a_simple_instrument.png',format='png',dpi=300)
plt.show()