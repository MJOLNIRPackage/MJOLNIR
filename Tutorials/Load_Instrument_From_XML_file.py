import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument,Detector,Analyser,Wedge
import matplotlib.pyplot as plt
import numpy as np
Instr = Instrument.Instrument(filename='SimpleInstrument.xml') # Load XML file


fig = plt.figure() # Create 3D figure
ax = fig.gca(projection='3d')

Instr.plot(ax) # Plot instrument

ax.set_xlim(0.0,1.5)
ax.set_ylim(-0.2,0.2)
ax.set_zlim(0.0,1.1)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
plt.tight_layout()
plt.savefig('SimpleInstrument.png',format='png',dpi=300)

Instr.initialize() # Initialize instrument

plt.figure()
for det in range(len(Instr.wedges[0].detectors)):
    plt.scatter(range(Instr.wedges[0].detectors[det].pixels),
				Instr.A4[0][det]*180.0/np.pi,zorder=10,s=3)

plt.grid('on')
plt.xlabel('Pixel')
plt.ylabel('A4 [deg]')
plt.savefig('SimpleInstrument_A4.png',format='png',dpi=300)

plt.figure()
for det in range(len(Instr.wedges[0].detectors)):
    plt.scatter(range(Instr.wedges[0].detectors[det].pixels),
				Instr.Ef[0][det],zorder=10,s=3)

plt.grid('on')
plt.xlabel('Pixel')
plt.ylabel('Ef [meV]')
plt.tight_layout()
plt.savefig('SimpleInstrument_Ef.png',format='png',dpi=300)
plt.show()
