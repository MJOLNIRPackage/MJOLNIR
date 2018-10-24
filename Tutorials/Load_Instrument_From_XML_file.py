import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Geometry import Instrument,Detector,Analyser,Wedge
def test_Load_XML(save=False):
    import matplotlib.pyplot as plt
    import numpy as np
    Instr = Instrument.Instrument(fileName='Tutorials/SimpleInstrument.xml') # Load XML file


    fig = plt.figure() # Create 3D figure
    ax = fig.gca(projection='3d')

    Instr.plot(ax) # Plot instrument

    ax.set_ylim(0.0,1.5)
    ax.set_xlim(-0.2,0.2)
    ax.set_zlim(0.0,1.1)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    plt.tight_layout()
    plt.savefig('SimpleInstrument.png',format='png',dpi=300)

    Instr.initialize() # Initialize instrument

    plt.figure()
    for det in range(len(Instr.wedges[0].detectors)):
        start,stop = Instr.wedges[0].detectors[det]._split[[0,-1]]
        plt.scatter(np.arange(start,stop),#range(Instr.wedges[0].detectors[det].pixels),
                    Instr.A4[0][det]*180.0/np.pi,zorder=10,s=3)

    plt.grid(True)
    plt.xlabel('Pixel')
    plt.ylabel('A4 [deg]')
    if save:
        plt.savefig('SimpleInstrument_A4.png',format='png',dpi=300)
    plt.figure()
    for det in range(len(Instr.wedges[0].detectors)):
        start,stop = Instr.wedges[0].detectors[det]._split[[0,-1]]
        plt.scatter(np.arange(start,stop),#Instr.wedges[0].detectors[det].pixels),
                    Instr.Ef[0][det],zorder=10,s=3)

    plt.grid(True)
    plt.xlabel('Pixel')
    plt.ylabel('Ef [meV]')
    plt.tight_layout()
    if save:
        plt.savefig('SimpleInstrument_Ef.png',format='png',dpi=300)
        plt.show()
        
if __name__ == '__main__':
    test_Load_XML(True)