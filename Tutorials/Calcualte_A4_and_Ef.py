import sys
sys.path.append('..')

from MJOLNIR.Geometry import Instrument,Detector,Analyser,Wedge

import numpy as np

Instr = Instrument.Instrument()

Det = Detector.TubeDetector1D(position=(1,0,1),direction=(1,0,0),pixels=10,split=[5],length=1.0)
Det2 = Detector.TubeDetector1D(position=(1,0.1,1),direction=(1,0,0),pixels=10,split=[5],length=1.0)
Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
Ana2 = Analyser.FlatAnalyser(position=(1,0,0),direction=(1,0,1))


wedge = Wedge.Wedge(position=(0,0,0),detectors=[Det,Det2],analysers=[Ana,Ana2])

Instr.append(wedge)


Instr.initialize()
print(Instr.A4)
print(Instr.Ef)