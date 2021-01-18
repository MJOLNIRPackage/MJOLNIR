"""This is the data module inteded to be used as a stand-alone data reduction and visualization. 
"""
import sys
sys.path.append('.')
#import Viewer3D
#import Viewer1D
try: # Imports needed for Windows generation of installers of MJOLNIRGui
    import numpy.random.common
    import numpy.random.bounded_integers
    import numpy.random.entropy
except ImportError:
    pass