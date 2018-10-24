import sys,os
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Geometry import Instrument
def test_Generate_normalization(plot=False):
    Instr = Instrument.Instrument(fileName='TestData/1024/CAMEA_Full.xml')
    Instr.initialize()

    VanNormFile = 'TestData/1024/EScanRunDoubleFocusHS.h5'
    Instr.generateCalibration(Vanadiumdatafile=VanNormFile,savelocation='TestData/test_',plot=plot,tables=[8]) 

if __name__ == '__main__':
    test_Generate_normalization(True)