import sys
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Geometry import Instrument
def Generate_normalization(plot=False):
    Instr = Instrument.Instrument(fileName='/home/lass/Dropbox/PhD/Software/MJOLNIR/Data/CAMEA_Updated.xml')#'TestData/1024/CAMEA_Full.xml')
    Instr.initialize()

    VanNormFile = '/home/lass/Dropbox/PhD/CAMEAData/camea2018n000119.hdf'#'/home/lass/Dropbox/PhD/CAMEAData/camea2018n000084.hdf'#'/TestData/1024/camea2018n000038.hdf'#'TestData/1024/EScanRunDoubleFocusHS.h5'
    Instr.generateCalibration(Vanadiumdatafile=VanNormFile,savelocation='/home/lass/Dropbox/PhD/CAMEAData/NormalizationUpdated/',plot=plot,tables=[1,3,8]) 

if __name__ == '__main__':
    Generate_normalization(False)
    
