import sys
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')
    
def Gaussian(x,A,mu,sigma,b):
    return A*np.exp(-np.power(mu-x,2.0)*0.5*np.power(sigma,-2.0))+b

from MJOLNIR.Geometry import Instrument
def Generate_normalization(plot=False):
    Instr = Instrument.Instrument(fileName='/TestData/1024/CAMEA_Full.xml')#'TestData/1024/CAMEA_Full.xml')
    Instr.initialize()

    VanNormFile = '/TestData/1024/camea2018n000038.hdf'#'TestData/1024/EScanRunDoubleFocusHS.h5'
    Instr.generateCalibration(Vanadiumdatafile=VanNormFile,savelocation='/TestData/1024/Normalization',plot=plot,tables=[1,3,8]) 

if __name__ == '__main__':
    Generate_normalization(True)
    
