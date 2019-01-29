import sys 
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from Tutorial_Class import Tutorial


def Tester():
    from MJOLNIR.Data import DataSet
    # Load and convert data
    fileName = ['/home/lass/Dropbox/PhD/CAMEAData/camea2018n000136.hdf','/home/lass/Dropbox/PhD/CAMEAData/camea2018n000137.hdf']
    ds = DataSet.DataSet(dataFiles=fileName)
    ds.convertDataFile(saveFile=False)

    # Usage of more keywords 
    ## TODO
    
title = 'Advanced View3D tutorial'

introText = 'Assuming that the Quick visualization from the Quick tutorials is understood, this tutorial seeks to '\
+'introduce more advanced features for the 3D viewer object.'

outroText = ''

introText = title+'\n'+'^'*len(title)+'\n'+introText


    
Example = Tutorial('View3D',introText,outroText,Tester,fileLocation = '/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Advanced')

def test_View3D():
    Example.test()

#if __name__ == '__main__':
Example.generateTutorial()
