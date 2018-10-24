import sys
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataSet
def test_Convert_Data(save=False):
    DataFile = ['TestData/1024/Magnon_ComponentA3Scan.h5']

    dataset = DataSet.DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='TestData/1024/',saveFile=save)


if __name__ == '__main__':
    test_Convert_Data(True)
