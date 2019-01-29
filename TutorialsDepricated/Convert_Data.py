import sys
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataSet
def test_Convert_Data(save=False):
    DataFile = ['Data/camea2018n000137.hdf']

    dataset = DataSet.DataSet(dataFiles=DataFile)
    dataset.convertDataFile(saveLocation='Data/',saveFile=save)


if __name__ == '__main__':
    test_Convert_Data(True)
