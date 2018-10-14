import sys
sys.path.append('..')

from MJOLNIR.Data import DataSet

DataFile = ['../TestData/1024/Magnon_ComponentA3Scan.h5.h5']

dataset = DataSet.DataSet(dataFiles=DataFile)
dataset.convertDataFile(saveLocation='../TestData/1024/')
