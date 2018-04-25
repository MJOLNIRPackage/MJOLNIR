import sys
sys.path.append('..')

from MJOLNIR.Data import DataSet

DataFile = ['../TestData/cameasim2018n000001.h5']

dataset = DataSet.DataSet(dataFiles=DataFile)
dataset.convertDataFile(saveLocation='../TestData/')
