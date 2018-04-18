import sys
sys.path.append('..')

from MJOLNIR.Data import DataSet

DataFile = ['../TestData/cameasim2018n000001.h5']

dataset = DataSet.DataSet(datafiles=DataFile)
dataset.ConvertDatafile(savelocation='../TestData/')
