import os
import sys
import pickle

settingsFile = os.path.join(*[x for x in os.path.split(__file__)[:-1]])
settingsFile = os.path.join(settingsFile,'.settings')

rawFileFormats = ['.hdf']
convertedFileFormats = ['.nxs']
fileFormats = rawFileFormats+convertedFileFormats



def loadSetting(string):
    if Exists():
        with open(settingsFile,"rb") as pickle_in:
            loaded_dict = pickle.load(pickle_in)
        returnValue = loaded_dict[string]
        return returnValue
    else:
        return None

def updateSetting(name,value):
    if Exists():
        with open(settingsFile,"rb") as pickle_in:
            loaded_dict = pickle.load(pickle_in)
    else:
        loaded_dict = {}
    loaded_dict[name]=value
    with open(settingsFile,"wb") as pickle_out:
        pickle.dump(loaded_dict, pickle_out)

def Exists(file = settingsFile):
    return os.path.isfile(settingsFile)