import os
import sys

settingsFile = os.path.join(*[x for x in os.path.split(__file__)[:-1]])
settingsFile = os.path.join(settingsFile,'.settings')

rawFileFormats = ['.hdf']
convertedFileFormats = ['.nxs']
fileFormats = rawFileFormats+convertedFileFormats



def loadSetting(string):
    returnVal=None
    if Exists():
        with open(os.path.realpath(settingsFile),'r') as f:
            lines = f.readlines()
            returnVal = None
            for l in lines:
                if string in l.split(':')[0]:
                    returnVal = l.split(':')[-1].strip()
                    if returnVal[0] in ['(','[']:
                        returnVal=tuple([x[1:-1] for x in str(returnVal)[1:-1].split(', ')])
    return returnVal

def updateSetting(name,value):
    if Exists():
        os.rename(settingsFile,settingsFile+'Old')
        lineWritten=False
        with open(settingsFile,'w') as newF:
            with open(settingsFile+'Old','r') as oldF:
                for line in oldF:
                    if line.split(':')[0]==name:
                        newF.write(name+':'+str(value)+'\n')
                        lineWritten = True
                    else:
                        newF.write(line)
                if not lineWritten:
                    newF.write(name+':'+str(value)+'\n')
        os.remove(settingsFile+'Old')
    else:
        with open(os.path.realpath(settingsFile),'w') as f:
            f.write(name+':'+value+'\n')

def Exists(file = settingsFile):
    if not os.path.isfile(settingsFile):
        with open(os.path.realpath(settingsFile),'w') as f:
            f.write('')
        return False
    else:
        return True