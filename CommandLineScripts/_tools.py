import os
import sys
settingsFile = '.settings'



def loadSetting(string):
    returnVal=None
    if Exists():
        with open(os.path.realpath(settingsFile),'r') as f:
            lines = f.readlines()
            returnVal = None
            for l in lines:
                if string in l.split(':')[0]:
                    returnVal = l.split(':')[-1].strip()
    return returnVal

def updateSetting(name,value):
    if Exists():
        os.rename(settingsFile,settingsFile+'Old')
        lineWritten=False
        with open(settingsFile,'w') as newF:
            with open(settingsFile+'Old','r') as oldF:
                for line in oldF:
                    if line.split(':')[0]==name:
                        newF.write(name+':'+value+'\n')
                        lineWritten = True
                    else:
                        newF.write(line)
                if not lineWritten:
                    newF.write(name+':'+value+'\n')
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