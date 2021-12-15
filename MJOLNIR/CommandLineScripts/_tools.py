import os
import sys
import pickle
import re
import MJOLNIR._tools
from MJOLNIR.Data import DataFile
import numpy as np
from PyQt5.QtWidgets import QApplication, QWidget, QFileDialog
from PyQt5.Qt import QApplication
from PyQt5 import Qt

from os.path import expanduser
settingsFile = expanduser("~") # Use home folder for storing settings file

settingsFile = os.path.join(settingsFile,'.MJOLNIRsettings')

allowedRawFilesString = 'Raw ('+' '.join(['*.'+str(x) for x in DataFile.supportedRawFormats])+')'
allowedConvertedFilesString = 'Converted ('+' '.join(['*.'+str(x) for x in DataFile.supportedConvertedFormats])+')'
allowedAllFilesString = 'All Files (*)'

allowedString = ';;'.join([allowedRawFilesString,allowedConvertedFilesString,allowedAllFilesString])
global dataFilesLoaded
# All functions are tested through the test_coonadline.py script but are called through the terminal thus not recorded as tested.
def loadSetting(string): # pragma: no cover
    if Exists():
        with open(settingsFile,"rb") as pickle_in:
            loaded_dict = pickle.load(pickle_in)
        if string in loaded_dict:
            returnValue = loaded_dict[string]
        else:
            returnValue = None
        return returnValue
    else:
        return None

def updateSetting(name,value):# pragma: no cover
    if Exists():
        with open(settingsFile,"rb") as pickle_in:
            loaded_dict = pickle.load(pickle_in)
    else:
        loaded_dict = {}
    loaded_dict[name]=value
    with open(settingsFile,"wb") as pickle_out:
        pickle.dump(loaded_dict, pickle_out)

def Exists(file = settingsFile):# pragma: no cover
    return os.path.isfile(settingsFile)

def extractDataFiles(args,settingsName,oneFile = False):# pragma: no cover
    if not 'DataFile' in args:
        startingPath = loadSetting(settingsName)
        if not startingPath is None:
            if isinstance(startingPath,tuple):
                startingDir = '/'.join([x for x in list(startingPath)[0].split('/')[:-1]])
            else:
                startingDir = startingPath
        else:
            startingDir = None
        
        if 'reuse' in args:
            reuse = args.reuse
        else:
            reuse = False

        if reuse:
            files = startingPath
        else:
            files = loadFiles(startingDir)
        files = tuple(np.unique(files))
        if len(list(files))==0: # No file chosen
            sys.exit()

        directory = os.path.split(files[0])[0]


    else:
        pattern = re.compile(r'^(\d+[-,]?){1,}\d') # Generate pattern to recognize format: Start with digits, and then - or . or nothing repeated at least once, and ending in digit.
        res = pattern.match(args.DataFile[0])
        if not res is None:
            if res.start()==0 and res.end() == len(args.DataFile[0]): # Correct expression is given
                if len(args.DataFile) == 3: # Provided is string, directory and year
                    fileNumbers, directory, year = args.DataFile
                    year = int(year)
                    files = MJOLNIR._tools.fileListGenerator(fileNumbers,directory,year)
                elif len(args.DataFile) == 4:
                    raise NotImplementedError('Using non-CAMEA input is currently not supported... sorry :-/ ')
                else:
                    raise AttributeError('Provided data file arguments does not fit the scheme of being number list, directory, and year! (+ possibly instrument)')
        else:
            files = args.DataFile
    return files

class App(QWidget):
    def __init__(self,startingPath=None):
        super().__init__()
        self.title = 'MJOLNIR Command Line - Load Files'
        self.startingPath = startingPath
        self.initUI()
    
    def initUI(self):
        self.hide()
        self.setWindowTitle(self.title)
        self.openFileNamesDialog()
        
        
        Qt.QCoreApplication.quit()
    
    def openFileNamesDialog(self):
        #options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        #print(self,"MJOLNIR Command Line - Load Files", allowedString,self.startingPath)
        files, _ = QFileDialog.getOpenFileNames(self,"MJOLNIR Command Line - Load Files", self.startingPath,allowedString)
        global dataFilesLoaded 
        dataFilesLoaded = files
        self.close()
    

    
def loadFiles(startingPath=None):
    app = QApplication(sys.argv)
    ex = App(startingPath=startingPath)
    return dataFilesLoaded