import os
import sys
import pickle
import re
import MJOLNIR._tools
import numpy as np

settingsFile = os.path.join(*[x for x in os.path.split(__file__)[:-1]])
settingsFile = os.path.join(settingsFile,'.settings')

rawFileFormats = ' '.join([x for x in ['.hdf']])
convertedFileFormats = ' '.join([x for x in ['*.nxs']])
fileFormats = ' '.join([x for x in [rawFileFormats,convertedFileFormats]])


fileTypes = (('Data files',fileFormats),('Raw files',rawFileFormats),('Converted files',convertedFileFormats))


def loadSetting(string):
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

def extractDataFiles(args,settingsName,oneFile = False):
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
            try:
                import tkinter as tk
                from tkinter import filedialog
            except:
                import Tkinter as tk
                import tkFileDialog as filedialog
            
            root = tk.Tk()
            root.withdraw()
            if oneFile:
                files = filedialog.askopenfilename(initialdir=startingDir, title = 'Select file',filetypes=fileTypes)
            else:
                files = filedialog.askopenfilenames(initialdir=startingDir, title = 'Select file(s)',filetypes=fileTypes)

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