#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 14:46:53 2018

@author: Jakob Lass

History generator tool for looking through data files and create history.
"""


import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._pylab_helpers
import os
import sys
import _tools
settingsName = 'HistoryDir'

sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataFile

import warnings



parser = argparse.ArgumentParser(description="History tool for displaying files and command for all data files in folder.")
parser.add_argument("-f", "--Folder", nargs ='?', default=argparse.SUPPRESS, type=str,help="Folder for history to be created. If none provided file dialog will appear.")
parser.add_argument("-s", "--save", type=str, default= '',help="Location to which the generated history will be saved.")

args = parser.parse_args()


if not 'Folder' in args:
    startingPath = _tools.loadSetting(settingsName)

    try:
        import tkinter as tk
        from tkinter import filedialog
    except:
        import Tkinter as tk
        import tkFileDialog as filedialog
    
    root = tk.Tk()
    root.withdraw()
    
    
    directory = filedialog.askdirectory(initialdir=startingPath, title = 'Select folder for history generation')
    
    if len(directory)==0: # No file chosen
        sys.exit()
    
    

    _tools.updateSetting(settingsName,directory)


else:
    directory = args.Folder

#
#if not args.save is '':
#    saveLocation = args.save
#else:
#    saveLocation = directory
    
 

extensions = ['.hdf','.nxs']

onlyfiles = np.sort([f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and os.path.splitext(f)[-1] in extensions])

def identicalPart(fileNames):
    from difflib import SequenceMatcher
    strings = np.array([x.split('.')[0] for x in fileNames]).copy()
    
    while(not np.all([strings[0]==x for x in strings])):
        match = [SequenceMatcher(None, strings[0], x).find_longest_match(0, len(strings[0]), 0, len(x)) for x in strings]
        strings = [strings[0][x.a: x.a + x.size] for x in match][::-1]
    return strings[0],match[0]




#string,lastMatch = identicalPart(onlyfiles)
    
completePaths = [directory+'/'+x for x in onlyfiles]

returnText = ''
#returnText +='Files, replacing {} with ^:\n'.format(string)
for file in completePaths:
    try:
        f = DataFile.DataFile(file)
    except:
        returnText+=file.split('/')[-1]+' not correct format\n'
    else:
        returnText+=f.name+': '+f.scanCommand[0].decode()+'\t'+f.title[0].decode()+'\n'

saveFile = args.save

if saveFile == '':
    print(returnText)
else:
    if saveFile.split('.')[-1]!='txt':
        saveFile+='.txt'
    if os.path.isfile(saveFile):
        oldName = '_old.'.join([x for x in saveFile.split('.')])
        print('File already exists. Renaming old file to {}'.format(oldName))
        os.rename(saveFile,oldName)
    with open(saveFile,'w') as f:
        f.write(returnText)
            






