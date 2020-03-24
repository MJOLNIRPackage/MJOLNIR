#!/usr/bin/env python
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

settingsName = 'HistoryDir'

from MJOLNIR.Data import DataFile
from MJOLNIR.CommandLineScripts import _tools

import warnings



parser = argparse.ArgumentParser(description="History tool for displaying files and command for selected data files.")
parser.add_argument("DataFile", nargs ='*', default=argparse.SUPPRESS, type=str,help="Data file(s) to be used. If none provided file dialogue will appear. Using string format, directory and year is also possible. See documentation.")
parser.add_argument("-s", "--save", type=str, default= '',help="Location to which the generated history will be saved.")
parser.add_argument("-r", "--reuse", action='store_true', help='Set flag to reuse files from previous usage. Default false.')

args = parser.parse_args()
files = _tools.extractDataFiles(args,settingsName)
_tools.updateSetting(settingsName,files)


returnText = ''

saveFile = args.save

for file in files:
    try:
        f = DataFile.DataFile(file)
    except:
        textLine = file.split('/')[-1]+' not correct format\n'
        returnText+= textLine
    else:
        textLine = f.name+': '
        try:
            textLine+=f.scanCommand.decode()
        except:
            textLine+=str(f.scanCommand)
        try:
            textLine+='\t'+f.title.decode()+'\n'
        except:
            textLine+='\t'+str(f.title)+'\n'
        
        returnText+=textLine
    print(textLine[:-1]) # Print string excluding the added new line character. 


if saveFile != '':
    if saveFile.split('.')[-1]!='txt':
        saveFile+='.txt'
    if os.path.isfile(saveFile):
        oldName = '_old.'.join([x for x in saveFile.split('.')])
        print('File already exists. Renaming old file to {}'.format(oldName))
        os.rename(saveFile,oldName)
    with open(saveFile,'w') as f:
        f.write(returnText)
            

def main():
    pass