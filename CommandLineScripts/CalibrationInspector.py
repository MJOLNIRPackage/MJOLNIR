#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 14:46:53 2018

@author: Jakob Lass

Inspection tool to open the normalization tables saved in a given data file; both hd5 and nxs files.
"""


import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._pylab_helpers
import os
import sys
settingsFile = '.settings'

sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from MJOLNIR.Data import DataFile

import warnings

PlotType = {
        0: "A4",
        1: 'Normalization',
        2: 'Ef',
        3: "EfOverview",
    }


def switch(argument):
    for idx, plotType in PlotType.items():
        if plotType.upper() == argument.upper():
            return idx
    warnings.warn('Plot type "{}" not understood and will be skipped.'.format(argument))


parser = argparse.ArgumentParser(description="Inspection tool to visialize calibration tables in a data file.")
parser.add_argument("DataFile", nargs ='?', default=argparse.SUPPRESS, type=str,help="Data file from which calibration table is to be plotted. If none provided file dialog will appear.")
parser.add_argument("-s", "--save", type=str, default= '',help="Location to which the generated file will be saved.")
parser.add_argument('-p', "--plot", nargs = '*', dest = 'plotList', help = 'List of wanted plots to be generated. Should be "{}". Default all of them.'.format('","'.join([str(x) for x in list(PlotType.values())])), default = PlotType.values())
parser.add_argument("-b", "--binning", type=int, default= '8',help="Binning to be inspected. Default '8'")

args = parser.parse_args()


if not 'DataFile' in args:
    startingPath = None
    try:
        with open(os.path.realpath(settingsFile),'r') as f:
            lines = f.readlines()
            for l in lines:
                if 'CalibrationInspectorDir:' in l:
                    startingPath = l.split(':')[-1].strip()
    except:
        pass

    try:
        import tkinter as tk
        from tkinter import filedialog
    except:
        import Tkinter as tk
        import tkFileDialog as filedialog
    
    root = tk.Tk()
    root.withdraw()
    
    
    file = filedialog.askopenfilename(initialdir=startingPath, title = 'Select file for calibration plotting',filetypes=(('CAMEA Data files',('*.h5','*.nxs')),('Calibration Files','*.calib'),('All Files','*')))
    
    if len(file)==0: # No file chosen
        sys.exit()
    directory = os.path.split(file)[0]
    
    if os.path.isfile(settingsFile):
        if not startingPath == directory: # If directory of new file is different from olde one, save the new to settings
            os.rename(settingsFile,settingsFile+'Old')
            with open(settingsFile,'w') as newF:
                with open(settingsFile+'Old','r') as oldF:
                    for line in oldF:
                        if 'CalibrationInspectorDir:' in line:
                            newF.write('CalibrationInspectorDir:'+directory+'\n')
                        else:
                            newF.write(line)
            os.remove(settingsFile+'Old')
    else:
        with open(settingsFile,'w') as f:
            f.write('CalibrationInspectorDir:'+directory+'\n')


else:
    file = args.DataFile
plot = args.plotList
binning = args.binning



argsIdx = []

booleanList = np.zeros((len(PlotType)),dtype=bool)
for arg in plot:
    argsIdx.append(switch(arg))

argsIdx = np.unique([x for x in argsIdx if x is not None])

booleanList[np.array(argsIdx)]=True



if not args.save is '':
    saveFile = True
    saveLocation = args.save
else:
    saveFile = False

File = DataFile.DataFile(file)


for id in argsIdx:
    if PlotType[id] == 'A4':
        File.plotA4(binning = binning)
    
    
    if PlotType[id] == 'Normalization':
        File.plotNormalization(binning = binning)

    
    if PlotType[id] == 'Ef':
        File.plotEf(binning = binning)
    
    if PlotType[id] == 'EfOverview':
        File.plotEfOverview(binning = binning)


if saveFile==True:
    figures=[manager.canvas.figure
         for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    if len(argsIdx)>1:
        #print('__________________')
        #print(saveLocation)
        #print('__________________')
        if not os.path.isdir(saveLocation):
            raise AttributeError('Provided save location is not a directory. This is needed when multiple figures are created.')
        else:
            for f in figures:
                title = f.get_axes()[0].get_title().replace(' ','_')
                print('Save figure as '+saveLocation+title+'.png')
                f.savefig(saveLocation+title+'.png')
    else:
        title = figures[0].get_axes()[0].get_title().replace(' ','_')
        if '.' in saveLocation:
            if saveLocation.split('.')[-1]=='png':
                print('Save figure as '+saveLocation)
                f.savefig(saveLocation)
        else:
            print('Save figure as '+saveLocation+title+'.png')
            f.savefig(saveLocation+title+'.png')
            


plt.show()
    
 
