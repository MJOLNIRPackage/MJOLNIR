#!/usr/bin/env python
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

from MJOLNIR.CommandLineScripts import _tools
import sys,os
settingsName = 'CalibrationInspectorDir'

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


parser = argparse.ArgumentParser(description="Inspection tool to visualize calibration tables in a data file.")
parser.add_argument("DataFile", nargs ='?', default=argparse.SUPPRESS, type=str,help="Data file from which calibration table is to be plotted. If none provided file dialogue will appear.")
parser.add_argument("-s", "--save", type=str, default= '',help="Location to which the generated file will be saved.")
parser.add_argument('-p', "--plot", nargs = '*', dest = 'plotList', help = 'List of wanted plots to be generated. Should be "{}". Default all of them.'.format('","'.join([str(x) for x in list(PlotType.values())])), default = PlotType.values())
parser.add_argument("-b", "--binning", type=int, default= '8',help="Binning to be inspected. Default '8'")

args = parser.parse_args()

if not args.save is '':
    saveFile = True
    saveLocation = os.path.realpath(args.save.strip())
    if not os.path.exists(saveLocation):
        os.mkdir(saveLocation)
    else:
        if not os.path.isdir(saveLocation):
            raise FileNotFoundError("[Errno 2] No such file or directory: '{}'".format(saveLocation))
else:
    saveFile = False

file = _tools.extractDataFiles(args,settingsName,oneFile=True)[0]

directory = os.path.split(file)

_tools.updateSetting(settingsName,directory)

plot = args.plotList
binning = args.binning

argsIdx = []

booleanList = np.zeros((len(PlotType)),dtype=bool)
for arg in plot:
    argsIdx.append(switch(arg))

argsIdx = np.unique([x for x in argsIdx if x is not None])

booleanList[np.array(argsIdx)]=True



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
        if not os.path.isdir(saveLocation):
            raise AttributeError('Provided save location is not a directory. Recieved "{}"'.format(saveLocation))
        else:
            for f in figures:
                title = f.get_axes()[0].get_title().replace(' ','_')
                print('Save figure as '+saveLocation+title+'_'+str(binning)+'.png')
                f.savefig(os.path.join(saveLocation,title+'_'+str(binning)+'.png'))
    else:
        title = figures[0].get_axes()[0].get_title().replace(' ','_')
        if '.' in saveLocation:
            if saveLocation.split('.')[-1]=='png':
                print('Save figure as '+saveLocation)
                f.savefig(saveLocation)
        else:
            print('Save figure as '+saveLocation+title+'_'+str(binning)+'.png')
            f.savefig(os.path.join(saveLocation,title+'.png'))
            



def main():
    plt.show()