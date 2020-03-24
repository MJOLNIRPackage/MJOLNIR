#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 14:46:53 2018

@author: Jakob Lass

Conversion tool for converting output h5 files to nxs files.
"""


import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib._pylab_helpers
import os
import sys
from MJOLNIR.CommandLineScripts import _tools
settingsName = 'ConvertionDir'


from MJOLNIR.Data import DataSet


parser = argparse.ArgumentParser(description="Conversion tool for converting output h5 files to nxs files.")
parser.add_argument("DataFile", nargs ='*', default=argparse.SUPPRESS, type=str,help="Data file(s) to be used. If none provided file dialogue will appear. Using string format, directory and year is also possible. See documentation.")
parser.add_argument("-s", "--save", type=str, default= '',help="Location to which the generated file will be saved.")
parser.add_argument("-b", "--binning", type=int, default= '8',help="Binning performed. Default '8'")
parser.add_argument("-r", "--reuse", action='store_true', help='Set flag to reuse files from previous usage. Default false.')

args = parser.parse_args()

files = _tools.extractDataFiles(args,settingsName)
directory = os.path.split(files[0])[0]+os.path.sep

_tools.updateSetting(settingsName,files)
binning = args.binning

if not args.save is '':
    saveLocation = args.save.strip()
else:
    saveLocation = directory
    
 
for file in files:
    dataSet = DataSet.DataSet(dataFiles = file)
    dataSet.convertDataFile(binning=binning,saveLocation=saveLocation,saveFile=True)

def main():
    pass