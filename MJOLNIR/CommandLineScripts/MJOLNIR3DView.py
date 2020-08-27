#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 17:24:53 2019

@author: Jakob Lass

Tool for quick viewing of data with viewer3D
"""


import argparse
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
from MJOLNIR.CommandLineScripts import _tools
settingsName = 'ConvertionDir'


from MJOLNIR.Data import DataSet
import MJOLNIR._tools
import re

import warnings



parser = argparse.ArgumentParser(description="Conversion tool for quick visualization using the viewer3D.")
parser.add_argument("DataFile", nargs ='*', default=argparse.SUPPRESS, type=str,help="Data file(s) to be used. If none provided file dialogue will appear. Using string format, directory and year is also possible. See documentation.")
parser.add_argument("-r", "--reuse", action='store_true', help='Set flag to reuse files from previous usage. Default false.')
parser.add_argument("-b", "--binning", type=int, default= '8',help="Binning performed. Default '8'")
parser.add_argument("-d", "--dQxdQydE", nargs=3, type=float, default=[0.03,0.03,0.08], help="Binning used to plot in 3D, Default [0.03,0.03,0.08]")
parser.add_argument("-M", "--VMax", type=float, default=None, help='Maximal value for plotting, default max of data')
parser.add_argument("-m", "--VMin", type=float, default=None, help='Minimal value for plotting, default min of data')


args = parser.parse_args()

files = _tools.extractDataFiles(args,settingsName)

_tools.updateSetting(settingsName,files)
binning = args.binning
dQxdQydE = args.dQxdQydE

files = list(files)

dataSet = DataSet.DataSet(np.array(files).copy())
dataSet.convertDataFile(binning=binning,saveFile=False)
dQx, dQy, dE = dQxdQydE
V = dataSet.View3D(dQx=dQx,dQy=dQy,dE=dE,grid=-20,rlu=True)

VMax = args.VMax
VMin = args.VMin

if VMax is None:
    VMax = np.nanmax(V.masked_array.data[np.isfinite(V.masked_array.data)])

if VMin is None:
    VMin = np.nanmin(V.masked_array.data[np.isfinite(V.masked_array.data)])

V.caxis = [VMin,VMax]


def main():
    plt.show()