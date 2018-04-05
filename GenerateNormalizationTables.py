#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 16:49:34 2018

@author: lass
"""
from MJOLNIR.Geometry import Detector,Analyser,Wedge,Instrument
from MJOLNIR.Data import DataSet
import numpy as np
import h5py as hdf
import matplotlib.pylab as plt

Instr = Instrument.Instrument(filename='/home/lass/Dropbox/PhD/Software/CAMEA_Full.xml')
Instr.initialize()


NF = '/home/lass/Dropbox/PhD/Software/CAMEA_Test_Files/cameasim2018n000004.h5'#'/home/lass/Documents/PhD/McStas/HDFFullEi/RITA_FULL_AUTO_20171023_103342.hdf'


dataset = DataSet.DataSet(instrument=Instr,normalizationfiles=NF)

dataset.EnergyCalibration(NF,'Normalization/',plot=True
                          ,tables=['Single','PrismaticLowDefinition','PrismaticHighDefinition'])
