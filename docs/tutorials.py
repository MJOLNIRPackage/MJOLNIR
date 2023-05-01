#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 13:03:12 2021

@author: Jakob Lass

Tool to replace batch in order to fully port development to windows
"""

import os,sys,shutil

DIR = os.path.join(str(os.path.dirname(__file__)),'notebooks')
RESULTDIR = os.path.join(str(os.path.dirname(__file__)),r'..\build','notebooks')

files = [f.path for f in os.scandir(DIR) if f.is_file()]

for f in files:
    if os.path.splitext(f)[1] != '.ipynb':
        continue
    sourceDate = os.path.getmtime(f)
    name = os.path.split(f)[-1]
    resultFile = os.path.join(RESULTDIR,name)
    resultDate = os.path.getmtime(resultFile)

    if sourceDate>resultDate:
        os.system('jupyter nbconvert --execute --to notebook --inplace '+f)
    else:
        print('Skipping:',name)

