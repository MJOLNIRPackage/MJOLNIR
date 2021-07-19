#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 13:03:12 2021

@author: Jakob Lass

Tool to replace batch in order to fully port development to windows
"""

import os,sys,shutil

DIR = str(os.path.dirname(__file__))

files = [f.path for f in os.scandir(DIR) if f.is_file()]

for f in files:
    if not f[-2:] == 'py':
        continue
    if 'tutorials.py' in f:
        continue
    os.system('python '+f)

for folder in ['.cache','.pyteste_cache','__pycache']:
    shutil.rmtree(os.path.join(DIR,folder))
