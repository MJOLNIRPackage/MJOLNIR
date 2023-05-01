# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:36:22 2023

@author: lass_j
"""


import os,sys,shutil

DIR = r'docs\notebooks'

files = [f.path for f in os.scandir(DIR) if f.is_file() and f.split('.')[-1] == 'ipynb']

for I,f in enumerate(files):
    os.system('jupyter nbconvert --execute --to notebook --inplace '+f)
    print(f'{I} of {len(files)} done')
    

# for folder in ['.cache','.pyteste_cache','__pycache']:
#     shutil.rmtree(os.path.join(DIR,folder))
