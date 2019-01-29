#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 10:31:10 2019

@author: lass
"""

import sys

string  = sys.argv[1]
version = [int(x) for x in string.split('-')[1].split('.')[:-2]]
version[-1]+=1
print('.'.join([str(x) for x in version]))

