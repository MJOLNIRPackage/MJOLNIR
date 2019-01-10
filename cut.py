#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 10:31:10 2019

@author: lass
"""

import sys

string  = sys.argv[1]

print('.'.join([x for x in string.split('-')[1].split('.')[:-2]]))