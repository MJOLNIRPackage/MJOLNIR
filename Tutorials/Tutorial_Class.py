#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 08:09:38 2019

@author: lass
"""
#import pytest

import inspect
import pytest

import sys 
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

def formatCode(text,indentChar = '>>> ', skipHeader=1):
    """Formater for raw code strings to doct string formation"""
    
    text= text.split('\n')[skipHeader:]
    
    newText = []
    figCounter = 0
    for line in text:
        if line[:8] == ' '*8: # Replace initial blank with >>>
            line = '>>> ' + line[8:]
        elif line[:4] == ' '*4:
            line = '>>> ' + line[4:]
        elif len(line)==0:
            line = '>>>  '
            
        if line.find('savefig') != -1: # code contains a save function
            startP = line.find('(')
            endP = line.find(')')
            line = line[:startP] + "('figure{}.png',format='png')".format(figCounter) + line[endP+1:]
            figCounter+=1
            
            
        newText.append(line)
            
    return '\n'.join(newText)

class Tutorial(object):
    def __init__(self,name,introText,outroText,code,fileLocation=None):
        self.name = name
        self.introText = introText
        self.outroText = outroText
        self.code = code
        if not fileLocation is None:
            if fileLocation[-1] != '/':
                fileLocation += '/'
            self.fileLocation = fileLocation + self.name.replace(' ','_')+'.rst'
        else:
            self.fileLocation = fileLocation
        
    def test(self):
        self.code()
        
    def generateTutorial(self):
        
        # Test if code is running
        codeLocation = inspect.getsourcefile(self.code)
        codeFunctionName = 'test_'+self.name.replace(' ','_')
        #print('pytest '+'{}::{}'.format(codeLocation,codeFunctionName))
        print('Running tests for "{}" in file "{}"'.format(codeFunctionName,codeLocation.split('/')[-1]))
        result = pytest.main(args=['-q','{}::{}'.format(codeLocation,codeFunctionName)])
        
        if result != 0:
            return None
            #raise RuntimeError('An error occured while running pytest for "{}" defined in function "{}"'.format(codeFunctionName,codeLocation))    
        else:
            print('Test successful!')
        fileLocation = self.fileLocation
        
        introText = self.introText
        outroText = self.outroText
        codeText = inspect.getsource(self.code)
        
        code = formatCode(codeText)
        
        
        text = introText + '\n\n' + code + '\n\n'+ outroText
        
        if fileLocation is None:
            print(text)
            
        else:
            print('Saving code example to "{}".'.format(fileLocation))
            with open(fileLocation,'w') as f:
                f.write(text)
            
