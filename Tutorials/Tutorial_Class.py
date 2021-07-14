#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 08:09:38 2019

@author: lass
"""
#import pytest

import inspect
import pytest,types

import os,sys,numpy as np

def formatCode(text,indentChar = '   ', skipHeader=1):
    """Formater for raw code strings to doc string formation"""
    
    text= text.split('\n')[skipHeader:]
    
    newText = []
    newText.append('.. code-block:: python\n   :linenos:\n')
    figCounter = 0
    for line in text:
        if line[:8] == ' '*8: # Replace initial blank with >>>
            line = indentChar + line[8:]
        elif line[:4] == ' '*4:
            line = indentChar + line[4:]
        elif len(line)==0:
            line = indentChar
            
        if line.find('savefig') != -1: # code contains a save function
            startP = line.find('(')
            endP = line.find(')')
            line = line[:startP] + "('figure{}.png',format='png')".format(figCounter) + line[endP+1:]
            figCounter+=1
            
        if line.find('Instr = Instrument.Instrument(fileName=')!=-1:
            start = line.find('fileName=')
            end = line.find(')')
            line = line[:start] + 'fileName=SimpleInstrument.xml' + line[end:]

        while(line.find('/home/lass/Dropbox/PhD/CAMEAData/')!=-1):
            text = '/home/lass/Dropbox/PhD/CAMEAData/'
            line = line[:line.find(text)]+'/Path/To/Data/'+line[line.find(text)+len(text):]

        while(line.find('/home/lass/Dropbox/PhD/')!=-1):
            start = line.find('/home/lass/Dropbox/PhD/') # Find position in line
            end = line[start:].find("'") # Find next ' in line
            line = line[:start]+'/Path/To/Save/Folder/'+line[start+end:]
        newText.append(line)
            
    return '\n'.join(newText)

class Tutorial(object):
    def __init__(self,name,introText,outroText,code,fileLocation=None,dependentFiles = None):
        self.name = name
        self.introText = introText
        self.outroText = outroText
        self.code = code
        self.dependentFiles = dependentFiles
        if not fileLocation is None:
            if fileLocation[-1] != os.path.sep:
                fileLocation += os.path.sep
            self.fileLocation = fileLocation + self.name.replace(' ','_')+'.rst'
        else:
            self.fileLocation = fileLocation
        
    def test(self):
        if sys.platform == 'win32':
            Tester = self.code

            consts = Tester.__code__.co_consts
            tempList = np.array(list(consts),dtype=object)
            returnList = []
            for item in tempList:
                try:
                    if '/home/' in item:
                        if 'Software' in item:
                            item = item.replace('/home/lass/Dropbox/PhD/Software',r'C:/Users/lass_j/Documents/Software').replace('/','\\')
                        elif 'CAMEAData' in item:
                            item = item.replace('/home/lass/Dropbox/PhD/CAMEAData',r'C:/Users/lass_j/Documents/CAMEA2018').replace('/','\\')
                except:
                    pass
                returnList.append(item)

            Tester.__code__ = Tester.__code__.replace(co_consts=tuple(returnList))
            Tester()
        else:
            self.code()
        
    def generateTutorial(self):

        if not self.dependentFiles is None:
            
            folder = os.path.sep.join(self.fileLocation.split(os.path.sep)[:-1])+os.path.sep
            print("Copying files to {}".format(folder))
            from shutil import copyfile
            for f in list(self.dependentFiles):
                copyfile(f, folder+f.split(os.path.sep)[-1])
        
        # Test if code is running
        codeLocation = inspect.getsourcefile(self.code)
        codeFunctionName = 'test_'+self.name.replace(' ','_')
        #print('pytest '+'{}::{}'.format(codeLocation,codeFunctionName))
        print('Running tests for "{}" in file "{}"'.format(codeFunctionName,codeLocation.split(os.path.sep)[-1]))
        result = pytest.main(args=['-q','{}::{}'.format(codeLocation,codeFunctionName)])
        
        if result != 0:
            return None
            #raise RuntimeError('An error occurred while running pytest for "{}" defined in function "{}"'.format(codeFunctionName,codeLocation))    
        else:
            print('Test successful!')

        if sys.platform == 'win32':
            fileLocation = self.fileLocation.replace('/home/lass/Dropbox/PhD/Software',r'C:/Users/lass_j/Documents/Software').replace('/','\\')
        else:
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
            
