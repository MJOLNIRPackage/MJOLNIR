import sys
#import warnings
import numpy as np
from difflib import SequenceMatcher

def KwargChecker(func): # Function to check if given key-word arguments has right capitalization
    def newFunc(self,*args,**kwargs):
        def checkArgumentList(argList,kwargs):
            for key in kwargs:
                #if key.lower() in [x.lower() for x in argList] and key not in argList:
                #    expected = argList[[x.lower() for x in argList].index(key.lower())]
                    #warnings.warn(message='Did you mean "{}" in stead of "{}"?'.format(expected,key),category=SyntaxWarning,stacklevel=3)

                if key not in argList:
                    similarity = np.array([SequenceMatcher(None, key.lower(), x.lower()).ratio() for x in argList])
                    maxVal = np.max(similarity)
                    maxId = np.argmax(similarity)
                    #if(maxVal>0.8):
                    raise AttributeError('Key-word argument "{}" not understood. Did you mean "{}"?'.format(key,argList[maxId]))
                    #else:
                    #    raise AttributeError('Key-word argument "{}" not understood?'.format(key))
                    
        argList = newFunc._original.__code__.co_varnames
        checkArgumentList(argList,kwargs)
        returnval = func(self,*args,**kwargs)
        return returnval
    
    newFunc._original = func
    return newFunc
