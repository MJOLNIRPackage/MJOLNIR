import sys
#import warnings
import numpy as np
from difflib import SequenceMatcher
import functools

def KwargChecker(function=None,include=None):
    """Function to check if given key-word is in the list of accepted Kwargs. If not directly therein, checks capitalization. If still not match raises error
    with suggestion of closest argument.
    
    Args:
    
        - func (function): Function to be decorated.

    Raises:

        - AttributeError
    """
    
    def KwargCheckerNone(func):
        @functools.wraps(func)
        def newFunc(self,*args,**kwargs):
            N = func.__code__.co_argcount            
            argList = list(newFunc._original.__code__.co_varnames[:N])
            if not function is None:
                for arg in function.__code__.co_varnames[:function.__code__.co_argcount]:
                    argList.append(str(arg))
            if not include is None:
                for arg in include:
                    argList.append(str(arg))
            checkArgumentList(argList,kwargs)
            returnval = func(self,*args,**kwargs)
            return returnval
        newFunc._original = func
        return newFunc
    return KwargCheckerNone

def checkArgumentList(argList,kwargs):
    notFound = []
    for key in kwargs:
        if key not in argList:
            similarity = np.array([SequenceMatcher(None, key.lower(), x.lower()).ratio() for x in argList])
            maxVal = np.max(similarity)
            maxId = np.argmax(similarity)
            notFound.append('Key-word argument "{}" not understood. Did you mean "{}"?'.format(key,argList[maxId]))
    if len(notFound)>0:
        if len(notFound)>1:
            errorMsg = 'The following key-word arguments are not understood:\n'
            errorMsg+='\n'.join(notFound)
        else:
            errorMsg = notFound[0]
        raise AttributeError(errorMsg)

def my_timer_N(N=0): # pragma: no cover
    """Timer function to measure time consumbtion of function.

    Kwargs:

        - N (int): Number of itterations to perform.

    Raises:

        - AttributeError
    """
    if N<0:
        raise AttributeError('Number of runs need to be bigger or equal to 1 or equal to 0 for no timing, but {} given.'.format(N))
    def my_timer(func):
        import time
        def newFunc(*args,**kwargs):
            Time = []
            if N ==0:
                returnval = func(*args,**kwargs)
            else:
                for i in range(N):
                    startT = time.time()
                    returnval = func(*args,**kwargs)
                    stopT = time.time()
                    Time.append(stopT-startT)
                if N>1:
                    print('Function "{}" took: {}s (\\pm{}s)'.format(func.__name__,np.mean(Time),np.std(Time)/np.sqrt(N)))
                else:
                    print('Function "{}" took: {}s'.format(func.__name__,Time[0]))
            return returnval
        return newFunc
    return my_timer