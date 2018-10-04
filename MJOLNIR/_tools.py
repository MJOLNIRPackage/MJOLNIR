import sys
#import warnings
import numpy as np
from difflib import SequenceMatcher

def KwargChecker(func):
    """Function to check if given key-word is in the list of accepted Kwargs. If not directly therein, checks capitalization. If still not match raises error
    with suggestion of closest argument.
    
    Args:
    
        - func (function): Function to be decorated.

    Raises:

        - AttributeError
    """
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

    for attr in ['__module__', '__name__', '__qualname__', '__doc__', '__annotations__','__weakref__']:
        try:
            value = getattr(func, attr)
        except AttributeError:
            pass
        else:
            setattr(newFunc, attr, value)
    for attr in ['__dict__']:
        getattr(newFunc, attr).update(getattr(func, attr, {}))
    newFunc._original = func
    return newFunc


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