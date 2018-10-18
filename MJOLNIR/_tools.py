import sys
#import warnings
import numpy as np
from difflib import SequenceMatcher
import functools
import logging

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
            argList = extractArgsList(func,newFunc,function,include)
            checkArgumentList(argList,kwargs)
            returnval = func(self,*args,**kwargs)
            return returnval
        #newFunc.__dict__ =  func.__dict__
        newFunc._original = func
        newFunc._include = include
        newFunc._function = function
        
        
        return newFunc
    return KwargCheckerNone

def extractArgsList(func,newFunc,function,include):
    N = func.__code__.co_argcount            
    argList = list(newFunc._original.__code__.co_varnames[:N])
    if not function is None:
        for arg in function.__code__.co_varnames[:function.__code__.co_argcount]:
            argList.append(str(arg))
    if not include is None:
        for arg in include:
            argList.append(str(arg))
    return argList

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
        error = AttributeError(errorMsg)
        raise error

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

def beautifyArgs(args=(),kwargs={}):
    """Beautify arguments and keyword arguments. Returns formated string with arguments and 
    keyword argumenets seperated with commas as called in a function"""
    returnStr = ''
    if not args == ():
        args = list(args)
        returnStr =', '.join([str(x) if type(x)!=str else '"{}"'.format(x) for x in args])
        if not kwargs == {}:
            returnStr+=', '
    if not kwargs == {}:
        kwformat = ['{}={}'.format(key,kwargs[key]) if type(kwargs[key])!=str else '{}="{}"'.format(key,kwargs[key])  for key in kwargs]
        returnStr+=', '.join([str(x) for x in kwformat])
    return returnStr

def createLogger(self,name,stream=sys.stdout,level=logging.ERROR):
    self._log =  logging.getLogger(name)
    self._log.setLevel(level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    self._log.addHandler(ch)

def logMethod(self,original_function): 
    @functools.wraps(original_function)     
    def new_function(*args,**kwargs):
        self._log.info('Calling {}({})'.format(new_function.original_function.__name__,beautifyArgs(args,kwargs)))
        try:
            x = original_function(*args,**kwargs)
        except Exception as WrapperError:
            self._log.error('***** FAILED *****')
            raise WrapperError
        return x    
    new_function.original_function= original_function                                   
    return new_function  

def logAttribute(self,original_attribute):      
    self._log.info('Calling attribute {}'.format(str(original_attribute)))                                     
    

def logClass(parent=None,log=__name__,stream=sys.stdout,level=logging.CRITICAL):
    if parent is None:
        parent = object
    def track_all_class_methods(Cls):
        class NewCls(parent):
            
            def __init__(self,*args,**kwargs):
                createLogger(self,name=log,stream=stream,level=level)
                self._log.info('Calling {}({})'.format(Cls.__name__,beautifyArgs(args,kwargs)))
                self._Instance = Cls(*args,**kwargs)

                
            def __getattribute__(self,s):
                """
                this is called whenever any attribute of a NewCls object is accessed. This function first tries to 
                get the attribute off NewCls. If it fails then it tries to fetch the attribute from self.oInstance (an
                instance of the decorated class). If it manages to fetch the attribute from self.oInstance, and 
                the attribute is an instance method then `time_this` is applied.
                """
                try:    
                    x = super(NewCls,self).__getattribute__(s)
                    skip = True
                except AttributeError:      
                    x = self._Instance.__getattribute__(s)
                    skip = False
                else:
                    pass
                finally:
                    if skip and s[0]!='_':
                        if type(x) == type(self.__init__):# and x != self.__init__: # it is an instance method
                            return logMethod(self,x)                 # this is equivalent of just decorating the method with logMethod
                            pass
                        else:
                            logAttribute(self,s)
                    return x

            def _printLog(self):
                pass
                #print('\n'.join([x for x in self._log]))
        d = Cls.__dict__
        keys = d.keys()
        copy = [x for x in list(Cls.__dict__) if not x in list(NewCls.__dict__)]
        
        for key in copy+['__doc__']:
            if key=='__dict__' or key in ['__getattribute__','__init__']: #  Check if the value in dict is a function
                continue
            
            if type(d[key])==type(Cls.__init__):
                if hasattr(d[key],'_original'):

                    setattr(NewCls, key, KwargChecker(d[key]._function,d[key]._include)(d[key]._original))
                else:
                    setattr(NewCls, key, d[key])
                    for funckey in d[key].__dict__:
                        setattr(d[key],funckey,d[key].__dict__[funckey])
                        
            else:
                setattr(NewCls, key, d[key])
        return NewCls
    return track_all_class_methods
