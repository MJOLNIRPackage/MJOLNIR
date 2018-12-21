import sys
#import warnings
import numpy as np
from difflib import SequenceMatcher
import functools
import logging
MPLKwargs = ['agg_filter','alpha','animated','antialiased or aa','clip_box','clip_on','clip_path','color or c','contains','dash_capstyle','dash_joinstyle','dashes','drawstyle','figure','fillstyle','gid','label','linestyle or ls','linewidth or lw','marker','markeredgecolor or mec','markeredgewidth or mew','markerfacecolor or mfc','markerfacecoloralt or mfcalt','markersize or ms','markevery','path_effects','picker','pickradius','rasterized','sketch_params','snap','solid_capstyle','solid_joinstyle','transform','url','visible','xdata','ydata','zorder']

def KwargChecker(function=None,include=None):
    """Function to check if given key-word is in the list of accepted Kwargs. If not directly therein, checks capitalization. If still not match raises error
    with suggestion of closest argument.
    
    Args:
    
        - func (function): Function to be decorated.

    Raises:

        - AttributeError
    """
    #@functools.wraps(None)
    def KwargCheckerNone(func):
        @functools.wraps(func)
        def newFunc(*args,**kwargs):
            argList = extractArgsList(func,newFunc,function,include)
            checkArgumentList(argList,kwargs)
            returnval = func(*args,**kwargs)
            return returnval
        newFunc._original = func
        newFunc._include = include
        newFunc._function = function
        return newFunc
    return KwargCheckerNone

def extractArgsList(func,newFunc,function,include):
    N = func.__code__.co_argcount # Number of arguments with which the function is called
    argList = list(newFunc._original.__code__.co_varnames[:N]) # List of arguments
    if not function is None:
        if isinstance(function,(list,np.ndarray)): # allow function kwarg to be list or ndarray
            for f in function:
                for arg in f.__code__.co_varnames[:f.__code__.co_argcount]: # extract all arguments from function
                    argList.append(str(arg))
        else: # if single function
            for arg in function.__code__.co_varnames[:function.__code__.co_argcount]:
                argList.append(str(arg))
    if not include is None:
        if isinstance(include,(list,np.ndarray)):
            for arg in include:
                argList.append(str(arg))
        else:
            argList.append(str(include))
        argList = list(set(argList)) # Cast to set to remove duplicates
        argList.sort() #  Sort alphabetically
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

def beautifyArgs(args=(),kwargs={}): # pragma: no cover
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

def createLogger(self,name,stream=sys.stdout,level=logging.ERROR): # pragma: no cover
    self._log =  logging.getLogger(name)
    self._log.setLevel(level)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    self._log.addHandler(ch)

def logMethod(self,original_function):  # pragma: no cover
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
    

def logClass(parent=None,log=__name__,stream=sys.stdout,level=logging.CRITICAL): # pragma: no cover
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


def binEdges(values,tolerance):
    """Generate binning of values array with minimum bin size of tolerance. Binning starts at values[0]-tolerance/2.0 and ends at values[-1]+tolerance/2.0.
    
    Args:
        
        - values (array): 1D array to be binned.
        
        - tolerance (float): Minimum length of bin sizes.
        
    Returns:
        
        - bins (array)
    
    """
    values_array = np.array(values).ravel()
    unique_values = np.asarray(list(set(values_array)))
    unique_values.sort()
    if len(unique_values)==0:
        return []
    bin_edges = [unique_values[0] - tolerance / 2.0]
    add = 1
    current = 0
    while current<len(unique_values) - 2:
        add=1
        while unique_values[current+add] - unique_values[current] < tolerance:
            if current+add < len(unique_values) - 2:
                add+=1
            else:
                current=len(unique_values)-add-1
                break
        bin_edges.append((unique_values[current] + unique_values[current+add]) / 2)
        current+=add+1
    bin_edges.append(unique_values[-1] + tolerance / 2)
    return np.array(bin_edges)


def fileListGenerator(numberString,folder,year, format = '{:}camea{:d}n{:06d}.hdf'):
    """Function to generate list of data files.
    
    Args:
        
        - numberString (str): List if numbers seperated with comma and dashes for sequences.
        
        - folder (str): Folder of wanted data files.
        
        - year (int): Year if wanted data files
        
    Kwargs:
        
        - format (str): format of data files (default '{:}camea{:d}n{:06d}.hdf')
        
    returns:
        
        - list of strings: List containing the full file string for each number provided.
        
    Example:
        >>> numberString = '201-205,207-208,210,212'
        >>> files = fileListGenerator(numberString,'data/',2018)
        ['data/camea2018n000201.hdf', 'data/camea2018n000202.hdf', 
        'data/camea2018n000203.hdf', 'data/camea2018n000204.hdf', 
        'data/camea2018n000205.hdf', 'data/camea2018n000207.hdf', 
        'data/camea2018n000208.hdf', 'data/camea2018n000210.hdf', 
        'data/camea2018n000212.hdf']
    """
        
    splits = numberString.split(',')
    dataFiles = []
    for sp in splits:
        isRange = sp.find('-')!=-1
        
        if isRange:
            spSplits = sp.split('-')
            if len(spSplits)>2:
                raise AttributeError('Sequence "{}" not understood - too many dashes.'.format(sp))
            startNumber = int(spSplits[0])
            endNumber = int(spSplits[1])
            numbers = np.arange(startNumber,endNumber+1)    
        else:
            numbers = [int(sp)]
        if not folder[-1]=='/':
            folder+='/'
        
        dataFiles.append([format.format(folder,year,x) for x in numbers])
    return list(np.concatenate(dataFiles))


def test_binEdges():
    values = np.exp(np.linspace(-0.1,1,101)) # Generate non-linear points
    #values[0]-tolerance/2.0 and ends at values[-1]+tolerance/2.0.
    minBin = 0.1
    bins = binEdges(values,minBin)

    assert(np.isclose(bins[0],values[0]-0.5*minBin)) # First bin starts at values[0]-tolerance/2.0
    assert(np.isclose(bins[-1],values[-1]+0.5*minBin)) # Last bin ends at values[-1]+tolerance/2.0
    assert(np.all(np.diff(bins)>=minBin)) # Assert that all bins are at least of size minBin

def test_fileListGenerator():
    numberStr = '0,20-23-24,4000'
    year = 2018
    folder = '/home/camea'
    try:
        files = fileListGenerator(numberStr,folder,year)
        assert False # Too many dasches
    except AttributeError:
        assert True

    numberStr = '0,20-23,4000'

    files = fileListGenerator(numberStr,folder,year)
    filesCorrect = np.array([
        '/home/camea/camea2018n000000.hdf', '/home/camea/camea2018n000020.hdf', 
        '/home/camea/camea2018n000021.hdf', '/home/camea/camea2018n000022.hdf', 
        '/home/camea/camea2018n000023.hdf', '/home/camea/camea2018n004000.hdf'])
    assert(np.all(filesCorrect == np.array(files)))
