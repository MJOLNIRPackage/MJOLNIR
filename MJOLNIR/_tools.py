import sys
#import warnings
import numpy as np
from difflib import SequenceMatcher
import functools
import logging
import math
from MJOLNIR.Marray import *
import os
import inspect
import matplotlib

MPLKwargs = ['agg_filter','alpha','animated','antialiased or aa','clip_box','clip_on','clip_path','color or c','contains','dash_capstyle','dash_joinstyle','dashes','drawstyle','figure','fillstyle','gid','label','linestyle or ls','linewidth or lw','marker','markeredgecolor or mec','markeredgewidth or mew','markerfacecolor or mfc','markerfacecoloralt or mfcalt','markersize or ms','markevery','path_effects','picker','pickradius','rasterized','sketch_params','snap','solid_capstyle','solid_joinstyle','transform','url','visible','xdata','ydata','zorder']

#Unused
def cutObject(func): # pragma: no cover
    import MJOLNIR.Statistics.CutObject
    @functools.wraps(func)
    def newFunction(*args,**kwargs):
        if 'internal' in kwargs:
            if kwargs['internal'] is True:
                del kwargs['internal']
                return func(*args,**kwargs)
        fittingKwargs = ['fitFunction','p0','costFunction','lock']
        temp = {}
        for fitKwarg in fittingKwargs:
            if fitKwarg in kwargs:
                temp[fitKwarg] = kwargs[fitKwarg]
                del kwargs[fitKwarg]
            else:
                temp[fitKwarg] = None

        ds = args[0]
        if len(args)>1:
            otherArgs = args[1:]
        else:
            otherArgs = None
        co = MJOLNIR.Statistics.CutObject.CutObject(*func(*args,**kwargs),dataSet = ds,cutFunction = func,kwargs = kwargs,args = otherArgs)
        for fitKwarg in fittingKwargs:
            setattr(co,fitKwarg,temp[fitKwarg])
        return co
    return newFunction

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

def logAttribute(self,original_attribute): # pragma: no cover
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


def binEdges(values,tolerance,startPoint=None,endPoint=None):
    """Generate binning of values array with minimum bin size of tolerance. Binning starts at values[0]-tolerance/2.0 and ends at values[-1]+tolerance/2.0.
    
    Args:
        
        - values (array): 1D array to be binned.
        
        - tolerance (float): Minimum length of bin sizes.
        

    Kwargs:

        - startPoint (float): Minimum position from wicht to start (default None)

        - endPoint (float): Maximal end bin position (default None)

    Returns:
        
        - bins (array)
    
    """
    values_array = np.array(values).ravel().flatten()
    unique_values = np.asarray(list(set(values_array)))
    unique_values.sort()
    if len(unique_values)==0:
        return []
    bin_edges = [unique_values[0] - tolerance * 0.1]
    add = 1
    current = 0
    while current<len(unique_values) - 1:
        add=1
        broken = False
        while (unique_values[current+add]+unique_values[current+add-1])*0.5 - bin_edges[-1] < tolerance:
            if current+add < len(unique_values)-1:
                add+=1
            else:
                broken=True
                break
        if not broken:
            bin_edges.append((unique_values[current+add-1] + unique_values[current+add]) / 2)
        current+=add
    if unique_values[-1]-bin_edges[-1]< 1.1*tolerance:
        bin_edges.append(bin_edges[-1]+tolerance)
    else:
        bin_edges.append(unique_values[-1]+0.1*tolerance)
    
    bin_edges = np.array(bin_edges)

    if not endPoint is None:
        if endPoint-bin_edges[-1]<tolerance:
            bin_edges = np.concatenate([bin_edges[:np.sum(bin_edges<endPoint)],[endPoint]])
    if not startPoint is None:
        if bin_edges[0]-startPoint<tolerance or bin_edges[0]<startPoint:
            bin_edges = np.concatenate([[startPoint],bin_edges[np.sum(bin_edges<startPoint):]])
    return bin_edges


def without_keys(dictionary, keys): # Remove key word argument from kwargs
    return {x: dictionary[x] for x in dictionary if x not in keys}

@KwargChecker()
def fileListGenerator(numberString,folder,year=2018, format = None, instrument = 'CAMEA'):
    """Function to generate list of data files.
    
    Args:
        
        - numberString (str): List if numbers separated with comma and dashes for sequences.
        
        - folder (str): Folder of wanted data files.
        
    Kwargs:

        - year (int): Year of wanted data files (default 2018)

        - format (str): format of data files (default None, but CAMEA if instrument is provided)

        - instrument (str): Instrument to be used to determine format string (default CAMEA)
        
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
    if format is None: # If no user specified format is provided
        if instrument == 'CAMEA':
            format = 'camea{:d}n{:06d}.hdf'
        elif instrument == 'MultiFLEXX':
            format = '{1:06d}'
        else:
            raise AttributeError('Provided instrument "{}" not understood'.format(instrument))


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

        dataFiles.append([os.path.join(folder,format.format(year,x)) for x in numbers])
    return list(np.concatenate(dataFiles))




def RoundBinning(X,deltas,Data=None):
    """ Bin points to nearest delta value in D dimensions and reorder Data.
    
    Args:
        
         - X (np.array): Input points in shape (D,N) for D dimensions.
        
         - deltas (float or list): binning size(s)
        
    Kwargs:
        
         - Data (list/array): List of data or data list of shape (N) to be binned like X (default None).
        
    Returns:
        
         - BinnedPoints (list): ND list of binned unique positions of points
        
         - indices (list): Inverse indices from which points-array can be created from unique points
        
         - count (list): Number of points going into each bin
        
        (- returnData: Rebinned data as according to points being binned) 
        
    Algorithm takes the data points and rebins them into points closest to delta in N dimensions. If deltas=[0.05,0.05] and 2D points are given, points will be binned to closest 0.05 (...-0.1,-0.05,0.0,0.05,0.1...) in both directions. Data lists are also reshuffeled to match.
    
    
    """
    points = np.array(X)
    deltas = np.array([deltas]).squeeze()
    if len(deltas.shape)==0 and points.shape[0]!=1:
        binning = np.ones([points.shape[0],1])*deltas
    else:
        binning = np.array([deltas]).reshape(-1,1)
        if not points.shape[0]==binning.shape[0]:
            raise AttributeError('Shape mismatch between X and deltas. Expected {} deltas for X with shape {}, recieved {}'.format(points.shape[0],points.shape,deltas))
        
    binningLog = np.log10(binning)
    binningOrder = np.floor(binningLog)
    binningScale = np.power(10,(binningOrder-binningLog)).reshape(-1,1)
    
    Scale = np.multiply(points,binningScale)
    
    binPoints,indices,count = np.unique(np.round(Scale*np.power(10,-binningOrder)),axis=1,return_inverse=True,return_counts=True)
    
    BinnedPoints = binPoints*np.power(10,binningOrder)/binningScale
    
    if not Data is None:
        if isinstance(Data,list):
            returnData = [np.histogram(indices,bins=BinnedPoints.shape[1],weights=x)[0] for x in Data]
        else:
            returnData = np.histogram(indices,bins=BinnedPoints.shape[1],weights=Data)[0]
        return BinnedPoints,indices,count,returnData
    
    return BinnedPoints,indices,count
    

def invert(M):
    """Invert non-square matrices as described on https://en.wikipedia.org/wiki/Generalized_inverse.
    
    Args:
        
        - M (matrix): Matrix in question.
        
    Returns:
        
        - Left or right inverse matrix depending on shape of provided matrix.
    """
    s = M.shape
    if s[0]>s[1]:
        return np.dot(np.linalg.inv(np.dot(M.T,M)),M.T)
    else:
        return np.dot(M.T,np.linalg.inv(np.dot(M,M.T)))


def overWritingFunctionDecorator(overWritingFunction):
    def overWriter(func):
        return overWritingFunction
    return overWriter



def calRefVector(points):# pragma: no cover
    """ Calcualte reference vector as vector pointing from mean point to geometric center. For half moon shape this is anti-radially.
    
    Args:
        
        - points (list): list of points for which reference vector is calcualted, shape is 2xN
        
    Returns:
        
        vector: Reference vector
    
    """
    center = np.mean(points,axis=1).reshape(2,1)
    argMinDist = np.argmin(np.linalg.norm(center-points,axis=0))
    return center-points[:,argMinDist].reshape(2,1),center


def minMax(x,axis=None):
    """Return minimal and maximal of list.
    
    Args:
        
        - x (list): Object from which min and max is to be found.
        
    Kwargs:
        
        - axis (int): Axis or axes along which to operate (default 0)
      
    Returns:
        
        - min: Minimal value
        
        - max: Maximal value
        
    """
    return np.min(x,axis=axis),np.max(x,axis=axis)


def unitVector(v):
    """Returns vector of unit length"""
    return v/np.linalg.norm(v)

def rotate2X(v):
    if np.isclose(v[2]/np.linalg.norm(v),1): # v is along z
        return rotMatrix([0,1,0],np.pi/2,deg=False)
    # Find axis perp to v and proj v into x-y plane -> rotate 2 plane and then to x
    vRotInPlane = np.array([-v[1],v[0],0])
    vPlan = np.array([v[0],v[1],0])
    ## TODO: Check this!
    #if np.isclose(np.dot(v,vPlan)/(np.linalg.norm(v)*np.linalg.norm(vPlan)),1.0):
    #    return rotMatrix([1,0,0],0.0)
    theta = np.arccos(np.dot(v,vPlan)/(np.linalg.norm(v)*np.linalg.norm(vPlan)))
    R = rotMatrix(vRotInPlane,theta,deg=False)
    v2 = np.dot(R,v)
    theta2 = np.arccos(np.dot(v2,np.array([1,0,0]))/np.linalg.norm(v2))
    R2 = rotMatrix(np.array([0,0,1.0]),-theta2,deg=False)
    
    Rotation = np.dot(R2,R)
    return Rotation

def LengthOrder(v):
    
    nonZeroPos = np.logical_not(np.isclose(v,0.0))
    if np.sum(nonZeroPos)==1:
        Rv = v/np.linalg.norm(v)
        return Rv
    if np.sum(nonZeroPos)==0:
        raise AttributeError('Provided vector is zero vector!')
    
    if np.sum(nonZeroPos)==3:
        v1 = Norm2D(v[:2])
        ratio = v1[0]/v[0]
        v2 = Norm2D(np.array([v1[0],v[2]*ratio]))
        ratio2 = v2[0]/v1[0]
        Rv = np.array([v2[0],v1[1]*ratio2,v2[1]])
    else:
        Rv = np.zeros(3)
        nonZeros = v[nonZeroPos]
        Rv[nonZeroPos] = Norm2D(nonZeros)
    
    if not np.isclose(np.dot(Rv,v)/(np.linalg.norm(Rv)*np.linalg.norm(v)),1.0):
        raise AttributeError('The found vector is not parallel to original vector: {}, {}',format(Rv,v))
    return Rv



def Norm2D(v):
    reciprocal = np.abs(1/v)
    if np.isclose(reciprocal[0],reciprocal[1]):
        return v*reciprocal[0]
    
    ratio = np.max(reciprocal)/np.min(reciprocal)
    if np.isclose(np.mod(ratio,1),0.0) or np.isclose(np.mod(ratio,1),1.0):
        return v*np.min(reciprocal)*ratio
    else:
        return v


def rotMatrix(v,theta,deg=True):
    """ Generalized rotation matrix.
    
    Args:
        
        - v (list): Rotation axis around which matrix rotates
        
        - theta (float): Rotation angle (by default in degrees)
        
    Kwargs:
        
        - deg (bool): Whether or not angle is in degrees or radians (Default True)
        
    Returns:
        
        - 3x3 matrix rotating points around vector v by amount theta.
    """
    if deg==True:
        theta = np.deg2rad(theta.copy())
    v/=np.linalg.norm(v)
    m11 = np.cos(theta)+v[0]**2*(1-np.cos(theta))
    m12 = v[0]*v[1]*(1-np.cos(theta))-v[2]*np.sin(theta)
    m13 = v[0]*v[2]*(1-np.cos(theta))+v[1]*np.sin(theta)
    m21 = v[0]*v[1]*(1-np.cos(theta))+v[2]*np.sin(theta)
    m22 = np.cos(theta)+v[1]**2*(1-np.cos(theta))
    m23 = v[1]*v[2]*(1-np.cos(theta))-v[0]*np.sin(theta)
    m31 = v[0]*v[2]*(1-np.cos(theta))-v[1]*np.sin(theta)
    m32 = v[1]*v[2]*(1-np.cos(theta))+v[0]*np.sin(theta)
    m33 = np.cos(theta)+v[2]**2*(1-np.cos(theta))
    return np.array([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

def Rot(theta,deg=True):
    """Create 2D rotation matrix
    
    Args:
        
        - theta (float): Rotation angle (by default in degrees)
        
    Kwargs:
        
        - deg (bool): Whether or not number provided is degree or radian (default True)
      
    Returns:
        
        - 2x2 rotation matrix definde with -sin in row 0 column 1
        
    """
    if deg==True:
        theta = np.deg2rad(theta)
    return np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])



def clockwiseangle_and_distance(point,origin=[0,0],refvec = [0,1]): # pragma: no cover
    """Sort points clockwise. Taken from https://stackoverflow.com/questions/41855695/sorting-list-of-two-dimensional-coordinates-by-clockwise-angle-using-python
    
    Args:
        
        - point (list): List of points in 2D of size 2xN
        
    Kwargs:
        
        - origin (list): Location of origin from which the points are to be sorted (default [0,0])
        
        - refvec (list): Vector direction for definition of zero point (default [0,1])
        
    """
    # Vector between point and the origin: v = p - o
    vector = [point[0]-origin[0], point[1]-origin[1]]
    # Length of vector: ||v||
    lenvector = math.hypot(vector[0], vector[1])
    # If length is zero there is no angle
    if lenvector == 0:
        return -math.pi, 0
    # Normalize vector: v/||v||
    normalized = [vector[0]/lenvector, vector[1]/lenvector]
    dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
    diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
    angle = math.atan2(diffprod, dotprod)
    # Negative angles represent counter-clockwise angles so we need to subtract them 
    # from 2*pi (360 degrees)
    if angle < 0:
        return 2*math.pi+angle, lenvector
    # I return first the angle because that's the primary sorting criterium
    # but if two vectors have the same angle then the shorter distance should come first.
    return angle, lenvector



@KwargChecker()
def rotationMatrix(alpha,beta,gamma,format='deg'):
    if format=='deg':
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    Rx = np.array([[1,0,0],[0,np.cos(alpha),-np.sin(alpha)],[0,np.sin(alpha),np.cos(alpha)]])
    Ry = np.array([[np.cos(beta),0,np.sin(beta)],[0,1,0],[-np.sin(beta),0,np.cos(beta)]])
    Rz = np.array([[np.cos(gamma),-np.sin(gamma),0],[np.sin(gamma),np.cos(gamma),0],[0,0,1]])
    return np.dot(Rz,np.dot(Ry,Rx))

@KwargChecker()
def vectorAngle(V1,V2):
    """calculate angle between V1 and V2.
    
    Args:
    
        - V1 (list): List or array of numbers
        
        - V2 (list): List or array of numbers
        
    Return:
        
        - theta (float): Angle in degrees between the two vectors
    """
    return np.arccos(np.dot(V1,V2.T)/(np.linalg.norm(V1)*np.linalg.norm(V2)))


# TODO: Write a test/update this function
def writeToSpinWFile(file,position,spinWaveEnergy,spinWaveWidth,spinWaveAmplitude,EMin,EMax,spinWaveEnergyErr=None): # pragma: no cover
    """Write fitted values for spin wave(s) into a SpinW readable format.
    
    Args:
        
        - files (string): File into which the spin waves is to be saved
        
        - position (3D vector, 3 x m): HKL position of spin wave(s)
        
        - spinWaveEnergy (array, n x m): Array with energy position of spin wave. For multiple spin waves fill with 0 if wave not found
        
        - spinWaveWidth (array, n x m): Standard deviation of spin wave(s). For multiple spin waves fill with 0 if wave not found
        
        - spinWaveAmplitude (array, n x m): Amplitude of spin wave(s). For multiple spin waves fill with 0 if wave not found
        
        - EMin (float): Lowest energy measured in data [meV]
        
        - EMin (float): Highest energy measured in data [meV]
        
        
    n is the number of spin waves
    m is the number of data points measured        
    """
    spinWaveEnergy = np.asarray(spinWaveEnergy)
    if not spinWaveEnergyErr is None:
        spinWaveEnergyErr = np.asarray(spinWaveEnergyErr)
        if not spinWaveEnergy.shape == spinWaveEnergyErr.shape:
            raise AttributeError('Arrays for spinWaveEnergy(shape: {}), spinWaveEnergyErr(shape: {}) have to have same shape.'.format(spinWaveEnergy.shape,spinWaveEnergyErr.shape))
   
        
    spinWaveWidth = np.asarray(spinWaveWidth)
    spinWaveAmplitude = np.asarray(spinWaveAmplitude)
    position = np.asarray(position)
    
    if not np.all([spinWaveEnergy.shape == spinWaveWidth.shape,spinWaveEnergy.shape == spinWaveAmplitude.shape]):
        raise AttributeError('Arrays for spinWaveEnergy(shape: {}), spinWaveWidth(shape: {}), and spinWaveAmplitude(shape: {}) have to have same shape.'.format(spinWaveEnergy.shape,spinWaveWidth.shape,spinWaveAmplitude.shape))
        
    if len(spinWaveEnergy.shape) == 1:
        spinWaveEnergy.shape = (1,-1)
        spinWaveWidth.shape = (1,-1)
        spinWaveAmplitude.shape = (1,-1)
        
    spinWaves,dataPoints = spinWaveEnergy.shape
    
    if not position.shape == (3,dataPoints):
        raise AttributeError('With provided spin wave parameters, expected position of shape (3,{}) but recieved {}.'.format(dataPoints,position.shape))
    
    
    if not spinWaveEnergyErr is None:
        columns = 4
    else:
        columns = 3
    dataMatrix = np.zeros((5+columns*spinWaves,dataPoints))
    
    dataMatrix[:3,:] = position
    dataMatrix[3,:] = EMin
    dataMatrix[4,:] = EMax
    
    energyIndices = np.array([6+x*columns for x in range(spinWaves)],dtype=int)
    widthIndices = np.array([7+x*columns for x in range(spinWaves)],dtype=int)
    amplitudeIndices = np.array([5+x*columns for x in range(spinWaves)],dtype=int)
    
    if not spinWaveEnergyErr is None:
        errorIndices = np.array([8+x*columns for x in range(spinWaves)],dtype=int)
        dataMatrix[errorIndices,:] = spinWaveEnergyErr
    
    dataMatrix[energyIndices,:] = spinWaveEnergy
    dataMatrix[widthIndices,:] = spinWaveWidth
    dataMatrix[amplitudeIndices,:] = spinWaveAmplitude
    
    headers = ['QH','QK','QL','ENlim1','ENlim2']
    if not spinWaveEnergyErr is None:
        headers+=list(np.concatenate([['I{}'.format(x+1),'EN{}'.format(x+1),'sigma{}'.format(x+1),'EErr{}'.format(x+1)] for x in range(spinWaves)]))
    else:
        headers+=list(np.concatenate([['I{}'.format(x+1),'EN{}'.format(x+1),'sigma{}'.format(x+1)] for x in range(spinWaves)]))
    titles = ''.join(['{:>14s} '.format(x) for x in headers])

    
    np.savetxt(file,dataMatrix.T,fmt='%+.7e',delimiter = ' ', header=titles,comments='')



def generateLabel(vec,labels=['H','K','L']):
    """Format a scattering vector with individual letters.
    
    Args:
    
        - vec (array): Vector to be formated.

    Kwargs:
        
        - lables (list): Letters to use for formating (default ['H','K','L']) 
        
    """
    # Individual letters
    integers = np.isclose(np.abs(vec)-np.floor(np.abs(vec)),0.0,atol=1e-4)
    signs = np.sign(vec)
    vec = np.abs(vec)
    zeros = np.isclose(vec,0.0,atol=1e-4)
    ones = np.isclose(vec,1.0,atol=1e-4)
    
    label = []
    for l,i,s,z,o,v in zip(labels,integers,signs,zeros,ones,vec):
        sign= '-' if s==-1 else ''
        if not i: # not an integer
            label.append(''.join([sign,str(v),l]))
        else: # It is an integer
            if z: # if zero
                label.append('0')
            elif o: # if one
                label.append(''.join([sign,l]))
            else: # if integer different from 0 or +- 1
                label.append(''.join([sign,str(int(v)),l]))
                
    returnLabel = ''.join(['(',', '.join([x for x in label]),')'])
    return returnLabel
