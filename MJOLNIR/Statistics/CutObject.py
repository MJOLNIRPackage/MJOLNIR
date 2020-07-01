# -*- coding: utf-8 -*-
import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')

import numpy as np
import pickle as pickle
import matplotlib.pyplot as plt
import matplotlib.axes
from matplotlib.collections import PatchCollection,PolyCollection
from MJOLNIR.Data import Viewer3D,RLUAxes
import MJOLNIR.Data.DataSet
import MJOLNIR.Data.Sample
from MJOLNIR import _tools
import pytest
from scipy.ndimage import filters
import scipy.optimize
from scipy.spatial import Voronoi,ConvexHull,KDTree
import time
import warnings

pythonVersion = sys.version_info[0]



class CutObject(object):# pragma: no test
    def __init__(self,*returnValues, **kwargsCut): # dataSet = None, args = None, kwargs = None, cutFunction = None
        dataSet = kwargsCut.get('dataSet',None)
        args = kwargsCut.get('args',None)
        kwargs = kwargsCut.get('kwargs',None)
        cutFunction = kwargsCut.get('cutFunction',None)
        if isinstance(dataSet,MJOLNIR.Data.DataSet.DataSet):
            self.dataSet = dataSet
        elif not dataSet is None:
            raise AttributeError('dataSet attribute has to be of type MJOLNIR.Data.DataSet.')

        if issubclass(type(returnValues[0]),matplotlib.axes.Axes):
            self._ax,data = returnValues[0],returnValues[1:]
        else:
            data = returnValues
        
        if len(data)>2:
            #print(self,'HERE')
            dataList,binList,centerPosition,binDistance = data
            self.x = centerPosition
            #print(len(dataList))
        else:
            #print(self,'Not there')
            dataList,binList = data
            centerPosition,binDistance = None,None
            X = np.array(binList[0])
            self.x = 0.5*(X[1:]+X[:-1])
        self._data = []
      

        # Extract fitting related Kwargs
        
        self.parameters = ParameterDictionary()
        if isinstance(cutFunction,type(lambda:None)):
            if hasattr(cutFunction,'_original'):
                self.cutFunction = cutFunction._original
            else:
                self.cutFunction = cutFunction

            if not args is None:
                for param,value in zip(self.cutFunction.__code__.co_varnames[1:],args):
                    self.parameters[param] = value

            self.parameters.update(kwargs)
            if hasattr(self,'_ax'):
                self.parameters.update({'ax':self._ax})
            self.__parameters = self.parameters.copy()
        elif not cutFunction is None:
            raise AttributeError('cutFunction has to be of type function.')
        
        
        self.dataList = dataList
        self.binList = binList
        functionData = [self.dataList,self.binList]
        

        if not centerPosition is None:
            functionData.append(centerPosition)
            self.centerPosition = centerPosition
        if not binDistance is None:
            functionData.append(binDistance)
            self.binDistance = binDistance
        self._data = functionData
        
        

        if self.data is None:
            self.recalculate = True # Flag to check if anything has changed since last calculation of data
        else:
            self.recalculate = False

    @property
    def data(self):
        return self._data
    
    @data.getter
    def data(self):
        if hasattr(self.dataSet,'isChanged'):
            self.recalculate |= self.dataSet.isChanged
        if not hasattr(self,'recalculate'):
            self.recalculate = False
        if self.recalculate or self.parameters._updated:
            self.recalculateData()
            self.recalculate = False

        return self._data

        
    @data.setter
    def data(self,data):
        self._data = data
        self.dataList = self._data[0]
        self.binList = self._data[1]
        if len(data)>2:
            self.centerPosition = self._data[2]
            if len(data)>3:
                self.binDistance = self._data[3]
        

    @property
    def I(self):
        return self._I
    
    @I.getter
    def I(self):
        if hasattr(self,'dataList'):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if len(self.dataList) == 1:
                    return np.divide(self.dataList[0][0]*self.dataList[0][-1], self.dataList[0][1]*self.dataList[0][2]).flatten()
                else:
                    return np.divide(self.dataList[0]*self.dataList[-1], self.dataList[1]*self.dataList[2]).flatten()
        else:
            raise AttributeError('CutObject does not have data...')

    @property
    def fitFunction(self):
        return self._fitFunction

    @fitFunction.getter
    def fitFunction(self):
        return self._fitFunction
    
    @fitFunction.setter
    def fitFunction(self,function):
        self._fitFunction = function
        
        """if hasattr(self,'binCenter'):
            X = self.binCenter
            print(np.array(X).shape)
        elif hasattr(self,'binList'):
            pass
        else:
            raise AttributeError('CutObject without binCenter or binList cannot be used to fit...')"""

    @property
    def lock(self):
        return self._lock

    @lock.getter
    def lock(self):
        return self._lock
    
    @lock.setter
    def lock(self,lock):
        if lock is None and hasattr(self,'p0'):
            if not self.p0 is None:
                self._lock = np.zeros((len(self.p0)),dtype=bool)
        self._lock = lock


    @property
    def p0(self):
        return self._p0

    @p0.getter
    def p0(self):
        return self._p0
    
    @p0.setter
    def p0(self,p0):
        
        self._p0= p0


    def recalculateData(self):
        if hasattr(self,'_ax'):
            if hasattr(self._ax,'colorbar'):
                self._ax.colorbar.remove()
            self._ax.cla()
            if hasattr(self._ax,'_button_press_event'):
                self._ax.figure.canvas.mpl_disconnect(self._ax._button_press_event)
            returnValues = self.cutFunction(self=self.dataSet,**self.parameters._dict)
            self.ax = returnValues[0]
            self.data = returnValues[1:]
        else:
            self.data = self.cutFunction(self=self.dataSet,**self.parameters._dict)
        self.recalculate = False
        self.parameters._updated = False
        self.__parameters = self.parameters.copy()
        


    def __getitem__(self,index):
        try:
            if hasattr(self,'_ax'):
                if index == 0:
                    return self._ax
                else:
                    self.data[index]
            else:
                return self.data[index]
        except IndexError:
            raise IndexError('Provided index {} is out of bounds for DataSet with length {}.'.format(index,len(self)))

    def __len__(self):
        return len(self.data)
    
    def __iter__(self):
        self._index=0
        return self
    
    def __next__(self):
        
        if hasattr(self,'_ax'):
            if self._index >= len(self)+1:
                raise StopIteration
            if self._index == 0:
                result = self._ax
            else:
                result = self.data[self._index-1]
        else:
            if self._index >= len(self):
                raise StopIteration
            result = self.data[self._index]
        self._index += 1
        return result

    def next(self):
        return self.__next__()

    def __delitem__(self,index):
        raise NotImplementedError('It is not possible to delete items in a CutObject.')


    def fittingIterator(self):
        X = np.concatenate(self.centerPosition,axis=0)
        IDs = np.concatenate([[i]*len(CP) for i,CP in enumerate(self.centerPosition)])

        DATA = []
        for D in self.dataList:
            DataLocal = []
            for E in range(len(D[0])):
                DataLocal.append([D[0][E],D[1][E],D[2][E],D[3][E]])
            DATA.append(DataLocal)


        DATA = np.concatenate(DATA,axis=0)
        QDistance = np.concatenate(self._ax.edgeQDistance)
        return zip(IDs,X,DATA,QDistance)
    


class ParameterDictionary(dict): # pragma: no test
    def __init__(self,*args,**kwargs):
        self._dict = dict(*args,**kwargs)
        self._updated =False
    def __getitem__(self, key):
        if key not in self._dict.keys():
            raise KeyError

        return self._dict[key]

    def __setitem__(self, key, value):
        if key not in self._dict.keys():
            self._updated = True
            self._dict[key] = value
        else:
            try:
                equal = self[key] != value
                equal = np.all(np.isclose(self._dict[key],value))
            except:
                equal = False
            if not equal:
                self._dict[key] = value
                self._updated = True
    
    def __repr__(self):
        return repr(self._dict)

    def __len__(self):
        return len(self._dict)

    def __delitem__(self, key):
        del self._dict[key]

    def clear(self):
        return self._dict.clear()

    def copy(self):
        return self._dict.copy()

    def has_key(self, k):
        return k in self._dict

    def update(self, *args, **kwargs):
        return self._dict.update(*args, **kwargs)

    def keys(self):
        return self._dict.keys()

    def values(self):
        return self._dict.values()

    def items(self):
        return self._dict.items()

    def pop(self, *args):
        return self._dict.pop(*args)

    def __cmp__(self, dict_):
        return self.__cmp__(self._dict, dict_)

    def __contains__(self, item):
        return item in self._dict

    def __iter__(self):
        return iter(self._dict)

    def __unicode__(self):
        return unicode(repr(self._dict))
        

def test_ParameterDictionary_initialization():
    PD = ParameterDictionary()
    assert(len(PD)==0)
    pars = {'a':10,'b':[1,2,3]}
    PD = ParameterDictionary(pars)
    assert(pars==pars)

    assert(PD['a']==10)
    print(PD)
    print(str(PD))
    print(repr(PD))

    for key,value in PD.items():
        print(key,value)
    assert(PD.has_key('b'))

    assert(np.all([x in PD.keys() for x in ['a','b']]))
    assert(len(PD)==2)

    