import copy
import functools
import os
import sys
from operator import sub

sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from MJOLNIR import _tools
from mpl_toolkits.axisartist import SubplotHost
from mpl_toolkits.axisartist.grid_helper_curvelinear import \
    GridHelperCurveLinear





def createRLUAxes(self,figure=None,ids=[1, 1, 1],nbinsx=None,nbinsy=None,basex=None,basey=None):
    """Create a reciprocal lattice plot for a given DataSet object.
    
    Args:
        
        - Dataset (DataSet): DataSet object for which the RLU plot is to be made.

    Kwargs:

        - figure: Matplotlib figure in which the axis is to be put (default None)

        - ids (array): List of integer numbers provided to the SubplotHost ids attribute (default [1,1,1])
        
        - nbinsx (int): Number of bins used along the x-axis (default around 5 - depends on angle between projection vectors)

        - nbinsy (int): Number of bins used along the y-axis (default None)

        - basex (float): Ticks are positioned at multiples of this value along x (default None)

        - basey (float): Ticks are positioned at multiples of this value along y (default None)

    Returns:
        
        - ax (Matplotlib axes): Axes containing the RLU plot.

    .. note::
        When rlu axis is created, the orientation of Qx and Qy is assumed to be rotated as well. 
        This is to be done in the self.View3D method call!

    .. note::
        The number of ticks and their location cannot be changed after the initiation of the axis due to current stage of experimental development
        of the GridHelperCurveLinear object. Provide either nbinsy, (nbinsx,nbinsy), basex, or (basex,basey). If none is provided 0.25 is used as base 
        in both directions.
        
    """

    sample = copy.deepcopy(self.sample)
    
    sample.convert = np.einsum('ij,j...->i...',sample.RotMat,sample.convert)
    sample.convertinv = np.linalg.inv(sample.convert) # Convert from Qx, Qy to projX, projY

    sample.convertHKL = np.einsum('ij,j...->i...',sample.RotMat,sample.convert)
    sample.convertHKLINV = _tools.invert(sample.convertHKL) # Convert from Qx, Qy to HKL

    sample.orientationMatrixINV = np.linalg.inv(np.dot(sample.RotMat3D,sample.orientationMatrix))
    

    if figure is None:
        fig = plt.figure(figsize=(7, 4))
    else:
        fig = figure
    def calculateTicks(ticks,angle,round=True):
        val = ticks/np.tan(angle/2.0)
        if round:
            return np.array(np.round(val),dtype=int)
        else:
            val


    if not nbinsx is None or not nbinsy is None: # Either nbinsx or nbinsy is provided (or both)
        if not nbinsy is None and nbinsx is None:
            nbinsx = int(np.round(nbinsy/np.sin(sample.projectionAngle)))
        elif not nbinsx is None and nbinsy is None:
            nbinsy = int(np.round(nbinsx*np.sin(sample.projectionAngle)))

        grid_locator1 = MaxNLocator(nbins=nbinsx)
        grid_locator2 = MaxNLocator(nbins=nbinsy)
        baseGrid = False
        locatorGrid = True
    
    elif not basex is None or not basey is None: # Either basex or basey is provided (or both)
        if basex is None:
            basex = calculateTicks(basey,sample.projectionAngle,round=False)
        elif basey is None:
            basey = basex/calculateTicks(1.0,sample.projectionAngle,round=False)
        baseGrid = True
        locatorGrid = False
        grid_locator1 = MultipleLocator(base=basex)
        grid_locator2 = MultipleLocator(base=basey)
    else:
        basex = 0.25
        basey = 0.25
        baseGrid = True
        locatorGrid = False
        grid_locator1 = MultipleLocator(base=basex)
        grid_locator2 = MultipleLocator(base=basey)
        
    grid_helper = GridHelperCurveLinear((sample.inv_tr, sample.tr),grid_locator1=grid_locator1,grid_locator2=grid_locator2)
    
    ax = SubplotHost(fig, *ids, grid_helper=grid_helper)
    ax.sample = sample
    
    ax.baseGrid = baseGrid
    ax.locatorGrid = locatorGrid
    if baseGrid:
        ax.basex = basex
        ax.basey = basey
    else:
        ax.yticks = nbinsy
        ax.xticks = nbinsx

    def set_axis(ax,v1,v2,*args):
        if not args is ():
            points = np.concatenate([[v1,v2],[x for x in args]],axis=0)
        else:
            points = np.array([v1,v2])
            
        if points.shape[1] == 3:
            points = ax.sample.calculateHKLtoProjection(points[:,0],points[:,1],points[:,2]).T
        boundaries = np.array([ax.sample.inv_tr(x[0],x[1]) for x in points])
        ax.set_xlim(boundaries[:,0].min(),boundaries[:,0].max())
        ax.set_ylim(boundaries[:,1].min(),boundaries[:,1].max())
        ax.forceGridUpdate()


    fig.add_subplot(ax)
    ax.set_aspect(1.)
    ax.grid(True, zorder=0)
    
    if not np.isclose(sample.projectionAngle,np.pi/2.0,atol=0.001):
        ax.axis["top"].major_ticklabels.set_visible(True)
        ax.axis["right"].major_ticklabels.set_visible(True)

    ax.format_coord = sample.format_coord
    ax.set_axis = lambda v1,v2,*args: set_axis(ax,v1,v2,*args)
    ax.set_xlabel('{} [RLU]'.format(', '.join([str(x) for x in sample.projectionVector1.astype(int)])))
    ax.set_ylabel('{} [RLU]'.format(', '.join([str(x) for x in sample.projectionVector2.astype(int)])))

    def forceGridUpdate(self):
        self._grid_helper._force_update = True
        self.pchanged()
        self.stale = True

    ax.calculateTicks = lambda value:calculateTicks(value,sample.projectionAngle)
    ax.forceGridUpdate = lambda:forceGridUpdate(ax)
    ax._oldXlimDiff = np.diff(ax.get_xlim())
    ax._oldYlimDiff = np.diff(ax.get_ylim())

    def get_aspect(ax):
        figW, figH = ax.get_figure().get_size_inches()
        _, _, w, h = ax.get_position().bounds
        disp_ratio = (figH * h) / (figW * w)
        data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())
        return (disp_ratio / data_ratio)

    ax.get_aspect_ratio = lambda: get_aspect(ax)

    def axisChanged(axis,forceUpdate=False,direction='both'):
        """Function to recalcualte the base number for RLU axis"""
        s = axis.sample
        xlim = axis.get_xlim()
        ylim = axis.get_ylim()
        xlimDiff = np.diff(xlim)
        ylimDiff = np.diff(ylim)
        if axis.baseGrid and not np.any([np.isclose(axis._oldXlimDiff,xlimDiff),
            np.isclose(axis._oldYlimDiff,ylimDiff)]) or forceUpdate:# Check if gridding with base and if new image size

            Q1 = np.array([s.tr(xlim[0],ylim[0])])
            Q2 = np.array([s.tr(xlim[1],ylim[1])])
            xspan,yspan = np.abs(Q1.T-Q2.T)

            if direction.lower() == 'y' or direction.lower()=='both':
                if not hasattr(axis,'yticks'):
                    yticks = 7
                else:
                    yticks = axis.yticks
                axis._oldYlimDiff = ylimDiff
                l2 = axis._grid_helper.grid_finder.grid_locator2

                if not hasattr(axis,'ybase'):
                    ybase = calculateBase(l2,yspan,yticks)
                else:
                    ybase = ax.ybase
                l2._base._base = ybase

            if direction.lower() == 'x' or direction.lower()=='both':
                if not hasattr(axis,'xticks'):
                    try:
                        yticks
                    except:
                        yticks = 7
                    xticks = int(0.60*axis.calculateTicks(yticks)/(axis.get_aspect_ratio()))
                else:
                    xticks = axis.xticks
                axis._oldXlimDiff = xlimDiff
                
                l1 = axis._grid_helper.grid_finder.grid_locator1
                if not hasattr(axis,'xbase'):
                    xbase = calculateBase(l1,xspan,xticks)
                else:
                    xbase = ax.xbase
                l1._base._base = xbase

            
        elif axis.locatorGrid:
            return
        else:
            return 
        
        # force an update
        axis.forceGridUpdate()
    ax.callbacks.connect('xlim_changed', axisChanged)
    ax.callbacks.connect('ylim_changed', axisChanged)
    ax.callbacks.connect('draw_event',lambda ax: axisChanged(ax,forceUpdate=True))
    ax.axisChanged = lambda direction='both': axisChanged(ax,forceUpdate=True,direction=direction)
    
    @updateAxisDecorator(ax=ax,direction='x')
    def set_xticks_number(xTicks,ax=ax):
        """Setter of the approximate number of x ticks to be used for plotting

        Args:

            - xTicks (int): Approximate number of ticks

        """
        ax.xticks = xTicks
        if hasattr(ax,'xbase'):
            del ax.xbase
    ax.set_xticks_number = set_xticks_number


    @updateAxisDecorator(ax=ax,direction='x')
    def set_xticks_base(xBase,ax=ax):
        """Setter of the base x ticks to be used for plotting

        Args:

            - xBase (float): Base of the tick marks

        """
        ax.xbase = xBase
        if hasattr(ax,'xticks'):
            del ax.xticks
    ax.set_xticks_base = set_xticks_base
    
    @updateAxisDecorator(ax=ax,direction='y')
    def set_yticks_number(yTicks,ax=ax):
        """Setter of the approximate number of x ticks to be used for plotting

        Args:

            - xTicks (int): Approximate number of ticks

        """
        ax.yticks = yTicks
        if hasattr(ax,'ybase'):
            del ax.ybase
    ax.set_yticks_number = set_yticks_number

    @updateAxisDecorator(ax=ax,direction='y')
    def set_yticks_base(yBase,ax=ax):
        """Setter of the base y ticks to be used for plotting

        Args:

            - yBase (float): Base of the tick marks

        """
        ax.ybase = yBase
        if hasattr(ax,'yticks'):
            del ax.yticks
    ax.set_yticks_base = set_yticks_base


    return ax




def createQEAxes(DataSet=None,axis=0,figure = None, projectionVector1 = None, projectionVector2 = None):
    """Function to create Q E plot

    Kwargs:

        - DataSet (DataSet): If provided and no projections vectors creates QE axis for main direction (default None)

        - axis (int): Whether to create axis 0 or 1 (projection vector 0 or orthogonal to this, default 0)

        - figure (figure): If provided, this is used to create the axis withing (default None)

        - projectionVector1 (vec): Projection vector along wich data is plotted. If not provided sample vector is used (default None)

        - projectionVector2 (vec): Projection vector orthogonal to data. If not provided sample vector is used (default None)


    """
    
    if projectionVector1 is None or projectionVector2 is None:
        v1 = DataSet.sample.projectionVector1
        v2 = DataSet.sample.projectionVector2
        angle = DataSet.sample.projectionAngle
        orientationMatrix = DataSet.sample.orientationMatrix
    else:
        v1 = np.asarray(projectionVector1)
        v2 = np.asarray(projectionVector2)
        if not np.all([x[0].shape==3 for x in [v1,v2]]) and not np.all([len(x.shape)==1 for x in [v1,v2]]):
            raise AttributeError('Provided vector(s) is not 3D: projectionVector1.shape={} or projectionVector2.shape={}'.format(v1.shape,v2.shape))
        angle = np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
        orientationMatrix = np.ones(3)

    v2Length = np.linalg.norm(v2)
    
    
    projectionMatrix = np.linalg.inv(np.array([[1,0],[np.cos(angle)*v2Length,np.sin(angle)*v2Length]]).T)
    
    projectionVectorQX = np.dot(np.dot(projectionMatrix,[1,0]),np.array([v1,v2]))
    projectionVectorQY = np.dot(np.dot(projectionMatrix,[0,1]),np.array([v1,v2]))
    projectionVectorQX = _tools.LengthOrder(projectionVectorQX)
    projectionVectorQY = _tools.LengthOrder(projectionVectorQY)
    projectionVectorQXLength = np.linalg.norm(np.dot(orientationMatrix,projectionVectorQY))
    projectionVectorQYLength = np.linalg.norm(np.dot(orientationMatrix,projectionVectorQX))
    projectionVectorQXFormated = ', '.join(['{:.3f}'.format(x) for x in projectionVectorQX])
    projectionVectorQYFormated = ', '.join(['{:.3f}'.format(x) for x in projectionVectorQY])
    
    if axis == 0:
        projectionVectorLength = projectionVectorQYLength
        projectionVectorLengthORthogonal = projectionVectorQXLength
        projectionVectorFormated = projectionVectorQXFormated
        projectionVector = projectionVectorQX
        projectionVectorOrthogonal = projectionVectorQY
    elif axis == 1:
        projectionVectorLength = projectionVectorQXLength
        projectionVectorFormated = projectionVectorQYFormated
        projectionVectorLengthORthogonal = projectionVectorQYLength
        projectionVector = projectionVectorQY
        projectionVectorOrthogonal = projectionVectorQX
    else:
        raise AttributeError('Provided axis of {} is not allowed. Should be eiter 0 or 1.'.format(axis))

    if figure is None:
        
        figure = plt.figure(figsize=(7, 4))
    else:
        figure.clf()
    def inv_tr(l,x,y):
        return x*l,y
    
    def tr(l,x,y):
        return x/l,y
    
    
    
    grid_helper = GridHelperCurveLinear((lambda x,y:inv_tr(projectionVectorLength,x,y), 
                                        lambda x,y:tr(projectionVectorLength,x,y)))
    
    ax = SubplotHost(figure, 1, 1, 1, grid_helper=grid_helper)
    
    figure.add_subplot(ax)
    #ax.set_aspect(1.)
    ax.grid(True, zorder=0)
    
    def calculateRLU(l,v1,x,y,v,step):
        return np.asarray(x)/l*v1+v*step, np.asarray(y)

    def format_coord(x,y): # pragma: no cover # x is H,K,L and y is  energy
        xformated = ', '.join(['{} = {}'.format(Y[0],Y[1]) for Y in zip(['h','k','l'],['{:.4f}'.format(X) for X in x])])
        return '{}, E, {:.4f}'.format(xformated,y)
    
    
    ax.set_xlabel('{} [RLU]'.format(projectionVectorFormated))
    
    ax.set_ylabel('E [meV]')
    ax._length = projectionVectorLengthORthogonal
    ax._projectionVector = projectionVector 
    ax._projectionVectorOrthogonal = projectionVectorOrthogonal
    ax._step = 0.0
    ax.calculateRLU = lambda x,y: calculateRLU(projectionVectorLength,ax._projectionVector,x,y,ax._projectionVectorOrthogonal,ax._step)
    ax.format_coord = lambda x,y: format_coord(*ax.calculateRLU(x,y))
    return ax







def calculateBase(l,span,ticks):
    """Calcualte the tick mark base suitable for current span and ticks

    Args:

        - l (grid locator): Matplotlib grid locator

        - span (float): Width of view

        - ticsk (int): Number of ticks wanted

    Returns:

        - base (float): Closest base number accorting to l.multiples
    """
    ytickorder = np.ceil(np.log10(span/ticks))
    minimalMultiplesy = np.argmin(np.abs(np.power(10,-ytickorder)*span/ticks-l.multiples))
    return l.multiples[minimalMultiplesy]*np.power(10,ytickorder)[0]
    


def updateAxisDecorator(ax,direction='both'):
    def axisDecorator(func):
        @functools.wraps(func)
        def newFunc(*args,**kwargs):
            returnval = func(*args,**kwargs)
            ax.axisChanged(direction=direction)
            return returnval
        return newFunc
    return axisDecorator


class MaxNLocator(mticker.MaxNLocator):
    def __init__(self, nbins=10, steps=None,
                 trim=True,
                 integer=False,
                 symmetric=False,
                 prune=None):
        # trim argument has no effect. It has been left for API compatibility
        mticker.MaxNLocator.__init__(self, nbins, steps=steps,
                                     integer=integer,
                                     symmetric=symmetric, prune=prune)
        self.create_dummy_axis()
        self._factor = None

    def __call__(self, v1, v2):
        if self._factor is not None:
            self.set_bounds(v1*self._factor, v2*self._factor)
            locs = mticker.MaxNLocator.__call__(self)
            return np.array(locs), len(locs), self._factor
        else:
            self.set_bounds(v1, v2)
            locs = mticker.MaxNLocator.__call__(self)
            return np.array(locs), len(locs), None

    def set_factor(self, f):
        self._factor = f


class MultipleLocator(mticker.MultipleLocator):
    def __init__(self,base=0.25):
        
        mticker.MultipleLocator.__init__(self, base)
        self.create_dummy_axis()
        self._factor = None
        self._multiplerVals = np.array([1,2,4,5,10])
        self.multiples = 1.0/self.multiplerVals

    @property
    def multiplerVals(self):
        return self._multiplerVals

    @multiplerVals.getter
    def multiplerVals(self):
        return self._multiplerVals

    @multiplerVals.setter
    def multiplerVals(self,multiplerVals):
        self._multiplerVals = multiplerVals
        self.multiples = 1.0/multiplerVals


    def __call__(self, v1, v2):
        if self._factor is not None:
            self.set_bounds(v1*self._factor, v2*self._factor)
            locs = mticker.MultipleLocator.__call__(self)
            return np.array(locs), len(locs), self._factor
        else:
            self.set_bounds(v1, v2)
            locs = mticker.MultipleLocator.__call__(self)
            return np.array(locs), len(locs), None

    def set_factor(self, f):
        self._factor = f
