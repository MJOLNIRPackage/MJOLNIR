import copy
import functools
import os
import sys
from operator import sub

pythonVersion = sys.version_info[0]

sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')

import matplotlib.pyplot as plt
#if pythonVersion == 3: # Only for python 3
import matplotlib.ticker as mticker

import numpy as np
from MJOLNIR import _tools
from matplotlib.patches import Ellipse
from mpl_toolkits.axisartist import SubplotHost#,axislines
import warnings
try:
    from mpl_toolkits.axisartist.grid_helper_curvelinear import \
        GridHelperCurveLinear
except ImportError:
    from mpl_toolkits.axisartist.floating_axes.grid_helper_curvelinear import \
        GridHelperCurveLinear



def calculateBase(l,span,ticks):
    """Calcualte the tick mark base suitable for current span and ticks

    Args:

        - l (grid locator): Matplotlib grid locator

        - span (float): Width of view

        - ticks (int): Number of ticks wanted

    Returns:

        - base (float): Closest base number according to l.multiples
    """
    ytickorder = np.ceil(np.log10(span/ticks))
    minimalMultiplesy = np.argmin(np.abs(np.power(10,-ytickorder)*span/ticks-l.multiples))
    return l.multiples[minimalMultiplesy]*np.power(10,ytickorder)
    


def updateAxisDecorator(ax,direction='both'):
    def axisDecorator(func):
        @functools.wraps(func)
        def newFunc(*args,**kwargs):
            returnval = func(*args,**kwargs)
            ax.axisChanged(direction=direction)
            return returnval
        return newFunc
    return axisDecorator

def updateXAxisDecorator(ax):
    def axisDecorator(func):
        @functools.wraps(func)
        def newFunc(*args,**kwargs):
            returnval = func(*args,**kwargs)
            ax.xAxisChanged()
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
        super(MaxNLocator,self).__init__(nbins, steps=steps,
                                    integer=integer,
                                    symmetric=symmetric, prune=prune)
        self.create_dummy_axis()
        self._factor = 1#0.0#None

    def __call__(self, v1, v2): # pragma: no cover
        if self._factor is not None:
            self.set_bounds(v1*self._factor, v2*self._factor)
            locs = mticker.MaxNLocator.__call__(self)
            return np.array(locs), len(locs), self._factor
        else:
            self.set_bounds(v1, v2)
            locs = mticker.MaxNLocator.__call__(self)
            return np.array(locs), len(locs)

    def set_factor(self, f):
        self._factor = f


class MultipleLocator(mticker.MultipleLocator):
    def __init__(self,base=None):
        if base is None:
            base = 0.25
        super(MultipleLocator,self).__init__(base)
        self.create_dummy_axis()
        self._factor = 1#0.0
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


    def __call__(self, v1, v2): # pragma: no cover
        if self._factor is not None:
            self.axis.set_view_interval(v1*self._factor, v2*self._factor)
            locs = mticker.MultipleLocator.__call__(self)
            return np.array(locs), len(locs), self._factor
        else:
            self.axis.set_view_interval(v1, v2)
            locs = mticker.MultipleLocator.__call__(self)
            return np.array(locs), len(locs), None

    def set_factor(self, f):
        self._factor = f

def forceGridUpdate(self):
    self._grid_helper._force_update = True
    self.pchanged()
    self.stale = True

def get_aspect(ax):
    figW, figH = ax.get_figure().get_size_inches()
    _, _, w, h = ax.get_position().bounds
    disp_ratio = (figH * h) / (figW * w)
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())
    return (disp_ratio / data_ratio)

def axisChanged(event,axis=None,forceUpdate=False,direction='both'):
    """Function to recalculate the base number for RLU axis"""
    # Two arguments are needed; event for canvas callback and axis is needed for ax.axisChanged
    if axis is None:
        axis = event
    s = axis.sample
    xlim = axis.get_xlim()
    ylim = axis.get_ylim()
    xlimDiff = np.diff(xlim)
    ylimDiff = np.diff(ylim)
    
    if direction == 'both':
        direction = ['x','y']
        difference = [xlimDiff,ylimDiff]
        locators = [axis._grid_helper.grid_finder.grid_locator1,axis._grid_helper.grid_finder.grid_locator2]
    else:
        if direction == 'x':
            difference = [xlimDiff]  
            locators = [axis._grid_helper.grid_finder.grid_locator1]
        else:
            difference = [ylimDiff]  
            locators = [axis._grid_helper.grid_finder.grid_locator2]
        direction = [direction]
    for d,diff,locator in zip(direction,difference,locators):
        forceUpdate = True
        if isinstance(locator,MultipleLocator):
            if not np.isclose(getattr(axis,'_old{}limDiff'.format(d.capitalize())),diff) \
            or forceUpdate:# if new image size

                points = np.array(np.meshgrid(xlim,ylim)).reshape(-1,2) # Generate the 4 corner points
                Qs = np.array([s.tr(p[0],p[1]) for p in points])
                if d=='x':
                    span = np.max(Qs[:,0])-np.min(Qs[:,0])
                else:
                    span = np.max(Qs[:,1])-np.min(Qs[:,1])
                if not hasattr(axis,'{}ticks'.format(d)):
                    setattr(axis,'{}ticks'.format(d),7)
                
                ticks = getattr(axis,'{}ticks'.format(d))
                setattr(axis,'_old{}limDiff'.format(d.capitalize()),diff)
            

            if not hasattr(axis,'{}base'.format(d)):
                base = calculateBase(locator,span,ticks)
            else:
                base = getattr(axis,'{}base'.format(d))

            locator.set_params(base=base)
        
        elif isinstance(locator,MaxNLocator):
            if hasattr(axis,'{}ticks'.format(d)):
                ticks = getattr(axis,'{}ticks'.format(d))
            else:
                ticks = 7
            locator.set_params(nbins = ticks)
        else:
            return 
    
    # force an update
    axis.forceGridUpdate()



def createQAxis(self,rlu = True, withoutOnClick = False, figure=None,ids=[1, 1, 1],basex=None,basey=None):#,ax=None):
    """Create a reciprocal lattice plot for a given DataSet object.
    
    Args:
        
        - Dataset (DataSet): DataSet object for which the RLU plot is to be made.

    Kwargs:

        - rlu (bool): If true, plot in reciprocal lattice units, else in Qx, Qy (default True)

        - withoutOnClick (bool): If False, utilize the updated onClick method (default False)

        - figure: Matplotlib figure in which the axis is to be put (default None)

        - ids (array): List of integer numbers provided to the SubplotHost ids attribute (default [1,1,1])

        - basex (float): Ticks are positioned at multiples of this value along x (default None)

        - basey (float): Ticks are positioned at multiples of this value along y (default None)

    Returns:
        
        - ax (Matplotlib axes): Axes containing the RLU plot.

    .. note::
        When rlu axis is created, the orientation of Qx and Qy is assumed to be rotated as well. 
        This is to be done in the self.View3D method call!

    .. note::
        When using python 2 the changing of tick marks is not supported due to limitations in matplotlib. However, if python 3 is used, the number 
        of ticks and their location can be change after the initialization using the set_xticks_number, set_yticks_number chaning the wanted number 
        of tick marks, or the set_xticks_base or set_yticks_base to change the base number, see RLU tutorial under Tools. As default a sufficient base
        number is found and will update when zooming.
        
    """
    
    sample = copy.deepcopy(self.sample)
    if rlu:
        for samp in sample:
            samp.convert = np.einsum('ij,j...->i...',samp.RotMat,samp.convert)
            #sample.convert = np.einsum('ij,j...->i...',sample.RotMat,sample.convert)
            samp.convertinv = np.linalg.inv(samp.convert) # Convert from Qx, Qy to projX, projY

            samp.orientationMatrix = np.dot(samp.RotMat3D,samp.orientationMatrix)
            samp.orientationMatrixINV = np.linalg.inv(samp.orientationMatrix)
            samp.theta = 0.0

    if figure is None:# and ax is None:
        fig = plt.figure(figsize=(7, 4))
    else:
        fig = figure#ax.get_figure()
    if rlu:
        def calculateTicks(ticks,angle,round=True):
            val = ticks/np.tan(angle/2.0)
            if round:
                return np.array(np.round(val),dtype=int)
            else:
                return val

        
        if  not basex is None or not basey is None: # Either basex or basey is provided (or both)
            if basex is None:
                basex = calculateTicks(basey,sample[0].projectionAngle,round=False)
            elif basey is None:
                basey = basex/calculateTicks(1.0,sample[0].projectionAngle,round=False)

            grid_locator1 = MultipleLocator(base=basex)
            grid_locator2 = MultipleLocator(base=basey)
        else:
            basex = 0.5
            basey = 0.5

            grid_locator1 = MultipleLocator(base=basex)
            grid_locator2 = MultipleLocator(base=basey)
            
        grid_helper = GridHelperCurveLinear((sample[0].inv_tr, sample[0].tr),grid_locator1=grid_locator1,grid_locator2=grid_locator2)
        

        #if not ax is None:
        #    if not isinstance(ax,axislines.Axes):
        #        subargs=ax.get_subplotspec()
        #        fig.delaxes(ax)
        #        ax = SubplotHost(fig, subargs, grid_helper=grid_helper)
        #        
        #else:
        #    #ax2=ax1
        ax = SubplotHost(fig, *ids, grid_helper=grid_helper)
        #fig.add_subplot(ax)
    else:
        ax = fig.add_subplot(111)
    ax.sample = sample[0]
    ax.rlu = rlu
    
    
    if ax.rlu:
        
        ax.basex = basex
        ax.basey = basey

        def set_axis(ax,v1,v2,*args):
            if args != ():
                points = np.concatenate([[v1,v2],[x for x in args]],axis=0)
            else:
                points = np.array([v1,v2])
                
            if points.shape[1] == 3:
                points = ax.sample.calculateHKLtoProjection(points[:,0],points[:,1],points[:,2]).T
            boundaries = np.array([ax.sample.inv_tr(x[0],x[1]) for x in points])
            ax.set_xlim(boundaries[:,0].min(),boundaries[:,0].max())
            ax.set_ylim(boundaries[:,1].min(),boundaries[:,1].max())
            if pythonVersion == 3: # Only possible in python 3
                ax.forceGridUpdate()


        fig.add_subplot(ax)
        ax.set_aspect(1.)
    ax.grid(True, zorder=0)
    
    if ax.rlu:
        if not np.isclose(ax.sample.projectionAngle,np.pi/2.0,atol=0.001):
            ax.axis["top"].major_ticklabels.set_visible(True)
            ax.axis["right"].major_ticklabels.set_visible(True)

    if rlu:
        ax.format_coord = ax.sample.format_coord
    else:
        def format_coord(x,y):
            return ', '.join([label+' = {:.3f}'.format(val) for val,label in zip([x,y],['qx','qy'])])

        ax.format_coord = format_coord
    
    
    if ax.rlu:
        ax.set_axis = lambda v1,v2,*args: set_axis(ax,v1,v2,*args)

        def beautifyLabel(vec):
            Vec = [x.astype(int) if np.isclose(x.astype(float)-x.astype(int),0.0) else x.astype(float) for x in vec]
            return '{} [RLU]'.format(', '.join([str(x) for x in Vec]))

        ax.set_xlabel(beautifyLabel(ax.sample.projectionVector1))
        ax.set_ylabel(beautifyLabel(ax.sample.projectionVector2))

    else:
        ax.set_xlabel('Qx [1/AA]')
        ax.set_ylabel('Qy [1/AA]')
    
    def onclick(ax, event,Qx,Qy,data,outputFunction): # pragma: no cover
        if event.xdata is not None and ax.in_axes(event):
            try:
                C = ax.get_figure().canvas.cursor().shape() # Only works for pyQt5 backend
            except:
                pass
            else:
                if C != 0:# or ax.drawState != States.INACTIVE:
                    return
        if not ax.in_axes(event):
            return
        printString = ''
        printString+=ax.format_coord(event.xdata, event.ydata)+', '
        QX = np.array(Qx[0])
        QY = np.array(Qy[0])
        
        Qx = 0.25*(QX[1:,1:]+QX[:-1,:-1]+QX[1:,:-1]+QX[:-1,1:])
        Qy = 0.25*(QY[1:,1:]+QY[:-1,:-1]+QY[1:,:-1]+QY[:-1,1:])

        arg = np.argmin(np.linalg.norm(np.array([Qx,Qy])-np.array([event.xdata,event.ydata]).reshape(2,1,1),axis=0))
        arg2D = np.unravel_index(arg,Qx.shape)
        
        cts = data[0][0][arg2D[0],arg2D[1]]
        Norm = data[1][0][arg2D[0],arg2D[1]]
        Mon = data[2][0][arg2D[0],arg2D[1]]
        NC = data[3][0][arg2D[0],arg2D[1]]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            Intensity = np.divide(cts*NC,Norm*Mon)

        localSize = np.linalg.norm(np.array([QX[arg2D[0],arg2D[1]]-Qx[arg2D[0],arg2D[1]],QY[arg2D[0],arg2D[1]]-Qy[arg2D[0],arg2D[1]]]))

        if not np.isfinite(Intensity) or np.linalg.norm(np.array([Qx[arg2D[0],arg2D[1]]-event.xdata,Qy[arg2D[0],arg2D[1]]-event.ydata]))>localSize:
            printString+='I = NaN'
        else:
            printString+='I = {:.4E}'.format(Intensity)
            printString+=', Cts = {:d}, Norm = {:.3f}, Mon = {:d}, NormCount = {:d}'.format(cts,Norm,int(Mon),NC)
        if not ax.suppressPrint:
            outputFunction(printString)

    if not withoutOnClick:
        ax.onClick = lambda x: onclick(ax,x,ax.Qx,ax.Qy,[ax.intensity,ax.monitorCount,ax.Normalization,ax.NormCount],outputFunction=ax.outputFunction)
        ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event', ax.onClick)
    


    if ax.rlu: # 
        ax.calculateTicks = lambda value:calculateTicks(value,ax.sample.projectionAngle)
        ax.forceGridUpdate = lambda:forceGridUpdate(ax)
        ax._oldXlimDiff = np.diff(ax.get_xlim())
        ax._oldYlimDiff = np.diff(ax.get_ylim())

        ax.get_aspect_ratio = lambda: get_aspect(ax)

        ax.callbacks.connect('xlim_changed', axisChanged)
        ax.callbacks.connect('ylim_changed', axisChanged)
        ax.get_figure().canvas.mpl_connect('draw_event',lambda event: axisChanged(event,axis=ax,forceUpdate=True))
        ax.axisChanged = lambda direction='both': axisChanged(None,ax,forceUpdate=True,direction=direction)
    
        @updateAxisDecorator(ax=ax,direction='x')
        def set_xticks_base(xBase,ax=ax):
            """Setter of the base x ticks to be used for plotting

            Args:

                - xBase (float): Base of the tick marks

            """
            if not isinstance(ax._grid_helper.grid_finder.grid_locator1,MultipleLocator):
                l1 = MultipleLocator(base=xBase)
                ax._grid_helper.update_grid_finder(grid_locator1=l1)

            ax.xbase = xBase

        @updateAxisDecorator(ax=ax,direction='y')
        def set_yticks_base(yBase,ax=ax):
            """Setter of the base y ticks to be used for plotting

            Args:

                - yBase (float): Base of the tick marks

            """
            if not isinstance(ax._grid_helper.grid_finder.grid_locator2,MultipleLocator):
                l2 = MultipleLocator(base=yBase)
                ax._grid_helper.update_grid_finder(grid_locator2=l2)
            ax.ybase = yBase

        @updateAxisDecorator(ax=ax,direction='x')
        def set_xticks_number(xNumber,ax=ax):
            """Setter of the number of x ticks to be used for plotting

            Args:

                - xNumber (int): Number of x tick marks

            """
            if not isinstance(ax._grid_helper.grid_finder.grid_locator1,MaxNLocator):
                l1 = MaxNLocator(nbins=xNumber)
                ax._grid_helper.update_grid_finder(grid_locator1=l1)
            ax.xticks = xNumber

        @updateAxisDecorator(ax=ax,direction='y')
        def set_yticks_number(yNumber,ax=ax):
            """Setter of the number of y ticks to be used for plotting

            Args:

                - yNumber (int): Number of y tick marks

            """
            if not isinstance(ax._grid_helper.grid_finder.grid_locator2,MaxNLocator):
                l2 = MaxNLocator(nbins=yNumber)
                ax._grid_helper.update_grid_finder(grid_locator2=l2)
            ax.yticks = yNumber

        ax.set_xticks_base = set_xticks_base
        ax.set_yticks_base = set_yticks_base
        ax.set_xticks_number = set_xticks_number
        ax.set_yticks_number = set_yticks_number

    return ax




def createQEAxes(DataSet=None,axis=0,figure = None, projectionVector1 = None, projectionVector2 = None, rlu = True, figurePos=[1,1,1]):
    """Function to create Q E plot

    Kwargs:

        - DataSet (DataSet): If provided and no projections vectors creates QE axis for main direction (default None)

        - axis (int): Whether to create axis 0 or 1 (projection vector 0 or orthogonal to this, default 0)

        - figure (figure): If provided, this is used to create the axis within (default None)

        - projectionVector1 (vec): Projection vector along which data is plotted. If not provided sample vector is used (default None)

        - projectionVector2 (vec): Projection vector orthogonal to data. If not provided sample vector is used (default None)
        
        - rlu (Bool): Boolean to control if rlu or qxqy (default True)


    """
        
    if projectionVector1 is None or projectionVector2 is None:
        if rlu is True:
            v1 = DataSet.sample[0].projectionVector1
            v2 = DataSet.sample[0].projectionVector2
            angle = DataSet.sample[0].projectionAngle
            orientationMatrix = DataSet.sample[0].orientationMatrix
        else:
            v1 = np.array([1.0,0.0])
            v2 = np.array([0.0,1.0])
            angle = np.pi/2.0
            orientationMatrix = np.ones(2)
    else:
        v1 = np.array(projectionVector1)
        v2 = np.array(projectionVector2)

        if not np.all([x.shape==(3,) for x in [v1,v2]])\
        or not np.all([len(x.shape)==1 for x in [v1,v2]])\
        or np.all([x.shape==(2,) for x in [v1,v2]]):
            raise AttributeError('Provided vector(s) is not 3D/2D: projectionVector1.shape={} or projectionVector2.shape={}'.format(v1.shape,v2.shape))
        angle = np.arccos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
        if rlu is True:
            orientationMatrix = np.ones(3)
        else:
            orientationMatrix = np.ones(2)


    sample = copy.deepcopy(DataSet.sample)


    v2Length = np.linalg.norm(v2)/np.linalg.norm(v1)
    projectionMatrix = np.linalg.inv(np.array([[1,0],[np.cos(angle)*v2Length,np.sin(angle)*v2Length]]).T)

    projectionVectorQX = np.dot(np.dot(projectionMatrix,[1,0]),np.array([v1,v2]))
    projectionVectorQY = np.dot(np.dot(projectionMatrix,[0,v2Length]),np.array([v1,v2]))

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
        raise AttributeError('Provided axis of {} is not allowed. Should be either 0 or 1.'.format(axis))

    if figure is None:

        figure = plt.figure(figsize=(7, 4))
    #else:
    #    figure.clf()
    def inv_tr(l,x,y):
        return x*l,y

    def tr(l,x,y):
        return x/l,y


    grid_locator1 = MultipleLocator(base=1.0) # Standard X ticks is multiple locator
    grid_helper = GridHelperCurveLinear((lambda x,y:inv_tr(projectionVectorLength,x,y), 
                                    lambda x,y:tr(projectionVectorLength,x,y)),grid_locator1=grid_locator1)

    if rlu:
        labels = ['h','k','l']
    else:
        labels = ['qx','qy']

    ax = SubplotHost(figure, *figurePos, grid_helper=grid_helper)
    
    ax.inv_tr = lambda x,y:inv_tr(projectionVectorLength,x,y)
    ax.tr = lambda x,y:tr(projectionVectorLength,x,y)
    
    ax.sample = sample[0]

    figure.add_subplot(ax)
    #ax.set_aspect(1.)
    ax.grid(True, zorder=0)

    def calculateRLU(l,v1,x,y,v,step):
        return np.asarray(x)/l*v1+v*step, np.asarray(y)

    def format_coord(x,y,rlu=rlu,labels=labels): # pragma: no cover # x is H,K,L and y is  energy
        xformated = ', '.join(['{} = {}'.format(Y[0],Y[1]) for Y in zip(labels,['{:.4f}'.format(X) for X in x])])
        return '{}, E={:.4f}'.format(xformated,y)


    if rlu:
        ax.set_xlabel('{} [RLU]'.format(projectionVectorFormated))
    else:
        ax.set_xlabel('{} [1/A]'.format(projectionVectorFormated))

    ax.set_ylabel('E [meV]')
    ax._length = projectionVectorLengthORthogonal
    ax._projectionVector = projectionVector 
    ax._projectionVectorOrthogonal = projectionVectorOrthogonal
    ax.plotDirection = projectionVector
    ax._step = 0.0
    ax.calculateRLU = lambda x,y: calculateRLU(projectionVectorLength,ax._projectionVector,x,y,ax._projectionVectorOrthogonal,ax._step)
    ax.format_coord = lambda x,y: format_coord(*ax.calculateRLU(x,y))


    ax.forceGridUpdate = lambda:forceGridUpdate(ax)
    ax.xticks = 7

    def xAxisChanged(event, axis=None, forceUpdate=False):
        # Two arguments are needed; event for canvas callback and axis is needed for ax.xAxisChanged
        if axis is None:
            axis = event
        locator = axis._grid_helper.grid_finder.grid_locator1
        xlim = axis.get_xlim()
        xlimDiff = np.diff(xlim)
        if isinstance(locator,MultipleLocator):
            if hasattr(axis,'xBase'):
                base = axis.xBase
            else:
                base = calculateBase(locator,xlimDiff,axis.xticks)
            locator.set_params(base=base)

        elif isinstance(locator,MaxNLocator):
            if hasattr(axis,'xTicks'):
                ticks = getattr(axis,'xTicks')
            else:
                ticks = 7
            locator.set_params(nbins = ticks)
        else:
            return
        axis.forceGridUpdate()

    ax.callbacks.connect('xlim_changed', xAxisChanged)

    ax.get_figure().canvas.mpl_connect('draw_event',lambda event: xAxisChanged(event,axis=ax,forceUpdate=True))
    ax.xAxisChanged = lambda: xAxisChanged(None,ax,forceUpdate=True)


    @updateXAxisDecorator(ax=ax)
    def set_xticks_base(xBase=None,ax=ax):
        """Setter of the base x ticks to be used for plotting

        Kwargs:

            - xBase (float): Base of the tick marks (default automatic)

        """


        if not isinstance(ax._grid_helper.grid_finder.grid_locator1,MultipleLocator):
            l1 = MultipleLocator(base=xBase)
            ax._grid_helper.update_grid_finder(grid_locator1=l1)

        if xBase is None:
            if hasattr(ax,'xBase'):
                delattr(ax,'xBase')
        else:
            ax.xBase = xBase

    @updateXAxisDecorator(ax=ax)
    def set_xticks_number(xNumber = None,ax=ax):
        """Setter of the number of x ticks to be used for plotting

        Kwargs:

            - xNumber (int): Number of x tick marks (default 7)

        """
        if xNumber is None:
            xNumber = 7

        if not isinstance(ax._grid_helper.grid_finder.grid_locator1,MaxNLocator):
            l1 = MaxNLocator(nbins=xNumber)
            ax._grid_helper.update_grid_finder(grid_locator1=l1)
        ax.xTicks = xNumber


    ax.set_xticks_base = set_xticks_base
    ax.set_xticks_number = set_xticks_number

    return ax






