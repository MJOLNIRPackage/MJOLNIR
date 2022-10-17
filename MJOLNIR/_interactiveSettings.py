"""Settings for the interactive plotting used in the MJOLNIR package. Collected here for a common reference and easier adaptability and future-proofing."""

import matplotlib as mpl
import numpy as np
from enum import Enum
from matplotlib.patches import Ellipse,Rectangle,Circle
from matplotlib.backend_bases import MouseButton
#try:
#    import matplotlib.backends.backend_qt5.cursors as cursors
#    mplVersion = None
#except ImportError:
from  matplotlib.backend_tools import Cursors as cursors
#    mplVersion = 3.5
mplVersion = float('.'.join(mpl.__version__.split('.')[:2]))
import MJOLNIR
import os
import PyQt5.QtCore
from PyQt5.QtGui import QCursor, QPixmap

import sys
from collections import defaultdict


## Viewer 3D

Viewer3DSettings = {'upwards':         ['+','up'  ,'right'],
                    'upwardsScroll':   ['up'],
                    'downwards':       ['-','down','left'],
                    'downwardsScroll': ['down'],
                    'home':            ['home'],
                    'end':             ['end'],
                    'QxE':             ['0'],
                    'QyE':             ['1'],
                    'QxQy':            ['2'],
#                    'resolution':      ['r']
                    }

# 1D rectangular cuts

interactive1DKeys = {
                 'cutting':  ['c'],
                 'resolution': ['r']
}

resolutionSettings = {
                        'return': ['b'],
                        'clear':['c'],
                     }

cut1DSettings = {
                 'new':  ['n'],
                 'manage': ['m'],
                 'cut':  ['c'],
                 'return': ['b'],
                }


def setupModes(ax):

    # Set up mode to be a private variable with value INACTIVE
    ax._mode = None#_interactiveSettings.Modes.INACTIVE
    
    def setMode(ax,modeName):
        if hasattr(ax,'isActive'):
            if not ax.isActive:
                return
        if not ax._mode is None: # A previous mode was set, deactivate it
            if ax._mode.upper() == modeName.upper(): # if new mode is the same as old, return
                return
            getattr(sys.modules[__name__],'deactivate'+ax._mode.upper())(ax)    
            ax.get_figure().canvas.draw_idle()
        
        if not modeName.upper() in [x.name for x in Modes]:
            raise AttributeError('Wanted state "{:}" not recognized. Must be one of: '.format(modeName)+','.join([x.name for x in Modes]))
        getattr(sys.modules[__name__],'initialize'+modeName.upper())(ax)
        ax._mode = modeName.upper()
        setCursor(ax,pointerType[modeName.upper()])
        ax.get_figure().canvas.draw_idle()
        

    ax.setMode = lambda stateName: setMode(ax,stateName)
    ax.setMode('inactive')
    
    return ax




def setCursor(ax,cursor):
    if mplVersion >= 3.5:
        pass#ax.get_figure().canvas.set_cursor(cursor)
    else:
        if isinstance(cursor,type(cursors.HAND)): # Standard Matplotlib cursor
            ax.get_figure().canvas.toolbar.set_cursor(cursor)
        else:
            ax.get_figure().canvas.setCursor(cursor)
    

interactive1DKeysReversed = {}
for key,value in interactive1DKeys.items():
    for v in value:
        interactive1DKeysReversed[v] = key


interactive1DKeysAll = list(np.concatenate(list(interactive1DKeys.values())))
# Concatenate all possible settings for 1D rectangular cuts
cut1DSettingsAll = list(np.concatenate(list(cut1DSettings.values())))

# kwargs for plotting of the rectangles
cut1DKkwargs = {'linewidth':1, 'edgecolor':'r', 'facecolor':'none','zorder':22}

# Global variable to give back generated 1DCuts to the user in scripting mode
global cut1DHolder 
cut1DHolder = []   

# Colors for selected and deselected cuts
selectedColor = (1.0,0.0,0.0,1)
deselectedColor = (0.0,1.0,0.0,1)




# Interactive Modes

modes = {}
# Possible modes are :
#      - inactive
#      - manager: Move or delete 1D cuts if any are present
#      - cutting: Create new 1D cut
#      - resolution: Overplot resolution ellipses

for i,t in enumerate(['inactive','cutting','resolution']):
    modes[t.upper().replace(' ','')]=i
Modes= Enum('Modes', modes)

#CutImageLocation = os.path.join(os.path.split(MJOLNIR.Data.__file__)[0],'scissors.png')
#CutCursor = QCursor(QPixmap(CutImageLocation))


## Cursor type mode
pointerType = defaultdict(lambda: cursors.POINTER)
# pointerType['RESOLUTION'] = resolutionCursor

pointerType['CUTTING_INACTIVE'] = cursors.SELECT_REGION#CutCursor
pointerType['CUTTING_EMPTY'] = cursors.SELECT_REGION#PyQt5.QtCore.Qt.ForbiddenCursor

if mplVersion >= 3.5:
    pointerType['CUTTING_INITIAL'] = cursors.SELECT_REGION
    pointerType['CUTTING_WIDTH'] = cursors.SELECT_REGION
    pointerType['CUTTING_DIRECTION'] = cursors.SELECT_REGION
    pointerType['CUTTING_MOVE'] = cursors.HAND
else:
    pointerType['CUTTING_INITIAL'] = PyQt5.QtCore.Qt.CrossCursor
    pointerType['CUTTING_WIDTH'] = PyQt5.QtCore.Qt.CrossCursor
    pointerType['CUTTING_DIRECTION'] = PyQt5.QtCore.Qt.CrossCursor
    pointerType['CUTTING_MOVE'] = PyQt5.QtCore.Qt.SizeAllCursor



def resetUsingKey(ax,keys):
    """make onKeyPress function which resets using any of the provided keys, always include escape"""
    if not 'escape' in keys:
        keys.append('escape')
    if hasattr(ax,'_key_press_event'): # if a previous onKeyPres has been defind remove it.
        ax.figure.canvas.mpl_disconnect(ax._key_press_event)
    def onKeyPress(event,ax):
        if not event.key in keys:
            return
        ax.setMode('inactive')
    ax._key_press_event = ax.figure.canvas.mpl_connect('key_press_event',lambda event: onKeyPress(event, ax) )


def initializeINACTIVE(ax):
    if hasattr(ax,'_key_press_event'):
        ax.figure.canvas.mpl_disconnect(ax._key_press_event)

    if not hasattr(ax,'suppressPrint'):
        ax.suppressPrint = False

    def onKeyPress(event,ax):
        possibleKeys = interactive1DKeysReversed.keys()
        if not event.key in possibleKeys:
            return
        wantedState = interactive1DKeysReversed[event.key]
        ax.setMode(wantedState)
    ax._key_press_event = ax.figure.canvas.mpl_connect('key_press_event',lambda event: onKeyPress(event, ax) )


    return

def deactivateINACTIVE(ax):
    return


def initializeRESOLUTION(ax):
    

    # Change cursor to signify new mode
    
    # Reset using same key as activated it
    resetUsingKey(ax,resolutionSettings['return'])
    # Setup colours needed for plotting
    if not hasattr(ax,'EiColors'):
        colors = np.array(['#ff7f0e',
                            '#2ca02c',
                            '#d62728',
                            '#9467bd',
                            '#8c564b',
                            '#e377c2',
                            '#7f7f7f',
                            '#bcbd22',
                            '#17becf',
                            '#1f77b4'])
        # Any energies close than 0.01 are regarded equal
        uniqueEi = np.unique(np.round([df.Ei for df in ax.ds],2))
        ax.EiColors = {}
        for i,ei in enumerate(uniqueEi):
            ax.EiColors[ei] = colors[np.mod(i,len(colors))]


    # List of ellipses to plot
    ax.resolutionEllipses = []


    ax.old_format_coord = ax.format_coord
    def format_coord(ax,x,y):
        if hasattr(ax,'isActive'):
            if not ax.isActive:
                return
        returnString = ax.old_format_coord(x,y)
        if ax.type in ['QE','QELine','QELineView3D']:
            Eis = ax.ds.findEi(y)
        elif ax.type == 'QPlane':
            deltaE = np.mean([ax.EMin,ax.EMax])
            Eis = ax.ds.findEi(deltaE)
            

        if ax.type in ['QE','QELine','QELineView3D','QPlane']:

            if ax.get_autoscale_on():
                ax.autoscale(False)
            if ax.type in ['QE','QELine','QELineView3D']:
                if ax.type == 'QELineView3D': # Currently only done with rlu
                    dirVector = ax.sample.calculateHKLToQxQy(*ax.calculateRLU(1,0)[0])-ax.sample.calculateHKLToQxQy(*ax.calculateRLU(0,0)[0])#ax.calculateRLU(0,0)[0]

                elif ax.type == 'QE':
                    if ax.rlu:
                        dirVector = ax.ds[0].sample.calculateHKLToQxQy(*ax.endPoint)-ax.ds[0].sample.calculateHKLToQxQy(*ax.startPoint)
                    else:
                        dirVector = ax.endPoint-ax.startPoint
                else:
                    dirVector = np.diff(ax.QPoints,axis=0)[ax.calculateIndex(x)[0]]

                if ax.type == 'QELineView3D':
                    pos = ax.calculateRLU(x,0)[0]
                else:
                    pos = ax.calculatePosition(x).flatten()
            else: # QPlane
                if ax.rlu:
                    pos = ax.sample.calculateQxQyToHKL(x,y).flatten()
                else:
                    pos = np.array([x, y])
            
            returnString+= ' [resolution]'
            # If the number of ellipses in ax.resolutionEllipses does not match needed for this energy, fix it!
            if len(Eis)!=len(ax.resolutionEllipses):

                for E in ax.resolutionEllipses:
                    E.set_visible(False)
                    E.remove()

                ax.resolutionEllipses = []
                for _ in range(len(Eis)):
                        ax.resolutionEllipses.append(Ellipse(xy=np.array([0,0]),width=0.0,height=0.0,edgecolor=[0.0,0.0,0.0,0.0],
                        facecolor=[0.0,0.0,0.0,0.0],zorder=22))
                        ax.add_patch(ax.resolutionEllipses[-1])
            #elif len(Eis)<len(ax.resolutionEllipses): # too many ellipses, delete some
            #    removeLength = len(ax.resolutionEllipses)-len(Eis)
            #    for E in ax.resolutionEllipses[-removeLength:]:
            #        patch = ax.patches[ax.patches.index(E)]
            #        E.remove()
            #        patch.remove()
            #        del ax.resolutionEllipses[-1]
                
            if ax.type in ['QELine','QELineView3D','QE']:
                if ax.type == 'QE':
                    P1 = np.array([*dirVector/np.linalg.norm(dirVector)**2,0.0,0.0])
                else:
                    P1 = np.array([*dirVector/np.linalg.norm(dirVector),0.0,0.0])
                # Energy direction
                P2 = np.array([0.0,0.0,0.0,1.0])
               

                P = np.zeros((3,5),dtype='float')
                P[0,:-1] = P1
                P[1,:-1] = P2
                P[-1,-1] = 1
                for Ei,E in zip(Eis,ax.resolutionEllipses):
                    Ef = Ei-y
                    M = ax.ds.calculateResolutionMatrix(pos,Ei,Ef,rlu=ax.rlu)
                    Q = np.diag([0,0,0,0.0,-1.0])
                    Q[:4,:4] = M
                    
                    # Calculate projection onto orthogonal directions
                    C = np.linalg.inv(np.dot(P,np.dot(np.linalg.inv(Q),P.T)))
                    #C[0]/=np.linalg.norm(dirVector)**2

                    eigenValues,eigenVectors = np.linalg.eig(C[:2,:2])
                    sigma = np.power(eigenValues,-0.5)
                    # M,eigenVectors,sigma = ax.ds.calculateResolutionMatrixAndVectors(pos,P1,P2,Ei,Ef,rlu=ax.rlu,rluAxis=True)
                    
                    ellipseColor = ax.EiColors[Ei]
                    angle = -np.rad2deg(np.arctan2(*eigenVectors[0,::-1])) 
                    E.set_center(np.array([x,y]))
                    E.set_width(2.0*sigma[0])
                    E.set_height(2.0*sigma[1])
                    E.set_angle(angle)
                    E.set_edgecolor(ellipseColor)
                    E.set_visible(True)
                ax.get_figure().canvas.draw_idle()
            
            elif ax.type == 'QPlane':
                
                if ax.rlu:
                    P1 = np.array([*ax.sample.calculateHKLToQxQy(*ax.sample.projectionVector1),0.0,0.0])
                    P1*=1.0/np.linalg.norm(P1)
                    
                    P2 = np.array([*ax.sample.calculateHKLToQxQy(*ax.sample.projectionVector2),0.0,0.0])
                    P2*=1.0/np.linalg.norm(P2)
                else:
                    P1 = np.array([1.0,0.0,0.0,0.0])
                    P2 = np.array([0.0,1.0,0.0,0.0])

                for Ei,E in zip(Eis,ax.resolutionEllipses):
                    Ef = Ei-deltaE
                    _,eigenVectors,sigma = ax.ds.calculateResolutionMatrixAndVectors(pos,P1,P2,Ei,Ef,rlu=ax.rlu,rluAxis=True)
                    ellipseColor = ax.EiColors[Ei]
                    angle = -np.rad2deg(np.arctan2(*eigenVectors[0,::-1])) 
                    E.set_center(np.array([x,y]))
                    E.set_width(2.0*sigma[0])
                    E.set_height(2.0*sigma[1])
                    E.set_angle(angle)
                    E.set_edgecolor(ellipseColor)
                    E.set_visible(True)
                #    E.draw(renderer)
                ax.get_figure().canvas.draw_idle()
                
        return returnString
    
    ax.format_coord = lambda x,y: format_coord(ax,x,y)

    #if not hasattr(ax,'_button_press_event'): # View3D
    #    ax.
    #else:
    ax.old_onClick = ax.onClick
    if not hasattr(ax,'_ellipseHolder'):
        ax._ellipseHolder = []

    if hasattr(ax,'_button_press_event'):
        ax.figure.canvas.mpl_disconnect(ax._button_press_event)
        
    if ax.rlu:
        ax.labels = ['H','K','L','E']
    else:
        ax.labels = ['Qx','Qy','E']

    ### Generate the onClick function
    def onClick(ax,event):
        if hasattr(ax,'isActive'):
            if not ax.isActive:
                return
        

        x,y = event.xdata,event.ydata
        printString = ax.old_format_coord(x,y)

        if ax.type in ['QE','QELine','QELineView3D']:

            Eis = ax.ds.findEi(y)
            
            if ax.type == 'QELineView3D': # Currently only done with rlu
                #if not hasattr(ax,'dirVector'):
                dirVector = ax.sample.calculateHKLToQxQy(*ax.calculateRLU(1,0)[0])-ax.sample.calculateHKLToQxQy(*ax.calculateRLU(0,0)[0])#ax.calculateRLU(0,0)[0]

            elif ax.type == 'QE':
                if ax.rlu:
                    dirVector = ax.ds[0].sample.calculateHKLToQxQy(*ax.endPoint)-ax.ds[0].sample.calculateHKLToQxQy(*ax.startPoint)
                else:
                    dirVector = ax.endPoint-ax.startPoint
            else:
                dirVector = np.diff(ax.QPoints,axis=0)[ax.calculateIndex(x)[0]]

            if ax.type == 'QELineView3D':
                pos = ax.calculateRLU(x,0)[0]
            else:
                pos = ax.calculatePosition(x).flatten()
            P = np.zeros((3,5))
            P[-1,-1] = 1
            ## Along cut
            if ax.type == 'QE': # In the QE axis, x is rescaled with np.linalg.norm(dirVector). Following takes that into account
                P[0,:2] = dirVector/np.linalg.norm(dirVector)**2
            else:
                P[0,:2] = dirVector/np.linalg.norm(dirVector)
            # Energy direction
            P[1,:-1] = np.array([0.0,0.0,0.0,1.0])
            
            
            for Ei in Eis:
                Ef = Ei-y
                M = ax.ds.calculateResolutionMatrix(pos,Ei,Ef,rlu=ax.rlu)

                Q = np.diag([0,0,0,0.0,-1.0])
                Q[:4,:4] = M
                
                # Calculate projection onto orthogonal directions
                C = np.linalg.inv(np.dot(P,np.dot(np.linalg.inv(Q),P.T)))
                
                eigenValues,eigenVectors = np.linalg.eig(C[:2,:2])
            
                sigma = np.power(eigenValues,-0.5)
                
                theta = np.linspace(0,np.pi*2,201)
                
                ellipseColor = ax.EiColors[Ei]


                # Calculate projection onto orthogonal directions
                C = np.linalg.inv(np.dot(P,np.dot(np.linalg.inv(Q),P.T)))
                
                
                #if ax.type == 'QE':
                #    C[0]/=np.linalg.norm(dirVector)
                eigenValues,eigenVectors = np.linalg.eig(C[:2,:2])
                sigma = np.power(eigenValues,-0.5)
                    
                ellipsePoints = eigenVectors[:,0].reshape(2,1)*sigma[0]*np.cos(theta).reshape(1,-1)+eigenVectors[:,1].reshape(2,1)*sigma[1]*np.sin(theta).reshape(1,-1)+np.array([x,y]).reshape(2,1)
                

                ax._ellipseHolder.append(ax.plot(*ellipsePoints,zorder=21,color=ellipseColor,label=', '.join([l+' = {:.2f}'.format(x) for l,x in zip(ax.labels,[*pos,y])])+'(Ei={:.2f})'.format(Ei))[0])
                ax.get_figure().canvas.draw()


                if ax.type == 'QE':  # Recalculate the projection due to plotting method of the QE axis
                    P[0,:2] = dirVector/np.linalg.norm(dirVector)
                    C = np.linalg.inv(np.dot(P,np.dot(np.linalg.inv(Q),P.T)))
                    eigenValues,eigenVectors = np.linalg.eig(C[:2,:2])
                    sigma = np.power(eigenValues,-0.5)
                printString+='\nReso (Ei={:.2f}):\nM = '.format(Ei)+str(M)+'\nIn projection (In units of 1/AA and meV):\nV1 with sigma={:.4f}, V1 = '.format(sigma[0])+str(eigenVectors[:,0].flatten())+\
                '\nV2 with sigma={:.4f}, V2 = '.format(sigma[1])+str(eigenVectors[:,1].flatten())

        elif ax.type == 'QPlane':
            if ax.rlu:
                pos = ax.sample.calculateQxQyToHKL(x,y).flatten()
            else:
                pos = np.array([x, y])
            
            if ax.rlu:
                P1 = np.array([*ax.sample.calculateHKLToQxQy(*ax.sample.projectionVector1),0.0,0.0])
                P1*=1.0/np.linalg.norm(P1)
                
                P2 = np.array([*ax.sample.calculateHKLToQxQy(*ax.sample.projectionVector2),0.0,0.0])
                P2*=1.0/np.linalg.norm(P2)
            else:
                P1 = np.array([1.0,0.0,0.0,0.0])
                P2 = np.array([0.0,1.0,0.0,0.0])
            
            deltaE = np.mean([ax.EMin,ax.EMax])
            Eis = ax.ds.findEi(deltaE)
            for Ei in Eis:
                Ef = Ei-deltaE
                M,eigenVectors,sigma = ax.ds.calculateResolutionMatrixAndVectors(pos,P1,P2,Ei,Ef,rlu=ax.rlu,rluAxis=True)

                theta = np.linspace(0,np.pi*2,201)
                
                ellipseColor = ax.EiColors[Ei]
                
                ellipsePoints = eigenVectors[:,0].reshape(2,1)*sigma[0]*np.cos(theta).reshape(1,-1)+eigenVectors[:,1].reshape(2,1)*sigma[1]*np.sin(theta).reshape(1,-1)+np.array([event.xdata,event.ydata]).reshape(2,1)
                ax._ellipseHolder.append(ax.plot(*ellipsePoints,zorder=21,color=ellipseColor,label=', '.join([l+' = {:.2f}'.format(x) for l,x in zip(ax.labels,[*pos,deltaE])])+'(Ei={:.2f})'.format(Ei))[0])
                ax.get_figure().canvas.draw()
                printString+='\nReso (Ei={:.2f}, Ef={:.2f}, deltaE={:.2f}):\nM = '.format(Ei,Ef,deltaE)+str(M)+'\nIn projection (In units of 1/AA and meV):\nV1 with sigma={:.4f}, V1 = '.format(sigma[0])+str(eigenVectors[:,0].flatten())+\
                '\nV2 with sigma={:.4f}, V2 = '.format(sigma[1])+str(eigenVectors[:,1].flatten())
        if not ax.suppressPrint:
            ax.outputFunction(printString)

    
    def on_key_press(self,event):# pragma: no cover
        if hasattr(self,'isActive'):
            if not self.isActive:
                return
        if event.key in resolutionSettings['clear']: # clear the ellipses
            if hasattr(ax,'_ellipseHolder'):
                for p in ax._ellipseHolder:
                    p.set_visible(False)
                    p.remove()
                else:
                    ax._ellipseHolder = []

    ax.ellipse_reset_cid = ax.figure.canvas.mpl_connect(
        'key_press_event', lambda event:on_key_press(ax,event))

    ax.onClick = lambda event: onClick(ax,event)
    ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event',ax.onClick)



def deactivateRESOLUTION(ax):
    
    ax.format_coord = ax.old_format_coord
    for E in list(ax.resolutionEllipses)+[x for x in ax.patches if isinstance(x,Ellipse)]:
        E.set_visible(False)
        try:
            E.remove()
        except (AttributeError,ValueError):
            pass

        try:
            ax.resolutionEllipses[ax.resolutionEllipses.index(E)].remove()
        except:
            pass
    

    ax.figure.canvas.mpl_disconnect(ax._button_press_event)
    ax.get_figure().canvas.draw_idle()
    

    if hasattr(ax,'ellipse_reset_cid'):
        ax.figure.canvas.mpl_disconnect(ax.ellipse_reset_cid)
    ax.onClick = ax.old_onClick
    ax._button_press_event = ax.figure.canvas.mpl_connect('button_press_event',ax.onClick)
    


    ax.get_figure().canvas.draw_idle()



def CuttingModeCursorManagerHover(axes,event):
    if hasattr(axes,'isActive'):
        if not axes.isActive:
            return
    if not hasattr(axes,'drawState'):
        setCursor(axes,pointerType['CUTTING_EMPTY'])
        return
    state = 'CUTTING_'+axes.drawState.name
    
    if state == 'CUTTING_INACTIVE' and axes.new>0:
        setCursor(axes,pointerType['CUTTING_WIDTH'])
    else:
        setCursor(axes,pointerType[state])
    
    

def initializeCUTTING(ax):
        
    resetUsingKey(ax,cut1DSettings['return'])
    

    if ax.type in ['QE','QELineView3D']:
        Draggables = [DraggableRectanglePerpendicular,DraggableRectangleHorizontal,DraggableRectangleVertical]

        
        def cut1DFunctionRectangleDefault(self,dr):
            global cut1DHolder
            parameters = extractCut1DPropertiesRectanglePerpendicular(dr.rect,self.ds.sample[0])
            
            # Convert center point into actual position in Q
            if ax.type == 'QELineView3D':
                middlePoint = ax.calculateRLU(*parameters['center'])[0]
                if ax.rlu:
                    dirVector = ax.calculateRLU(1,0)[0]-ax.calculateRLU(0,0)[0]
                    orthogonalVector = np.cross(self.ds.sample[0].planeNormal,dirVector)
                    orthogonalVector*=1/np.linalg.norm(orthogonalVector)
                else:
                    # TODO: fix
                    orthogonalVector = np.array([ax.plotDirection[1],-ax.plotDirection[0]])
            else:
                middlePoint = ax.calculatePosition(parameters['center'][0]).T
                if ax.rlu:
                    orthogonalVector = np.cross(self.ds.sample[0].planeNormal,ax.plotDirection.flatten())
                    orthogonalVector*=1/np.linalg.norm(orthogonalVector)
                else:
                    orthogonalVector = np.array([ax.plotDirection[1],-ax.plotDirection[0]])
            
            del parameters['center'] # remove the 'center' as it is not allowed in plotCut1D

            # transform the orthogonal vector if needed 
            
            
            parameters['q1'] = middlePoint
            parameters['q2'] = middlePoint+orthogonalVector

            parameters['minPixel'] = ax.minPixel
            cut1DHolder.append([self.ds.plotCut1D(**parameters,extend=True)])
        
        def cut1DFunctionRectangleHorizontalDefault(self,dr):
            global cut1DHolder
            parameters = extractCut1DPropertiesRectangleHorizontal(dr.rect,self.ds.sample[0])

            if ax.type == 'QELineView3D':
                parameters['q1'] = ax.calculateRLU(*parameters['q1'])[0]
                parameters['q2'] = ax.calculateRLU(*parameters['q2'])[0]
            else:
                # Convert center point into actual position in Q
                parameters['q1'] = ax.calculatePosition(parameters['q1'][0]).T
                parameters['q2'] = ax.calculatePosition(parameters['q2'][0]).T
                
            parameters['minPixel'] = ax.minPixel
            parameters['width'] = ax.width


            cut1DHolder.append([self.ds.plotCut1D(**parameters)])

        def cut1DFunctionRectangleVerticalDefault(self,dr):
            global cut1DHolder
            parameters = extractCut1DPropertiesRectangleVertical(dr.rect,self.ds.sample[0])
            
            
            # Convert center point into actual position in Q
            if ax.type == 'QELineView3D':
                parameters['q'] = ax.calculateRLU(parameters['q'],0.0)[0]
            else:
                # Convert center point into actual position in Q
                parameters['q'] = ax.calculatePosition(parameters['q']).T
            
            
            parameters['minPixel'] = ax.dE
            parameters['width'] = ax.width
            cut1DHolder.append([self.ds.plotCut1DE(**parameters)])
        

            

        if ax.cut1DFunctionRectanglePerpendicular is None:

            ax.cut1DFunctionRectanglePerpendicular = lambda dr: cut1DFunctionRectangleDefault(ax,dr)
        

        if ax.cut1DFunctionRectangleHorizontal is None:

            ax.cut1DFunctionRectangleHorizontal = lambda dr: cut1DFunctionRectangleHorizontalDefault(ax,dr)
        

        if ax.cut1DFunctionRectangleVertical is None:
            ax.cut1DFunctionRectangleVertical = lambda dr: cut1DFunctionRectangleVerticalDefault(ax,dr)


        DraggableFunctions = [ax.cut1DFunctionRectanglePerpendicular,ax.cut1DFunctionRectangleHorizontal,ax.cut1DFunctionRectangleVertical]

        prepareInteractiveCutting(ax,Draggables,DraggableFunctions)

    if ax.type == 'QPlane':
        Draggables = [DraggableCircle,DraggableRectangle]

        ## Define the cut functions to be called
        def cut1DFunctionRectangleDefault(self,dr):
            global cut1DHolder
            parameters = extractCut1DPropertiesRectangle(dr.rect,self.sample)
            
            EMin = self.EMin
            EMax = self.EMax
            cut1DHolder.append([self.ds.plotCut1D(**parameters,Emin=EMin,Emax=EMax)])


        def cut1DFunctionCircleDefault(self,dr):
            global cut1DHolder
            parameters = extractCut1DPropertiesCircle(dr.circ,self.sample)
            parameters['E1'] = self.ds.energy.min()
            parameters['E2'] = self.ds.energy.max()
            parameters['minPixel'] = self.EMax-self.EMin
            
            cut1DHolder.append([self.ds.plotCut1DE(**parameters)])

        if ax.cut1DFunctionRectangle is None:
            ax.cut1DFunctionRectangle = lambda dr: cut1DFunctionRectangleDefault(ax,dr)
        
        
        if ax.cut1DFunctionCircle is None:
            ax.cut1DFunctionCircle = lambda dr: cut1DFunctionCircleDefault(ax,dr)
        
        

        DraggableFunctions = [ax.cut1DFunctionCircle,ax.cut1DFunctionRectangle]

        prepareInteractiveCutting(ax,Draggables,DraggableFunctions)


    if not hasattr(ax,'_old_format_coord'):
        ax._old_format_coord = ax.format_coord
    else:
        raise AttributeError('While initializeCUTTING, ax already has _old_format_coord')
    
    if not hasattr(ax,'draggableShapes'):
        def format_coord(ax,x,y):
            return '[CUTTING not possible] '+ ax._old_format_coord(x,y)
        axes = ax
        setCursor(ax,pointerType['CUTTING_EMPTY'])
    else:
        def format_coord(ax,x,y):
            return '[CUTTING] '+ ax._old_format_coord(x,y)

        if hasattr(ax,'ax'):
            axes = ax.ax
        else:
            axes = ax
    


    ax.hoverID = ax.get_figure().canvas.mpl_connect('motion_notify_event', lambda event: CuttingModeCursorManagerHover(axes,event))

    ax.format_coord = lambda x,y: format_coord(ax,x,y)


    
    

    # if not hasattr(ax,'annotateBox'):
    #     ax.annotateBox = ax.annotate("Cutting", xy=(0,0), xytext=(0,3),textcoords="offset points",
    #                     bbox=dict(boxstyle='round4', fc='linen',ec='k',lw=1),
    #                     #arrowprops=dict(arrowstyle='-|>')
    #                     )
    
    # ax.annotateBox.set_visible(True)

    # def onMove(ax,event):
    #     if not ax.in_axes(event):
    #         ax.annotateBox.set_visible(False)
    #         return
    #     ax.annotateBox.set_visible(True)
        
    #     x = event.xdata
    #     y = event.ydata
        
    #     ax.annotateBox.xy = (x,y)
    
    #     ax.get_figure().canvas.draw_idle()
    
    # ax.onMoveAnnotate = lambda event: onMove(ax,event)
    # ax.onMoveAnnotateID = ax.get_figure().canvas.mpl_connect('motion_notify_event', ax.onMoveAnnotate)


def deactivateCUTTING(ax):
    # ax.get_figure().canvas.mpl_disconnect(ax.onMoveAnnotateID)
    # ax.annotateBox.set_visible(False)
    
    ax.get_figure().canvas.mpl_disconnect(ax.hoverID)
    
    if hasattr(ax,'draggableShapes'): # it was possible to perform cuts
        # if in a cutting state, revert it!
        if hasattr(ax,'ax'):
            axes = ax.ax
        else:
            axes = ax
        cancel(ax,axes)

        fig = ax.get_figure()
        del ax.resetFormatCoord
        del ax.new
        #??ax.cidmove = None
        del ax.drawState#= States.INACTIVE
        
        #ax.shapes = []
        #ax.drs = []
        #ax.draggableShapes = None#draggables
        #ax.draggableFunctions = None# draggableFunctions
        #del ax.ax
        
        #del ax.selectedDr# = None
        #del ax.lock# = True
        
        ax.suppressPrint = False
        
        fig.canvas.mpl_disconnect(fig._key_press_event)
        fig.canvas.mpl_disconnect(fig._button_press_event)

    
    ax.format_coord = ax._old_format_coord
    del ax._old_format_coord


## Cutting states
states = {}
# States are:
#    - inactive: Not possible to draw shape
#    - initial: No point has been selected
#    - direction: Start point selected, define direction
#    - width: Start and end points selected, define width 
for i,t in enumerate(['inactive','initial','direction','width','move']):
    states[t.upper().replace(' ','')]=i
States= Enum('States', states)



def cancel(self,axes):# pragma: no cover
    self.resetFormatCoord()
    self.newShape = None
    self.new = False
    
    # Only when something has been drawn is cid set
    if not self.cidmove is None:
        axes.get_figure().canvas.mpl_disconnect(self.cidmove)
        self.cidmove = None
        
        self.patches[-1].remove()
        self.get_figure().canvas.draw_idle()

    # Reset state and return
    self.drawState = States.INACTIVE

    return 



def on_key_press(self,event):# pragma: no cover
    if hasattr(self,'isActive'):
        if not self.isActive:
            return
    
    if not hasattr(self,'drawState'):
        return
    if hasattr(self,'ax'):
        axes = self.ax
    else:
        axes = self
    if not event.key in cut1DSettingsAll:
        self.resetFormatCoord()
        self.suppressPrint = False

        return
    if event.key in cut1DSettings['manage']:
        if hasattr(self,'ax'):
            axes = self.ax
        else:
            axes = self
            
        cancel(self,axes)
            
        if self.lock is None:
            self.lock = True
            self.resetFormatCoord()
        else:
            self.lock = None
            self.format_coord_old = self.format_coord
            self.format_coord = lambda x,y:self.format_coord_old(x,y)+' [manage]'
            self.suppressPrint = True
            self.drawState = States.MOVE

        return
    elif event.key in cut1DSettings['cut'] and self.new is False:
        self.lock=True
        self.resetFormatCoord()
        self.drawState = States.INACTIVE
        dr = self.selectedDr
        if not dr is None:
            dr.Cut1DFunction(dr=dr)

            return
    
    elif event.key in cut1DSettings['new']:
        self.lock=True
        if hasattr(self,'ax'):
            axes = self.ax
        else:
            axes = self
            
        if not self.drawState == States.INACTIVE:
            cancel(self,axes)
            return
            
        # If new is False, no new shapes will be initialized
        # If new is not False, the shape is the draggableShape with index self.new-1
        # Modulus ensures that there is an overflow back to shape 0 instead of out of bounce
        self.new = np.mod(self.new+1,len(self.draggableShapes)+1)#
        
        
        
        if self.new>0:
            self.selectedShape = self.draggableShapes[self.new-1]
            self.selectedShape.shapeSelected(self,axes)
            self.suppressPrint = True
            setCursor(axes,pointerType['CUTTINGCUT'])
        else:
            self.resetFormatCoord()
            setCursor(axes,pointerType['CUTTING'])
            
        self.lock = True
    


def on_press(self,event):# pragma: no cover
    if hasattr(self,'isActive'):
        if not self.isActive:
            return
    if not hasattr(self,'drawState'):
        return
    if hasattr(self,'ax'):
        axes = self.ax
    else:
        axes = self

    if not self.new > 0:
        return
    if event.inaxes != axes:
        return

    
    # Right-clicking cancels any ongoing action and reset to inactive
    if event.button == MouseButton.RIGHT:
        cancel(self,axes)
        return 


    if not self.drawState is None:
        if not self.cidmove is None: # If exists, disconnect it
            axes.get_figure().canvas.mpl_disconnect(self.cidmove)
            self.cidmove = None
        stateName = self.drawState.name.lower()
        # call the designated state function on the selected shape
        on_move = getattr(self.selectedShape,stateName)(self,axes,event)
        self.suppressPrint = True
        # connect the returned move function if not None. None is return in the last function 
        
        if not on_move is None:
            self.cidmove = axes.get_figure().canvas.mpl_connect('motion_notify_event',lambda event:on_move(self,event))
        else:
            self.resetFormatCoord()
        
        
def reset(self):# pragma: no cover
    if not hasattr(self,'drawState'):
        return
    self.suppressPrint=False
    if hasattr(self,'format_coord_old'):
        self.format_coord = self.format_coord_old
        del self.format_coord_old


def prepareInteractiveCutting(ax,draggables, draggableFunctions):# pragma: no cover
    fig = ax.get_figure()
    ax.resetFormatCoord = lambda: reset(ax)
    ax.new=False
    ax.cidmove = None
    ax.drawState = States.INACTIVE
    
    if not hasattr(ax,'shapes'):
        ax.shapes = []
        ax.drs = []

    ax.draggableShapes = draggables
    ax.draggableFunctions = draggableFunctions
    if not hasattr(ax,'ax'):
        ax.ax = ax
    
    if not hasattr(ax,'selectedDr'):
        ax.selectedDr = None
    ax.lock = True
    
    ax.suppressPrint = False
    
    fig._key_press_event = fig.canvas.mpl_connect(
        'key_press_event', lambda event:on_key_press(ax,event))
    
    fig._button_press_event = fig.canvas.mpl_connect(
                'button_press_event', lambda event:on_press(ax,event))


def prepareInteractiveCuttingView3D(self,draggables, draggableFunctions):# pragma: no cover
    fig = self.ax.get_figure()
    self.ax.resetFormatCoord = lambda: reset(self.ax)
    self.ax.new=False
    self.ax.cidmove = None
    self.ax.drawState = States.INACTIVE
    
    self.ax.shapes = []
    self.ax.drs = []
    self.ax.draggableShapes = draggables
    self.ax.draggableFunctions = draggableFunctions
    self.ax.ax = self.ax
    
    self.ax.selectedDr = None
    self.ax.lock = True
    
    self.ax.suppressPrint = False
    
    fig.canvas.mpl_connect(
        'key_press_event', lambda event:on_key_press(self.ax,event))
    
    fig.canvas.mpl_connect(
                'button_press_event', lambda event:on_press(self.ax,event))



#### All the shapes... ############################################################################

class DraggableShape(): # pragma: no cover
    lock = None  # only one can be moved at a time
    selectedShape = None # Only one can be selected at one point
    cidkeypress = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def disconnect(self):
        pass
        #'disconnect all the stored connection ids'
        #for idx in ['cidpress','cidrelease','cidmotion']:
        #    try:
        #        self.rect.figure.canvas.mpl_disconnect(getattr(self,idx))
        #    except AttributeError:
        #        pass

        #self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        #self.rect.figure.canvas.mpl_disconnect(self.cidmotion)


            
class DraggableRectangle(DraggableShape):# pragma: no cover
    
    def __init__(self, rect,line,plottingObject,Cut1DFunction,figure):
        self.rect = rect
        self.line = line
        self.press = None
        self.background = None
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        self.figure = figure
        
    def remove(self):
        self.disconnect()
        fig = self.rect.figure
        self.rect.set_visible(False)
        self.rect.remove()
        self.line.set_visible(False)
        self.line.remove()
        if self.selected:
            self.figure.selectedDr = None
        #idx = self.figure.shapes.index(self)
        #self.figure.shapes[idx].remove()
        #del self.figure.shapes[idx]
        if not hasattr(fig,'canvas'):
            fig.get_figure().canvas.draw()
        else:
            fig.canvas.draw()


    @property
    def selected(self):
        return self._selected
    
    @selected.getter
    def selected(self):
        return self._selected
    
    @selected.setter
    def selected(self,value):
        if self._selected == bool(value):
            return
        if value:
            if not self.figure.selectedDr is None:
                self.figure.selectedDr.selected = False
            self.figure.selectedDr = self
            color = selectedColor
        else:
            color = deselectedColor
            
        self.rect.set_edgecolor(color)
        self.line.set_color(color)
        self.rect.figure.canvas.draw()
        self._selected = bool(value)



    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes: return
        if self.figure.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        
        figure = self.plottingObject
        if event.button == MouseButton.RIGHT and figure.drawState == States.MOVE:
            self.remove()
            return 
        self.selected = True
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        self.figure.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        self.line.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.figure.lock is not self:
            return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.rect.set_x(x0+dx)
        self.rect.set_y(y0+dy)

        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        # restore the background region
        canvas.restore_region(self.background)

        
        center = corner2center(*self.rect.get_xy(), self.rect.get_width(), self.rect.angle)
        
        height = self.rect.get_height()
        angle = np.deg2rad(self.rect.angle+90)
        dx,dy = height*np.array([np.cos(angle),np.sin(angle)])
        self.line.set_data([center[0],center[0]+dx],[center[1],center[1]+dy])
        
        # redraw line and rect with the current rectangle
        axes.draw_artist(self.rect)
        axes.draw_artist(self.line)
        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if self.figure.lock is not self:
            return

        self.press = None
        self.figure.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.line.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()


    ## Class method
    def shapeSelected(figure,axes):
        '''Called when figure has chosen this shape but not yet drawn it'''
        if hasattr(axes,'format_coord_old'):
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) + ' [rect]'
            
        else:
            axes.format_coord_old = axes.format_coord
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) + ' [rect]'

    def inactive(figure,axes,event):
        
        width = 0.0
        height = 0.05
        center = corner2center(event.xdata,event.ydata, width, -180.0)
        figure.newShape = Rectangle(center,width,height,**cut1DKkwargs) # 
        
        axes.add_patch(figure.newShape)
        
        axes.get_figure().canvas.draw()
        def on_move(self,event):
            if event.inaxes:
                rectx,recty = self.newShape.xy
                center = corner2center(rectx,recty, width, self.newShape.angle)
                mousex = event.xdata
                mousey = event.ydata
                dx,dy = center[0]-mousex, center[1]-mousey
                angle = np.arctan2(dy,dx)+np.pi/2
                
                self.newShape.set_xy(center2corner(*center, width, np.rad2deg(angle)))
                self.newShape.angle = np.rad2deg(angle)
                
                self.newShape.set_height(np.linalg.norm([dx,dy]))
                
                
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw just the current rectangle
                axes.draw_artist(self.newShape)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        figure.newShape.set_animated(True)
        axes.get_figure().canvas.draw()
        figure.background = figure.figure.canvas.copy_from_bbox(axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(figure.newShape)


        # and blit just the redrawn area
        axes.get_figure().canvas.blit(axes.bbox)
        
        figure.drawState = States.INITIAL

        return on_move

    def initial(figure,axes,_):
        center = corner2center(*figure.newShape.get_xy(), figure.newShape.get_width(), figure.newShape.angle)

        height = figure.newShape.get_height()
        angle = np.deg2rad(figure.newShape.angle+90)
        dx,dy = height*np.array([np.cos(angle),np.sin(angle)])
        
        
        figure.line = axes.plot([center[0],center[0]+dx],[center[1],center[1]+dy],cut1DKkwargs['edgecolor'],zorder=cut1DKkwargs['zorder']+1)[0]
                
        def on_move(self,event):
            if event.inaxes:
                
                rectx,recty = self.newShape.xy
                angle = np.deg2rad(self.newShape.angle+90)
                dirvec = np.array([np.cos(angle),np.sin(angle)])
                
                center = corner2center(rectx,recty, self.newShape.get_width(), self.newShape.angle)
                
                mousePoint = np.array([event.xdata-center[0],event.ydata-center[1]])
                
                ortho = mousePoint-np.dot(dirvec,mousePoint)*dirvec
                sign = np.sign(np.dot(ortho,np.array([dirvec[1],-dirvec[0]])))
                
                width = sign*np.linalg.norm(ortho)*2.0
                
                
                self.newShape.set_width(width)
                self.newShape.set_xy(center2corner(*center, width, self.newShape.angle))
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw both line and rectangle
                axes.draw_artist(self.newShape)
                axes.draw_artist(self.line)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        axes.get_figure().canvas.blit(axes.bbox)
        figure.drawState = States.WIDTH

        return on_move

    def width(figure,axes,_):
        
        rect = Rectangle(figure.newShape.get_xy(),width=figure.newShape.get_width(),
                         height=figure.newShape.get_height(),angle=figure.newShape.angle,
                         **cut1DKkwargs)
        
        #Cut1DFunction = lambda dr: figure.cut1DFunction(figure,dr)
        Cut1DFunction = figure.draggableFunctions[figure.draggableShapes.index(DraggableRectangle)]
        
        axes.add_patch(rect)
        dr = DraggableRectangle(rect,figure.line,figure,Cut1DFunction,figure)
        dr.connect()
        dr.selected = True
        
        figure.new = False
        figure.newShape = None
        figure.line = None
        #figure.shapes.append(rect)
        figure.drawState = States.INACTIVE
        

        figure.shapes.append(dr)

        return None

def corner2center(x,y,width,angle):# pragma: no cover
    # Rotate from rectangle space to real space
    angleRad = np.deg2rad(angle+90)
    # Unit vector along height of rectangle
    dirvec = np.array([np.cos(angleRad),np.sin(angleRad)])
    # Orthogonal with half width
    ortho = np.array([dirvec[1],-dirvec[0]])*width*0.5
    # center point is orthogonal to height and half the width
    center = np.array([x,y])+ortho
    return center

def center2corner(x,y,width,angle):# pragma: no cover
    # Rotate from rectangle space to real space
    angleRad = np.deg2rad(angle+90)
    # Unit vector along height of rectangle
    dirvec = np.array([np.cos(angleRad),np.sin(angleRad)])
    # Orthogonal with half width
    ortho = np.array([dirvec[1],-dirvec[0]])*width*0.5
    # corner point is orthogonal to height and half the width but reversed
    corner = np.array([x,y])-ortho
    return corner

def extractCut1DPropertiesRectangle(rect,sample=None, rounding = 4):# pragma: no cover
    # Extract parameters from a rectangle
    xy = rect.get_xy()
    width = rect.get_width()
    height = rect.get_height()
    angle = rect.angle
    
    angleRad = np.deg2rad(angle+90)
    # Unit vector along height of rectangle
    dirvec = np.array([np.cos(angleRad),np.sin(angleRad)])
    ortho = np.array([dirvec[1],-dirvec[0]])
    botCenter = np.array(xy)+ortho*0.5*width
    topCenter = np.array(xy)+ortho*0.5*width+dirvec*height
    
    if not sample is None:
        botCenter = sample.calculateQxQyToHKL(*botCenter)
        topCenter = sample.calculateQxQyToHKL(*topCenter)
    
    minPixel = 0.05
    rlu = not sample is None # RLU only when a sample is provided
    constantBins = False
    ufit = False
    params = {'q1':np.array([np.round(x,rounding) for x in botCenter]),
              'q2':np.array([np.round(x,rounding) for x in topCenter]),
              'width':np.round(np.abs(width),rounding),
              'minPixel':np.round(minPixel,rounding),
              'rlu':rlu,
              'constantBins':constantBins,
              'ufit':ufit}
    return params


def extractCut1DPropertiesRectanglePerpendicular(rect,sample=None, rounding = 4):# pragma: no cover
    # Extract parameters from a rectangle
    xy = rect.get_xy()
    width = rect.get_width()
    height = rect.get_height()
    center = np.asarray(xy)+0.5*np.array([width,height])
    
    rlu = not sample is None # RLU only when a sample is provided
    constantBins = False
    ufit = False
    params = {'center':np.array([np.round(x,rounding) for x in center]),
              'Emin':np.round(xy[1],rounding),
              'Emax':np.round(xy[1]+height,rounding),
              'width':np.round(width,rounding),
              'rlu':rlu,
              'constantBins':constantBins,
              'ufit':ufit}
    return params

def extractCut1DPropertiesRectangleHorizontal(rect,sample=None, rounding = 4):# pragma: no cover
    # Extract parameters from a rectangle
    xy = rect.get_xy()
    width = rect.get_width()
    height = rect.get_height()
    q1 = np.asarray(xy)+np.array([0.0,0.5*height])
    q2 = np.asarray(xy)+np.array([width,height])
    
    
    rlu = not sample is None # RLU only when a sample is provided
    constantBins = False
    ufit = False
    params = {'q1':np.array([np.round(x,rounding) for x in q1]),
              'q2':np.array([np.round(x,rounding) for x in q2]),
              'Emin':np.round(xy[1],rounding),
              'Emax':np.round(xy[1]+height,rounding),
              'rlu':rlu,
              'constantBins':constantBins,
              'ufit':ufit}
    return params

def extractCut1DPropertiesRectangleVertical(rect,sample=None, rounding = 4):# pragma: no cover
    # Extract parameters from a rectangle
    xy = rect.get_xy()
    width = rect.get_width()
    height = rect.get_height()
    q1 = np.asarray(xy)+np.array([0.5*width,0.0])
    
    Energies = [xy[1],xy[1]+height]
    EMin,EMax = np.min(Energies),np.max(Energies)
    
    rlu = not sample is None # RLU only when a sample is provided
    constantBins = False
    ufit = False
    params = {'q':np.round(q1[0],rounding),
              
              'E1':np.round(EMin,rounding),
              'E2':np.round(EMax,rounding),
              'rlu':rlu,
              'constantBins':constantBins,
              'ufit':ufit}
    return params

def extractCut1DPropertiesCircle(circ,sample=None, rounding = 4):# pragma: no cover
    # Extract parameters from a rectangle
    center = circ.center
    width = circ.radius*2
    
    if not sample is None:
        center = sample.calculateQxQyToHKL(*center)
    
    rlu = not sample is None # RLU only when a sample is provided
    constantBins = False
    ufit = False
    params = {'q':np.array([np.round(x,rounding) for x in center]),
              'width':np.round(np.abs(width),rounding),
              'rlu':rlu,
              'constantBins':constantBins,
              'ufit':ufit}
    return params


## Utility functions to make plots interactive

def clearBoxes(self):# pragma: no cover
    """Clear all generated draggable rectangles"""
    if hasattr(self,'ax'):
        axes = self.ax
    else:
        axes = self

    for dr in self.drs:
        dr.remove()
        
    self.drs = []
    axes.get_figure().canvas.draw()


class DraggableCircle(DraggableShape):# pragma: no cover
    
    def __init__(self, circ,plottingObject,Cut1DFunction,figure):
        self.circ = circ
        self.press = None
        self.background = None
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        self.figure = figure
        
    def remove(self):
        fig = self.figure
        self.disconnect()
        self.circ.set_visible(False)
        self.circ.remove()
        #del self.circ
        if self.selected:
            self.figure.selectedDr = None
        #idx = self.figure.shapes.index(self)
        #self.figure.shapes[idx].remove()
        # #del self.figure.shapes[idx]
        if not hasattr(fig,'canvas'): # This is an axis an not a figure!
            fig.get_figure().canvas.draw()        
        else:
            fig.canvas.draw()        

    @property
    def selected(self):
        return self._selected
    
    @selected.getter
    def selected(self):
        return self._selected
    
    @selected.setter
    def selected(self,value):
        if self._selected == bool(value):
            return
        if value:
            if not self.figure.selectedDr is None:
                self.figure.selectedDr.selected = False
            self.figure.selectedDr = self
            color = selectedColor
        else:
            color = deselectedColor
            
        self.circ.set_edgecolor(color)

        self.circ.figure.canvas.draw()
        self._selected = bool(value)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.circ.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.circ.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.circ.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)


    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.circ.axes: return
        if self.figure.lock is not None: return
        contains, attrd = self.circ.contains(event)
        if not contains: return
        
        figure = self.plottingObject
        if event.button == MouseButton.RIGHT and figure.drawState == States.MOVE:
            self.remove()
            return 

        self.selected = True
        x0, y0 = self.circ.center
        self.press = x0, y0, event.xdata, event.ydata
        self.figure.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.circ.figure.canvas
        axes = self.circ.axes
        self.circ.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.circ.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.circ)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.figure.lock is not self:
            return
        if event.inaxes != self.circ.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.circ.center = [x0+dx,y0+dy]

        canvas = self.circ.figure.canvas
        axes = self.circ.axes
        # restore the background region
        canvas.restore_region(self.background)


        # redraw line and rect with the current rectangle
        axes.draw_artist(self.circ)
        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if self.figure.lock is not self:
            return

        self.press = None
        self.figure.lock = None

        # turn off the rect animation property and reset the background
        self.circ.set_animated(False)
        self.background = None

        # redraw the full figure
        self.circ.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        try:
            self.circ.figure.canvas.mpl_disconnect(self.cidpress)
            self.circ.figure.canvas.mpl_disconnect(self.cidrelease)
            self.circ.figure.canvas.mpl_disconnect(self.cidmotion)
        except AttributeError:
            pass

            
            
    ### Class methods
    def shapeSelected(figure,axes):
        '''Called when figure has chosen this shape but not yet drawn it'''
        if hasattr(axes,'format_coord_old'):
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) +' [circle ⊥]'
        else:
            axes.format_coord_old = axes.format_coord
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) + ' [circle ⊥]'
        
    
    def inactive(figure,axes,event):
        radius = 0.01

        figure.newShape = Circle([event.xdata,event.ydata],radius,**cut1DKkwargs) # 
        
        axes.add_patch(figure.newShape)
        
        figure.background = axes.get_figure().canvas.copy_from_bbox(axes.bbox)

        
        figure.cidmove = None
        figure.drawState = States.WIDTH
        
        def on_move(self,event):
            if event.inaxes:
                
                center = self.newShape.center
                mousePoint = np.array([event.xdata-center[0],event.ydata-center[1]])
                
                
                self.newShape.set_radius(np.linalg.norm(mousePoint))
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw both line and rectangle
                axes.draw_artist(self.newShape)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        return on_move
    
    def width(figure,axes,_):
        figure.drawState = States.INACTIVE

        center,radius = figure.newShape.center,figure.newShape.get_radius()
        
        circ = Circle(center,radius,
                        **cut1DKkwargs)
        
        # Find corresponding function to generated DR
        func = figure.draggableFunctions[figure.draggableShapes.index(DraggableCircle)]
        Cut1DFunction = func
        
        
        dr = DraggableCircle(circ,figure,Cut1DFunction,figure)
        axes.add_patch(circ)
        figure.new = False
        figure.newShape.set_visible(False)
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None

class DraggableRectanglePerpendicular(DraggableShape):# pragma: no cover
    ### To be used in a QE plane only allowing a perpendicular cut, i.e. Q_perp (integrating over E and Q_para)
    def __init__(self, rect,plottingObject,Cut1DFunction,figure):
        self.rect = rect
        self.press = None
        self.background = None
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        self.figure = figure
        
    def remove(self):
        fig = self.rect.figure
        self.disconnect()
        self.rect.set_visible(False)
        self.rect.remove()
        if self.selected:
            self.figure.selectedDr = None
        #idx = self.figure.shapes.index(self)
        #self.figure.shapes[idx].remove()
        #del self.figure.shapes[idx]
        if not hasattr(fig,'canvas'): # This is an axis an not a figure!
            fig.get_figure().canvas.draw()        
        else:
            fig.canvas.draw()    

    @property
    def selected(self):
        return self._selected
    
    @selected.getter
    def selected(self):
        return self._selected
    
    @selected.setter
    def selected(self,value):
        if self._selected == bool(value):
            return
        if value:
            if not self.figure.selectedDr is None:
                self.figure.selectedDr.selected = False
            self.figure.selectedDr = self
            color = selectedColor
        else:
            color = deselectedColor
            
        self.rect.set_edgecolor(color)

        self.rect.figure.canvas.draw()
        self._selected = bool(value)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)


    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes: return
        if self.figure.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        
        figure = self.plottingObject
        if event.button == MouseButton.RIGHT and figure.drawState == States.MOVE:
            self.remove()
            return 

        self.selected = True
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        self.figure.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.figure.lock is not self:
            return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.rect.set_xy([x0+dx,y0+dy])

        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        # restore the background region
        canvas.restore_region(self.background)


        # redraw line and rect with the current rectangle
        axes.draw_artist(self.rect)
        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if self.figure.lock is not self:
            return

        self.press = None
        self.figure.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()

    def disconnect(self):
        pass

            
            
    ### Class methods
    def shapeSelected(figure,axes):
        '''Called when figure has chosen this shape but not yet drawn it'''
        if hasattr(axes,'format_coord_old'):
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) +' [rect ⊥]'
        else:
            axes.format_coord_old = axes.format_coord
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) + ' [rect ⊥]'
        
    
    def inactive(figure,axes,event):
        width = height = 0.01

        center = [event.xdata,event.ydata]
        lowerLeft = [center[0]-0.5*width, center[1]-0.5*height]
        figure.newShape = Rectangle(lowerLeft,width=width,height=height,**cut1DKkwargs) # 
        
        axes.add_patch(figure.newShape)
        
        figure.background = axes.get_figure().canvas.copy_from_bbox(axes.bbox)

        
        figure.cidmove = None
        figure.drawState = States.WIDTH
        
        def on_move(self,event):
            if event.inaxes:
                
                width = self.newShape.get_width()
                height = self.newShape.get_height()
                center = np.asarray(self.newShape.xy)+np.array([0.5*width,0.5*height])
                
                mousePoint = np.array([event.xdata-center[0],event.ydata-center[1]])
                
                
                newWidth, newHeight = np.abs(mousePoint)*2.0
                lowerLeft = [center[0]-0.5*newWidth, center[1]-0.5*newHeight]

                self.newShape.set_width(newWidth)
                self.newShape.set_height(newHeight)
                self.newShape.set_xy(lowerLeft)
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw both line and rectangle
                axes.draw_artist(self.newShape)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        return on_move
    
    def width(figure,axes,_):
        figure.drawState = States.INACTIVE

        xy,width,height = figure.newShape.xy,figure.newShape.get_width(),figure.newShape.get_height()
        
        rect = Rectangle(xy,width = width, height=height,
                        **cut1DKkwargs)
        
        # Find corresponding function to generated DR
        func = figure.draggableFunctions[figure.draggableShapes.index(DraggableRectanglePerpendicular)]
        Cut1DFunction = func
        
        axes.add_patch(rect)
        dr = DraggableRectanglePerpendicular(rect,figure,Cut1DFunction,figure)
        figure.new = False
        figure.newShape.set_visible(False)
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None


class DraggableRectangleHorizontal(DraggableShape):# pragma: no cover
    ### To be used in a QE plane only allowing a QCut for constant energy
    def __init__(self, rect,line,plottingObject,Cut1DFunction,figure):
        self.rect = rect
        self.press = None
        self.background = None
        self.line = line
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        self.figure = figure
        
    def remove(self):
        fig = self.rect.figure
        self.disconnect()
        fig = self.rect.figure
        self.rect.set_visible(False)
        self.rect.remove()
        self.line.set_visible(False)
        self.line.remove()
        if self.selected:
            self.figure.selectedDr = None
        #idx = self.figure.shapes.index(self)
        #try:
        #    self.figure.shapes[idx].remove()
        #except ValueError:
        #    pass
        #del self.figure.shapes[idx]
        fig.canvas.draw()     

    @property
    def selected(self):
        return self._selected
    
    @selected.getter
    def selected(self):
        return self._selected
    
    @selected.setter
    def selected(self,value):
        if self._selected == bool(value):
            return
        if value:
            if not self.figure.selectedDr is None:
                self.figure.selectedDr.selected = False
            self.figure.selectedDr = self
            color = selectedColor
        else:
            color = deselectedColor
            
        self.rect.set_edgecolor(color)
        self.line.set_color(color)

        self.rect.figure.canvas.draw()
        self._selected = bool(value)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)


    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes: return
        if self.figure.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        
        figure = self.plottingObject
        if event.button == MouseButton.RIGHT and figure.drawState == States.MOVE:
            self.remove()
            return 

        self.selected = True
        x0, y0 = self.rect.xy
        lineXData = self.line.get_xdata()
        lineYData = self.line.get_ydata()
        self.press = x0, y0, event.xdata, event.ydata, lineXData, lineYData
        self.figure.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        self.line.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.figure.lock is not self:
            return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress, lineXData, lineYData = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.rect.set_xy([x0+dx,y0+dy])
        
        self.line.set_xdata(lineXData+dx)
        self.line.set_ydata(lineYData+dy)

        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        # restore the background region
        canvas.restore_region(self.background)


        # redraw line and rect with the current rectangle
        axes.draw_artist(self.rect)
        axes.draw_artist(self.line)
        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if self.figure.lock is not self:
            return

        self.press = None
        self.figure.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.line.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()

    def disconnect(self):
        pass



            
            
    ### Class methods
    def shapeSelected(figure,axes):
        '''Called when figure has chosen this shape but not yet drawn it'''
        if hasattr(axes,'format_coord_old'):
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) +' [rect constant E]'
        else:
            axes.format_coord_old = axes.format_coord
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) + ' [rect constant E]'
        
    
    def inactive(figure,axes,event):
        width = height = 0.0

        center = [event.xdata,event.ydata]
        lowerLeft = [center[0]-0.5*width, center[1]-0.5*height]
        figure.newShape = Rectangle(lowerLeft,width=width,height=height,**cut1DKkwargs) # 
        
        axes.add_patch(figure.newShape)

        
        
        figure.background = axes.get_figure().canvas.copy_from_bbox(axes.bbox)

        
        figure.cidmove = None
        figure.drawState = States.DIRECTION
        
        
        def on_move(self,event):
            if event.inaxes:
                
                
                # New width is given by distance from mouse to the left corner
                newWidth = event.xdata-self.newShape.get_x()
                

                self.newShape.set_width(newWidth)
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw both line and rectangle
                axes.draw_artist(self.newShape)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        return on_move
    
    def direction(figure,axes,event):
        
        newWidth = event.xdata-figure.newShape.get_x()
        figure.newShape.set_width(newWidth)

        center = figure.newShape.get_xy()+np.array([0.5*figure.newShape.get_height(),0.0])
        dx = figure.newShape.get_width()
        dy = 0.0
        
        figure.line = axes.plot([center[0],center[0]+dx],[center[1],center[1]+dy],cut1DKkwargs['edgecolor'],zorder=cut1DKkwargs['zorder']+1)[0]


        figure.cidmove = None

        
        figure.drawState = States.WIDTH

        
        def on_move(self,event):
            if event.inaxes:
                
                lowerLeft = self.newShape.get_y()
                center = lowerLeft+self.newShape.get_height()*0.5
                # New width is given by distance from mouse to the left corner
                newHeight = np.abs(event.ydata-center)*2.0
                

                self.newShape.set_height(newHeight)
                self.newShape.set_y(center-0.5*newHeight)
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw rectangle
                axes.draw_artist(self.newShape)
                axes.draw_artist(self.line)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        return on_move

    def width(figure,axes,event):
        figure.drawState = States.INACTIVE
        xy,width,height = figure.newShape.xy,figure.newShape.get_width(),figure.newShape.get_height()
        
        rect = Rectangle(xy,width = width, height=height,
                        **cut1DKkwargs)
        
        # Find corresponding function to generated DR
        func = figure.draggableFunctions[figure.draggableShapes.index(DraggableRectangleHorizontal)]
        Cut1DFunction = func
        
        axes.add_patch(rect)
        dr = DraggableRectangleHorizontal(rect,figure.line,figure,Cut1DFunction,figure)
        figure.new = False
        figure.newShape.set_visible(False)
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None


class DraggableRectangleVertical(DraggableShape):# pragma: no cover
    ### To be used in a QE plane only allowing a QCut for constant q
    def __init__(self, rect,line,plottingObject,Cut1DFunction,figure):
        self.rect = rect
        self.press = None
        self.background = None
        self.line = line
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        self.figure = figure
        
    def remove(self):
        fig = self.rect.figure
        self.disconnect()
        fig = self.rect.figure
        self.rect.set_visible(False)
        self.rect.remove()
        self.line.set_visible(False)
        self.line.remove()
        if self.selected:
            self.figure.selectedDr = None
        #idx = self.figure.shapes.index(self)
        #self.figure.shapes[idx].remove()
        #del self.figure.shapes[idx]
        if not hasattr(fig,'canvas'):
            fig.get_figure().canvas.draw()     
        else:
            fig.canvas.draw()     

    @property
    def selected(self):
        return self._selected
    
    @selected.getter
    def selected(self):
        return self._selected
    
    @selected.setter
    def selected(self,value):
        if self._selected == bool(value):
            return
        if value:
            if not self.figure.selectedDr is None:
                self.figure.selectedDr.selected = False
            self.figure.selectedDr = self
            color = selectedColor
        else:
            color = deselectedColor
            
        self.rect.set_edgecolor(color)
        self.line.set_color(color)

        self.rect.figure.canvas.draw()
        self._selected = bool(value)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)


    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes: return
        if self.figure.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        
        figure = self.plottingObject
        if event.button == MouseButton.RIGHT and figure.drawState == States.MOVE:
            self.remove()
            return 

        self.selected = True
        x0, y0 = self.rect.xy
        lineXData = self.line.get_xdata()
        lineYData = self.line.get_ydata()
        self.press = x0, y0, event.xdata, event.ydata, lineXData, lineYData
        self.figure.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        self.rect.set_animated(True)
        self.line.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.rect.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.rect)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.figure.lock is not self:
            return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress, lineXData, lineYData = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.rect.set_xy([x0+dx,y0+dy])
        
        self.line.set_xdata(lineXData+dx)
        self.line.set_ydata(lineYData+dy)

        canvas = self.rect.figure.canvas
        axes = self.rect.axes
        # restore the background region
        canvas.restore_region(self.background)


        # redraw line and rect with the current rectangle
        axes.draw_artist(self.rect)
        axes.draw_artist(self.line)
        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        if self.figure.lock is not self:
            return

        self.press = None
        self.figure.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.line.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()

    def disconnect(self):
        pass

            
            
    ### Class methods
    def shapeSelected(figure,axes):
        '''Called when figure has chosen this shape but not yet drawn it'''
        if hasattr(axes,'format_coord_old'):
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) +' [rect constant Q]'
        else:
            axes.format_coord_old = axes.format_coord
            axes.format_coord = lambda x,y: axes.format_coord_old(x,y) + ' [rect constant Q]'
        
    
    def inactive(figure,axes,event):
        width = height = 0.0

        center = [event.xdata,event.ydata]
        lowerLeft = [center[0]-0.5*width, center[1]-0.5*height]
        figure.newShape = Rectangle(lowerLeft,width=width,height=height,**cut1DKkwargs) # 
        
        axes.add_patch(figure.newShape)

        figure.background = axes.get_figure().canvas.copy_from_bbox(axes.bbox)

        
        figure.cidmove = None
        figure.drawState = States.DIRECTION
        
        def on_move(self,event):
            if event.inaxes:
                
                
                # New width is given by distance from mouse to the left corner
                newHeight = event.ydata-self.newShape.get_y()
                self.newShape.set_height(newHeight)
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw both line and rectangle
                axes.draw_artist(self.newShape)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        return on_move
    
    def direction(figure,axes,event):
        
        newHeight = event.ydata-figure.newShape.get_y()
        figure.newShape.set_height(newHeight)

        center = figure.newShape.get_xy()+np.array([0.5*figure.newShape.get_width(),0.0])
        dx = 0.0
        dy = newHeight
        
        figure.line = axes.plot([center[0],center[0]+dx],[center[1],center[1]+dy],cut1DKkwargs['edgecolor'],zorder=cut1DKkwargs['zorder']+1)[0]

        figure.cidmove = None

        
        figure.drawState = States.WIDTH
        
        def on_move(self,event):
            if event.inaxes:
                
                lowerLeft = self.newShape.get_x()
                center = lowerLeft+self.newShape.get_width()*0.5
                # New width is given by distance from mouse to the left corner
                newWidth = np.abs(event.xdata-center)*2.0
                

                self.newShape.set_width(newWidth)
                self.newShape.set_x(center-0.5*newWidth)
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw rectangle
                axes.draw_artist(self.newShape)
                axes.draw_artist(self.line)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        return on_move

    def width(figure,axes,event):
        figure.drawState = States.INACTIVE

        xy,width,height = figure.newShape.xy,figure.newShape.get_width(),figure.newShape.get_height()
        
        rect = Rectangle(xy,width = width, height=height,
                        **cut1DKkwargs)
        
        # Find corresponding function to generated DR
        func = figure.draggableFunctions[figure.draggableShapes.index(DraggableRectangleVertical)]
        Cut1DFunction = func
        
        axes.add_patch(rect)
        dr = DraggableRectangleVertical(rect,figure.line,figure,Cut1DFunction,figure)
        figure.new = False
        figure.newShape.set_visible(False)
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None