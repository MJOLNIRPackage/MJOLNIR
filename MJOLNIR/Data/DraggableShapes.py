import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib.patches import Rectangle,Circle
import numpy as np

#from MJOLNIR.Data import Viewer3D
from MJOLNIR._interactiveSettings import cut1DSettings,cut1DSettingsAll,selectedColor, \
deselectedColor, States, cut1DKkwargs

# draggable rectangle with the animation blit techniques; see
# http://www.scipy.org/Cookbook/Matplotlib/Animations


class DraggableShape():
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
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)


            
class DraggableRectangle(DraggableShape):
    
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
        self.rect.remove()
        self.line.remove()
        if self.selected:
            self.figure.selectedDr = None
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
        if event.button == MouseButton.RIGHT and figure.drawState == States.INACTIVE:
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
        # TODO: check for 3DViewer#print('here')
        #if hasattr(figure,'axes'): # only Viewer3D
        #    if figure.axis != 2:
        #        return figure,axes,event):

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
        figure.shapes.append(rect)
        figure.drawState = States.INACTIVE
       
        

        figure.drs.append(dr)

        return None

def corner2center(x,y,width,angle):
    # Rotate from rectangle space to real space
    angleRad = np.deg2rad(angle+90)
    # Unit vector along height of rectangle
    dirvec = np.array([np.cos(angleRad),np.sin(angleRad)])
    # Orthogonal with half width
    ortho = np.array([dirvec[1],-dirvec[0]])*width*0.5
    # center point is orthogonal to height and half the width
    center = np.array([x,y])+ortho
    return center

def center2corner(x,y,width,angle):
    # Rotate from rectangle space to real space
    angleRad = np.deg2rad(angle+90)
    # Unit vector along height of rectangle
    dirvec = np.array([np.cos(angleRad),np.sin(angleRad)])
    # Orthogonal with half width
    ortho = np.array([dirvec[1],-dirvec[0]])*width*0.5
    # corner point is orthogonal to height and half the width but reversed
    corner = np.array([x,y])-ortho
    return corner

def extractCut1DPropertiesRectangle(rect,sample=None, rounding = 4):
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


def extractCut1DPropertiesRectanglePerpendicular(rect,sample=None, rounding = 4):
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

def extractCut1DPropertiesRectangleHorizontal(rect,sample=None, rounding = 4):
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

def extractCut1DPropertiesRectangleVertical(rect,sample=None, rounding = 4):
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

def extractCut1DPropertiesCircle(circ,sample=None, rounding = 4):
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

def clearBoxes(self):
    """Clear all generated draggable rectangles"""
    if hasattr(self,'ax'):
        axes = self.ax
    else:
        axes = self

    for dr in self.drs:
        dr.remove()
        
    self.drs = []
    axes.get_figure().canvas.draw()


class DraggableCircle(DraggableShape):
    
    def __init__(self, circ,plottingObject,Cut1DFunction,figure):
        self.circ = circ
        self.press = None
        self.background = None
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        self.figure = figure
        
    def remove(self):
        fig = self.circ.figure
        self.disconnect()
        self.circ.remove()
        if self.selected:
            self.figure.selectedDr = None
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
        if event.button == MouseButton.RIGHT and figure.drawState == States.INACTIVE:
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
        self.circ.figure.canvas.mpl_disconnect(self.cidpress)
        self.circ.figure.canvas.mpl_disconnect(self.cidrelease)
        self.circ.figure.canvas.mpl_disconnect(self.cidmotion)

            
            
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
        
        axes.add_patch(circ)
        dr = DraggableCircle(circ,figure,Cut1DFunction,figure)
        figure.new = False
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None

class DraggableRectanglePerpendicular(DraggableShape):
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
        self.rect.remove()
        if self.selected:
            self.figure.selectedDr = None
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
        if event.button == MouseButton.RIGHT and figure.drawState == States.INACTIVE:
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
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)

            
            
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
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None


class DraggableRectangleHorizontal(DraggableShape):
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
        self.rect.remove()
        self.line.remove()
        if self.selected:
            self.figure.selectedDr = None
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
        if event.button == MouseButton.RIGHT and figure.drawState == States.INACTIVE:
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
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)

            
            
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
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None


class DraggableRectangleVertical(DraggableShape):
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
        self.rect.remove()
        self.line.remove()
        if self.selected:
            self.figure.selectedDr = None
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
        if event.button == MouseButton.RIGHT and figure.drawState == States.INACTIVE:
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
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)

            
            
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
        figure.newShape.remove()
        figure.newShape = None

        dr.selected = True
        dr.connect()
        figure.shapes.append(dr)
        
        return None


def cancel(self,axes):
    self.resetFormatCoord()
    self.newShape = None
    self.new = False
    
    # Only when something has been drawn is cid set
    if not self.cidmove is None:
        axes.get_figure().canvas.mpl_disconnect(self.cidmove)
        self.cidmove = None
        
        del self.patches[-1]
        plt.draw()  

    # Reset state and return
    self.drawState = States.INACTIVE
    return 



def on_key_press(self,event):
    if not hasattr(self,'drawState'):
        return
    if not event.key in cut1DSettingsAll:
        self.resetFormatCoord()
        self.suppressPrint = False

        return
    if event.key in cut1DSettings['move']:
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
        

        return
    elif event.key in cut1DSettings['cut'] and self.new is False:
        self.lock=True
        self.resetFormatCoord()
        dr = self.selectedDr
        if not dr is None:
            dr.Cut1DFunction(dr)

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
        else:
            self.resetFormatCoord()
        
        self.lock = True
    


def on_press(self,event):
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
        
        
def reset(self):
    if not hasattr(self,'drawState'):
        return
    self.suppressPrint=False
    if hasattr(self,'format_coord_old'):
        self.format_coord = self.format_coord_old
        del self.format_coord_old


def prepareInteractiveCutting(ax,draggables, draggableFunctions):
    fig = ax.get_figure()
    ax.resetFormatCoord = lambda: reset(ax)
    ax.new=False
    ax.cidmove = None
    ax.drawState = States.INACTIVE
    
    ax.shapes = []
    ax.drs = []
    ax.draggableShapes = draggables
    ax.draggableFunctions = draggableFunctions
    ax.ax = ax
    
    ax.selectedDr = None
    ax.lock = True
    
    ax.suppressPrint = False
    
    fig.canvas.mpl_connect(
        'key_press_event', lambda event:on_key_press(ax,event))
    
    fig.canvas.mpl_connect(
                'button_press_event', lambda event:on_press(ax,event))


def prepareInteractiveCuttingView3D(self,draggables, draggableFunctions):
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
