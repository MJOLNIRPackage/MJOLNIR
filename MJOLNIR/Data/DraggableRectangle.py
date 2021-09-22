import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib.patches import Rectangle
import numpy as np

from MJOLNIR.Data import Viewer3D
from MJOLNIR._interactiveSettings import cut1DSettings,cut1DSettingsAll,selectedColor, \
deselectedColor, States, cut1DKkwargs

# draggable rectangle with the animation blit techniques; see
# http://www.scipy.org/Cookbook/Matplotlib/Animations


class DraggableRectangle():
    lock = None  # only one can be moved at a time
    selectedDR = None # Only one can be selected at one point
    cidkeypress = None
    
    def __init__(self, rect,line,plottingObject,Cut1DFunction):
        self.rect = rect
        self.line = line
        self.press = None
        self.background = None
        self._selected = False
        self.plottingObject = plottingObject
        self.Cut1DFunction = Cut1DFunction
        
    def remove(self):
        self.rect.remove()
        self.line.remove()
        if self.selected:
            DraggableRectangle.selectedDR = None
        plt.draw()            

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
            if not DraggableRectangle.selectedDR is None:
                DraggableRectangle.selectedDR.selected = False
            DraggableRectangle.selectedDR = self
            color = selectedColor
        else:
            color = deselectedColor
            
        self.rect.set_edgecolor(color)
        self.line.set_color(color)
        plt.draw()
        self._selected = bool(value)

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        if DraggableRectangle.cidkeypress is None:
            DraggableRectangle.cidkeypress = self.rect.figure.canvas.mpl_connect(
                'key_press_event', self.on_key_press)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes: return
        if DraggableRectangle.lock is not None: return
        contains, attrd = self.rect.contains(event)
        if not contains: return
        self.selected = True
        figure = self.plottingObject
        if event.button == MouseButton.RIGHT and figure.drawState == States.INACTIVE:
            self.remove()
            return 
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata
        (DraggableRectangle).lock = self

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
        if DraggableRectangle.lock is not self:
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
        if DraggableRectangle.lock is not self:
            return

        self.press = None
        DraggableRectangle.lock = None

        # turn off the rect animation property and reset the background
        self.rect.set_animated(False)
        self.line.set_animated(False)
        self.background = None

        # redraw the full figure
        self.rect.figure.canvas.draw()
        
    def on_key_press(self,event):
       
        if not event.key in cut1DSettingsAll:
            return

        if event.key in cut1DSettings['move']:
            self.toggle_lock()
        
    def toggle_lock(self):
        if DraggableRectangle.lock is None:
            DraggableRectangle.lock = True
        else:
            DraggableRectangle.lock = None

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)
        if not DraggableRectangle.cidkeypress is None:
            self.rect.figure.canvas.mpl_disconnect(self.cidkeypress)
            DraggableRectangle.cidkeypress = None


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

def extractCut1DProperties(rect,sample=None, rounding = 4):
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

def on_key_press(self,event):
    if not event.key in cut1DSettingsAll:
        return
    if event.key in cut1DSettings['move']:
        self.new = False
        return
    if event.key in cut1DSettings['cut']:
        dr = DraggableRectangle.selectedDR
        if not DraggableRectangle.selectedDR is None:
            dr.Cut1DFunction(dr)

            return
    self.new = True
    DraggableRectangle.lock = True
    

def on_press(self,event):
    if hasattr(self,'ax'):
        axes = self.ax
    else:
        axes = self

    if not self.new is True:
        return
    if event.inaxes != axes:
        return

   
    # Right-clicking cancels any ongoing action and reset to inactive
    if event.button == MouseButton.RIGHT:
        self.newRect = None
        self.new = False
        
        # Only when something has been drawn is cid set
        if not self.cidmove is None:
            axes.get_figure().canvas.mpl_disconnect(self.cidmove)
            self.cidmove = None
            
            del self.patches[-1]
            plt.draw()  
            
        if not self.line is None:
            self.line.remove()
        # Reset state and return
        self.drawState = States.INACTIVE
        return 


    if self.drawState == States.INACTIVE:
        if isinstance(self,Viewer3D.Viewer3D):
            if self.axis != 2:
                return 
        width = 0.0
        height = 0.05
        center = corner2center(event.xdata,event.ydata, width, -180.0)
        self.newRect = Rectangle(center,width,height,**cut1DKkwargs) # 
        
        axes.add_patch(self.newRect)
        
        axes.get_figure().canvas.draw()
        def on_move(self,event):
            if event.inaxes:
                rectx,recty = self.newRect.xy
                center = corner2center(rectx,recty, width, self.newRect.angle)
                mousex = event.xdata
                mousey = event.ydata
                dx,dy = center[0]-mousex, center[1]-mousey
                angle = np.arctan2(dy,dx)+np.pi/2
                
                self.newRect.set_xy(center2corner(*center, width, np.rad2deg(angle)))
                self.newRect.angle = np.rad2deg(angle)
                
                self.newRect.set_height(np.linalg.norm([dx,dy]))
                
                
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw just the current rectangle
                axes.draw_artist(self.newRect)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        self.newRect.set_animated(True)
        axes.get_figure().canvas.draw()
        self.background = self.figure.canvas.copy_from_bbox(axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.newRect)

        # and blit just the redrawn area
        axes.get_figure().canvas.blit(axes.bbox)
        self.cidmove = axes.get_figure().canvas.mpl_connect('motion_notify_event',lambda event:on_move(self,event))
        self.drawState = States.INITIAL
    elif self.drawState == States.INITIAL:
        axes.get_figure().canvas.mpl_disconnect(self.cidmove)
        
        # Find position of center line
        center = corner2center(*self.newRect.get_xy(), self.newRect.get_width(), self.newRect.angle)

        height = self.newRect.get_height()
        angle = np.deg2rad(self.newRect.angle+90)
        dx,dy = height*np.array([np.cos(angle),np.sin(angle)])
        
        
        self.line = axes.plot([center[0],center[0]+dx],[center[1],center[1]+dy],cut1DKkwargs['edgecolor'],zorder=cut1DKkwargs['zorder']+1)[0]
        
        
        
        self.cidmove = None
        self.drawState = States.WIDTH
        
        def on_move(self,event):
            if event.inaxes:
                
                rectx,recty = self.newRect.xy
                angle = np.deg2rad(self.newRect.angle+90)
                dirvec = np.array([np.cos(angle),np.sin(angle)])
                center = corner2center(rectx,recty, self.newRect.get_width(), self.newRect.angle)
                mousePoint = np.array([event.xdata-center[0],event.ydata-center[1]])
                
                ortho = mousePoint-np.dot(dirvec,mousePoint)*dirvec
                sign = np.sign(np.dot(ortho,np.array([dirvec[1],-dirvec[0]])))
                
                width = sign*np.linalg.norm(ortho)*2.0
                
                
                self.newRect.set_width(width)
                self.newRect.set_xy(center2corner(*center, width, self.newRect.angle))
                
                canvas = axes.get_figure().canvas
                # restore the background region
                canvas.restore_region(self.background)
        
                # redraw both line and rectangle
                axes.draw_artist(self.newRect)
                axes.draw_artist(self.line)
        
                # blit just the redrawn area
                canvas.blit(axes.bbox)
                
        axes.get_figure().canvas.blit(axes.bbox)
        
        self.cidmove = axes.get_figure().canvas.mpl_connect('motion_notify_event',lambda event:on_move(self,event))
        
    elif self.drawState == States.WIDTH:

        axes.get_figure().canvas.mpl_disconnect(self.cidmove)
        self.cidmove = None
        self.drawState = States.INACTIVE

        rect = Rectangle(self.newRect.get_xy(),width=self.newRect.get_width(),
                         height=self.newRect.get_height(),angle=self.newRect.angle,
                         **cut1DKkwargs)
        
        
        Cut1DFunction = lambda dr: self.cut1DFunction(self,dr)
        
        axes.add_patch(rect)
        dr = DraggableRectangle(rect,self.line,self,Cut1DFunction)
        dr.connect()
        dr.selected = True
        self.drs.append(dr)
        self.new = False
        self.newRect = None
        self.rects.append(rect)
        plt.draw()
        