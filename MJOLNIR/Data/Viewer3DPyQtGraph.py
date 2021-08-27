import warnings
from pyqtgraph.Qt import QtGui, QtCore

import numpy as np
import pyqtgraph as pg
from MJOLNIR import _tools
import matplotlib.pyplot as plt

class Interactive3DViewer(QtGui.QWidget):
    def __init__(self,data,bins,sample,log=False,*args,**kwargs):
        annoyingFigure = plt.gca().get_figure()
        if not len(annoyingFigure.get_children()) <= 1: # If less than or equal to 1 child figure is empty
            plt.close(annoyingFigure) # Mega hack to save the day (close annoying mpl figure....... not a fan). MDT = magic don't touch
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        super().__init__(*args,**kwargs)
        self.l = QtGui.QGridLayout()
        self.setLayout(self.l)
        
        # Add plotItem to allow for axes
         
        self.imv1 = CustomImageView(view=pg.PlotItem())
        self.imv2 = pg.ImageView(view=pg.PlotItem())
        
        #Insert widgets into layout
        self.l.addWidget(self.imv1, 0, 0)
        self.l.addWidget(self.imv2, 1, 0)
        
        
        # Create the region-of-interest line segment
        self.roi = pg.ROI([1, 0], [1,1], pen='r', resizable=True)
        self.roi.setSize([1,0.1])
        self.roi.addScaleHandle([0.5, 1], [0.5, 0.5],)
        self.roi.addScaleHandle([0, 0.5], [0.5, 0.5])
        self.roi.addRotateHandle([1.0, 1.0], [0.5, 0.5])

        # Change color of roi markers
        handleColor = QtGui.QColor(255,102,0)
        for handle in self.roi.getHandles():
            handle.pen.setColor(handleColor)


        self.imv1.addItem(self.roi)
        
        
        self.Data,self.bins = data,bins
        
        if len(data)==4:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.data = np.divide(self.Data[0]*self.Data[-1],self.Data[1]*self.Data[2])
                self.data[self.Data[-1]==0] = -np.nanmin(self.data[self.data!=0])
        else:
            self.data = self.Data

        # set unmeasured areas to negative value
        if log:
            self.data = np.log10(self.data+1e-20)
        
        
        # transpose data to comply with viewer requirements
        self.data = self.data.transpose((2,0,1))
        
        # Extract the sample
        self.sample = sample
        # Calculate the projection directions
        self.axes = np.array([self.sample.calculateQxQyToHKL(*X) for X in self.sample.RotMat])
        
        
        
        ## Display the data
        self.imv1.setImage(self.data,xvals=self.bins[2][0,0,:],levels=(-1e-7, 1e-5),autoLevels=False)
        
        
        
        self.imv1.setHistogramRange(-8.5e-8, 1e-5)
        # unlock aspect ratio
        self.imv1.view.setAspectLocked(False)
        
        # Generate and add correct labels to axes
        self.xaxis = self.imv1.view.axes['bottom']['item']
        self.yaxis = self.imv1.view.axes['left']['item']
        
        self.xlabel = _tools.generateLabel(_tools.LengthOrder(self.axes[0]))
        self.ylabel = _tools.generateLabel(_tools.LengthOrder(self.axes[1]))
        
        self.xaxis.setLabel(text=self.xlabel,units='RLU')
        self.yaxis.setLabel(text=self.ylabel,units='RLU')
                   
        
        # Setup color map for main window
        self.setupColorscale(self.imv1)
        
        
        # Calcualte scaling and offsets between image view and RLU axes
        self.qxRange = np.diff(self.bins[0][[0,-1],0,0])
        self.qyRange = np.diff(self.bins[1][0,[0,-1],0])
        self.ERange = np.diff(self.bins[2][0,0,[0,-1]])
        
        self.qxScale = self.qxRange*np.linalg.norm(self.axes[0])/self.bins[0].shape[0]
        self.qyScale = self.qyRange*np.linalg.norm(self.axes[1])/self.bins[0].shape[1]
        self.EScale =  self.ERange/self.bins[0].shape[2]
        
        self.qxCenter = self.bins[0][0,0,0]*np.linalg.norm(self.axes[0])
        self.qyCenter = self.bins[1][0,0,0]*np.linalg.norm(self.axes[1])
        self.ECenter = self.bins[2][0,0,0]
        
        # Apply scaling and translation
        img1Item = self.imv1.getImageItem()
        img1Item.scale(self.qxScale,self.qyScale)
        img1Item.translate(self.qxCenter/self.qxScale,self.qyCenter/self.qyScale)
        # Un-invert yaxis
        self.imv1.view.getViewBox().invertY(False)
        self.imv1.view.autoRange(True)
        self.imv1.view.setAutoVisible(x=True,y=True)
        
        
        # Add labels colormap for cut window
        self.imv2.view.setLabel("left", "Energy", units='meV')
        self.imv2.view.setLabel("bottom", "HKL", units='RLU')
        self.setupColorscale(self.imv2)
        
        # Hide all unneeded menus and uninvert yaxes
        self.imv2.ui.roiBtn.hide()
        self.imv2.ui.menuBtn.hide()
        self.imv2.view.getViewBox().invertY(False)
        
        # Extract projection matrix used for position calculation along cut
        self.projectionMatrix = self.axes
        self.projectionMatrix[0]*=1/np.linalg.norm(self.projectionMatrix[0])
        self.projectionMatrix[1]*=1/np.linalg.norm(self.projectionMatrix[1])
        self.xaxis2 = self.imv2.view.axes['bottom']['item']
        self.yaxis2 = self.imv2.view.axes['left']['item']
        
        
        
        self.xaxis2.tickStrings = lambda values,scale,spacing: self.XtickStrings(self.roi,values,scale,spacing)
        self.yaxis2.tickStrings = lambda values,scale,spacing: self.YtickStrings(values,scale,spacing)
        
        # Create function to be called when cut changes
                
        # Connect update-function to correct slot
        self.roi.sigRegionChanged.connect(self.update)
        
        # Call update for initial position
        self.update()
        
    def XtickStrings(self,roi,values,scale,spacing):
        d2,XY = roi.getArrayRegion(self.data, self.imv1.imageItem, axes=(1,2),returnMappedCoords=True)
    
        d2 = np.sum(d2,axis=-1) # sum together data along axis 1
        XY = np.mean(XY,axis=2)
        
        P1 = XY[0]#*qxScale+qxCenter
        P2 = XY[1]#*qyScale+qyCenter
        P = np.array([P1,P2])
        
        QLength = np.linalg.norm(np.diff(P[:,[0,-1]],axis=1).reshape(-1)/np.array([np.linalg.norm(self.axes[0]),np.linalg.norm(self.axes[1])]))
        QScale = QLength/P.shape[1]
        
        # Calculate scale in Å^-1
        HKL = np.einsum('ji,jk->ik',P,self.projectionMatrix)
        
        direction = HKL[-1]-HKL[0]
        values = (np.array(values)/P.shape[1]).reshape(-1,1)*direction.reshape(1,-1)+HKL[0].reshape(1,-1)
            
        places = max(0, np.ceil(-np.log10(spacing/P.shape[1]*QScale)))
        
        
        fmtStr = '{:.'+str(int(places))+'f}'
        strings = ['\n'.join([fmtStr.format(x) for x in pos]) for pos in values]
    
        
        return strings
        
    def YtickStrings(self,values,scale,spacing):
            
        E = self.imv1.tVals[[0,-1]]
        energies = values*np.diff(E)/len(self.imv1.tVals)+E[0]
        places = max(0, np.ceil(-np.log10(spacing/len(self.imv1.tVals))))
        fmtStr = '{:.'+str(int(places))+'f}'
        strings = [fmtStr.format(x) for x in energies]
        return strings
        
    def setupColorscale(self,imv):
        
        imv.setPredefinedGradient('viridis')
        g = imv.ui.histogram.gradient
        first = list(g.ticks)[0]
        color = QtGui.QColor(first.color)
        first.color.setAlpha(0)
        first.color.setRed(0)
        first.color.setGreen(0)
        first.color.setBlue(0)
        g.addTick(0.01, color)
        
    def update(self):
        
        # Extract cut data and positions
        d2 = self.roi.getArrayRegion(self.data, self.imv1.imageItem, axes=(1,2),returnMappedCoords=False)
        
        d2 = np.sum(d2,axis=-1) # sum together data along axis 1
        
        #QScale=1
        
        # Add data to be shown and correct scale and translation
        self.imv2.setImage(d2.T,autoLevels=True)#,levels=(-5.5e-8, 1e-5))

        
        # Unlock aspect ratio
        self.imv2.view.setAspectLocked(False)
        self.yaxis2.setRange(-1,d2.shape[0]+1)
        self.xaxis2.setRange(-5,d2.shape[1]+5)

    def set_clim(self,vMin,vMax):
        self.imv1.setLevels(vMin,vMax)
        self.imv2.setLevels(vMin,vMax)

    def set_title(self,title):
        self.parent().window().setWindowTitle(title)

        

from pyqtgraph import ptime

class CustomImageView(pg.ImageView):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        

        self.timeLine.hide()
        self.timeLine.sigPositionChanged.disconnect(self.timeLineChanged)
        self.ui.roiPlot.removeItem(self.timeLine)
        self.removeItem(self.timeLine)
        
        self.timeLine2 = pg.LinearRegionItem((0.0,1.0), movable=True)
        self.currentWidth = 1
        self.ui.roiPlot.addItem(self.timeLine2)
        self.ui.roiPlot.setMinimumHeight(60)
        self.ui.roiPlot.hideButtons()
        
        plotItem = self.ui.roiPlot.getPlotItem()
        Eaxis = plotItem.getAxis('bottom')
        
        
        Eaxis.setLabel('Energy',units='meV')
        
        
        self.timeLine2.sigRegionChanged.connect(self.timeLineChanged)
        
        self.ui.roiBtn.hide()
        self.ui.menuBtn.hide()
        
    def setImage(self,*args,**kwargs):
        super().setImage(*args,**kwargs)
        
        
        self.currentIndex = 0
        self.currentWidth = 1
        
        lowerPos = self.tVals[self.currentIndex]
        upperPos = self.tVals[self.currentIndex+self.currentWidth]
        
        
        self.timeLine2.setBounds(self.timeLine.bounds())
        self.timeLine2.setRegion((lowerPos,upperPos))
        
        
        
    def setCurrentIndex(self, ind):
        """Set the currently displayed frame index."""
        index = np.clip(ind, 0, self.getProcessedImage().shape[self.axes['t']]-1)
        self.ignorePlaying = True
        # Implicitly call timeLineChanged
        
        *_, width = self.timeIndex(self.timeLine2)
        self.currentWidth = width
        maxindex = np.min([len(self.tVals)-1,index+width]) # Don't go out of bounce
        self.timeLine2.setRegion((self.tVals[index],self.tVals[maxindex]))
        
        #self.timeLine.setValue(self.tVals[index])
        self.ignorePlaying = False
        
    def timeLineChanged(self):
        if not self.ignorePlaying:
            self.play(0)

        ind, time, width = self.timeIndex(self.timeLine2)
        if ind != self.currentIndex or width != self.currentWidth:
            self.currentIndex = ind
            self.currentWidth = width
            self.updateImage()
        self.sigTimeChanged.emit(ind, time)
        
        
    def timeIndex(self, slider):
        ## Return the time and frame index indicated by a slider
        if self.image is None:
            return [0,0,1]
        
        t = slider.getRegion()

        xv = self.tVals
        if xv is None:
            ind = int(t[0])
            width = int(t[1])-ind
        else:
            if len(xv) < 2:
                return [0,0,1]
            #totTime = xv[-1] + (xv[-1]-xv[-2])
            inds = np.argwhere(xv <= t[0])
            
            if len(inds) < 1:
                return [0,0,1]
            ind = inds[-1,0]
            width = np.argwhere(xv <= t[1])[-1][0]-ind
            width = np.max([1,width])
        return [ind, t, width]
    
    
    def updateImage(self, autoHistogramRange=True):
        ## Redraw image on screen
        if self.image is None:
            return
            
        image = self.getProcessedImage()
        
        if autoHistogramRange:
            self.ui.histogram.setHistogramRange(self.levelMin, self.levelMax)
        
        # Transpose image into order expected by ImageItem
        if self.imageItem.axisOrder == 'col-major':
            axorder = ['t', 'x', 'y', 'c']
        else:
            axorder = ['t', 'y', 'x', 'c']
        axorder = [self.axes[ax] for ax in axorder if self.axes[ax] is not None]
        image = image.transpose(axorder)
            
        # Select time index
        if self.axes['t'] is not None:
            self.ui.roiPlot.show()
            image = np.sum(image[self.currentIndex:self.currentIndex+int(self.currentWidth)],axis=0)
            
        self.imageItem.updateImage(image)

    def evalKeyState(self):
        if len(self.keysPressed) == 1:
            key = list(self.keysPressed.keys())[0]
            if key == QtCore.Qt.Key_Right:
                self.play(20)
                self.jumpFrames(1)
                self.lastPlayTime = ptime.time() + 0.2  ## 2ms wait before start
                                                        ## This happens *after* jumpFrames, since it might take longer than 2ms
            elif key == QtCore.Qt.Key_Left:
                self.play(-20)
                self.jumpFrames(-1)
                self.lastPlayTime = ptime.time() + 0.2
            elif key == QtCore.Qt.Key_Up:
                self.play(100)
            elif key == QtCore.Qt.Key_Down:
                self.play(-100)
            elif key == QtCore.Qt.Key_PageUp:
                self.play(1000)
            elif key == QtCore.Qt.Key_PageDown:
                self.play(-1000)
        else:
            self.play(0)