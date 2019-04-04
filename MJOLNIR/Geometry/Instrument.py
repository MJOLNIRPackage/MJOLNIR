import sys
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import numpy as np
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector,Wedge
from MJOLNIR import _tools
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import scipy.optimize
import h5py as hdf
import datetime
import pytest

NumberOfSigmas= 3 # Defining the active area of a peak on a detector as \pm n*sigma



class Instrument(GeometryConcept.GeometryConcept):
    @_tools.KwargChecker(include=['Author','Instrument','Date','Initialized']) # Not used as excess kwargs are put into settings
    def __init__(self, position=(0,0,0),wedges=[],fileName='',**kwargs):
        """Instrument object used to calculated analytic scattering coverage. 
        Based on the GeometryConcept object it contains all needed information about the setup used in further calculations.

        Kwargs:
            
            - position (float 3d): Position of the instrument always at origin (default (0,0,0))

            - wedges (list of wedges or single wedge): Wedge or list of wedges which the instrument consists of (default empty)

            - fileName (string): Filename of xml file (ending in xml). To load binary files use self.load(filename).

        Raises:
            
            - AttributeError
        
        """
        
        self._wedges = []
        
        self._settings = {}
        if fileName !='':
            if(fileName.split('.')[-1]=='xml'):
                parseXML(self,fileName)
            else:
                raise ValueError('File not of type XML.')
        else:
            super(Instrument,self).__init__(position)
            
            for key in kwargs:
                self.settings[key]=kwargs[key]
            self._settings['Initialized']=False
            self.append(wedges)

    @property
    def wedges(self):
        return self._wedges

    @wedges.getter
    def wedges(self):
        return self._wedges

    @wedges.setter
    def wedges(self,wedges):
        if len(self.wedges)!=0:
            warnings.warn('The list of wedges is not empty! Appending new wedges(s)',stacklevel=2)
        if isinstance(wedges, list):
            for ana in wedges:
                if not issubclass(type(ana),Wedge.Wedge):
                    raise AttributeError('Object is not an wedge or a simple list of these')
                self._wedges.append(ana)
        else:
            if not issubclass(type(wedges),Wedge.Wedge):
                raise AttributeError('Object is not an analyser or a simple list of these')
            self._wedges.append(wedges)
    
    def append(self,wedge):
        """Append wedge(s) to instrument.

        Args
            
            - wedge (Wedge(s)): Single wedge or list of wedges
        """
        if isinstance(wedge,list):
            for obj in wedge:
                if issubclass(type(obj),Wedge.Wedge):
                    self._wedges.append(obj)
                else:
                    raise AttributeError('Object not wedge or a simple list of wedges')
        else:
            if issubclass(type(wedge),Wedge.Wedge):
                    self._wedges.append(wedge)
            else:
                raise AttributeError('Object not wedge or a simple list of wedges')

    def plot(self,ax):
        """Recursive plotting routine."""
        for wedge in self.wedges:
            wedge.plot(ax,offset=self.position)

    @property
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')

    def __str__(self):
        string = '{} with settings:\n'.format(self.__class__)
        for attrib in self.settings:
            string+='{}:\t{}\n'.format(attrib,self.settings[attrib])
        string+='\n'    
        string+='Containing the following wedges:\n'
        for wedge in self.wedges:
            string+=str(wedge)+'\n'
        return string

    def initialize(self):
        """Method to initialize and perform analytical calculations of scattering quantities. 
        Initializes:

            -  A4: Matrix holding pixel A4. Shape (len(Wedges),len(detectors),pixels)
            
            -  Ef: Matrix holding pixel Ef. Shape (len(Wedges),len(detectors),pixels)
        """
        factorLambdasqrtE = 9.0445678

        if(len(self.wedges)==0):
            raise ValueError('Instrument does not contain any wedges and can thus not be initialized.')
        self._A4 = []
        self._Ef = []
        beamDirection = np.array([0.0,1.0,0.0])

        for wedge in self.wedges:
            detectorPixelPositions,analyserPixelPositions = wedge.calculateDetectorAnalyserPositions()

            A4 = [-np.arccos(np.divide(np.dot(AnalyserPos,beamDirection),
                np.linalg.norm(AnalyserPos,axis=1)))*np.sign(np.cross(AnalyserPos,beamDirection)[:,-1]) for AnalyserPos in analyserPixelPositions]

                
            relPos = [detectorPixelPositions[i]-analyserPixelPositions[i] for i in range(len(analyserPixelPositions))]

            A6 = [np.arccos(np.divide(np.einsum('ij,ij->i',analyserPixelPositions[i],relPos[i]),
                np.linalg.norm(analyserPixelPositions[i],axis=1)*np.linalg.norm(relPos[i],axis=1))) for i in range(len(analyserPixelPositions))]
            Ef = [np.power(factorLambdasqrtE/(wedge.analysers[0].d_spacing*2.0*np.sin(A6Sub/2.0)),2.0) for A6Sub in A6] 
            self._A4.append(A4)
            self._Ef.append(Ef)





        self.settings['Initialized']=True

    @property
    def A4(self):
        return self._A4
    @A4.getter
    def A4(self):
        if(self.settings['Initialized']==False):
            raise RuntimeError('Instrument is not initialized.')
        return self._A4
    @A4.setter
    def A4(self,*args,**kwargs):
        raise NotImplementedError('A4 cannot be overwritten.')

    @property
    def Ef(self):
        return self._Ef
    @Ef.getter
    def Ef(self):
        if(self.settings['Initialized']==False):
            raise RuntimeError('Instrument is not initialized.')
        return self._Ef
    @Ef.setter
    def Ef(self,*args,**kwargs):
        raise NotImplementedError('Ef cannot be overwritten.')

    
    def saveXML(self,fileName):
        """Method for saving current file as XML in fileName."""
        XMLString = '<?xml version="1.0"?>\n'
        XMLString+= '<Instrument '
        for attrib in self.settings:
            XMLString+="{}='{}' ".format(attrib,self.settings[attrib])

        XMLString+="position='"+','.join([str(x) for x in self.position])+"'"
        XMLString+='>\n'
            
        for wedge in self.wedges:
            XMLString+="\t<Wedge "
            for attrib in wedge.settings:
                XMLString+="{}='{}' ".format(attrib,wedge.settings[attrib])

            XMLString+="position='"+','.join([str(x) for x in wedge.position])+"'"
            XMLString+='>\n'

            for item in wedge.analysers + wedge.detectors:
                itemClass = str(item.__class__).split('.')[-1][:-2]
                XMLString+="\t\t<{}".format(itemClass)
                for key in item.__dict__:
                    value = item.__getattribute__(key)
                    if isinstance(value,type(np.array([0,0,0]))):
                        valueStr = ','.join([str(x) for x in item.__getattribute__(key)])
                    else:
                        valueStr = str(value)
                    XMLString+=" {}='{}'".format(str(key)[1:],valueStr)
                XMLString+="></{}>\n".format(itemClass)
                
            
            XMLString+="\t</Wedge>\n"
        XMLString+="</Instrument>\n"
    
        f = open(fileName,'w')
        f.write(XMLString)
        f.close()

    def generateCAMEAXML(self,fileName):
        """Generate CAMEA XML file to be used as instrument file.

        Args:

            - fileName: Name of file to be saved (required)

        """
        ang_1 = np.zeros((7,))
        ang_2 = np.zeros((6,))

        ang_1[6]=-3.3#3
        ang_1[5]=-2.2#2
        ang_1[4]=-1.1#1
        ang_1[3]=0
        ang_1[2]=1.1#1
        ang_1[1]=2.2#2
        ang_1[0]=3.3#3

        ang_2[5]=-2.75#75
        ang_2[4]=-1.65#65
        ang_2[3]=-0.55#5
        ang_2[2]=0.55#5
        ang_2[1]=1.65#65
        ang_2[0]=2.75#75


        z_an = np.zeros((8,))
        z_an[0]=0.9300
        z_an[1]=0.9939
        z_an[2]=1.0569
        z_an[3]=1.1195
        z_an[4]=1.1827
        z_an[5]=1.2456
        z_an[6]=1.3098
        z_an[7]=1.3747
    
        

        H1 = 0.7
        H2 = 0.71

        det_cen = 1.2
        wedges=8

        offset =0.0# -4.835960288880082# offset such that last pixel of detector 0 is at 0


        string = "<?xml version='1.0'?>\n<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>\n"
        for W in -np.arange(wedges):
            
            string+="\t<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>\n"
            
            Anaposx = -np.sin((W*8.0+offset)*np.pi/180)*z_an
            Anaposy = np.cos((W*8.0+offset)*np.pi/180)*z_an
            
            for i in range(len(z_an)):
                XX = Anaposx[i]/(np.sqrt(2)*z_an[i])
                YY = Anaposy[i]/(np.sqrt(2)*z_an[i])
                string+="\t\t<FlatAnalyser position='"+str(Anaposx[i])+','+str(Anaposy[i])+",0.0' direction='"+str(YY)+","+str(XX)+",0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>\n"
            
            
            detx_1 = -np.sin((ang_1+W*8.0+offset)*np.pi/180)*det_cen
            detz_1 = np.cos((ang_1+W*8.0+offset)*np.pi/180)*det_cen
            
            
            detx_2 = -np.sin((ang_2+W*8.0+offset)*np.pi/180)*det_cen
            detz_2 = np.cos((ang_2+W*8.0+offset)*np.pi/180)*det_cen
            for i in range(7):
                string+="\t\t<TubeDetector1D position='"+str(detx_1[i])+','+str(detz_1[i])+','+str(H2 if np.mod(W,2) else H1)+"' direction='"+str(detx_1[i])+','+str(detz_1[i])+",0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"
                if i<6:
                    string+="\t\t<TubeDetector1D position='"+str(detx_2[i])+','+str(detz_2[i])+','+str(H1 if np.mod(W,2) else H2)+"' direction='"+str(detx_2[i])+','+str(detz_2[i])+",0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"
            string+="\t</Wedge>\n"
            
        string+="</Instrument>"

        if fileName.split('.')[-1]!='xml':
            fileName+='.xml'

        with open(fileName,'w') as f:
            f.write(string)

    @_tools.KwargChecker()
    def generateCalibration(self,Vanadiumdatafile,A4datafile=False,savelocation='calibration/',tables=['Single','PrismaticLowDefinition','PrismaticHighDefinition'],plot=False,mask=True):
        """Method to generate look-up tables for normalization. Saves calibration file(s) as 'Calibration_Np.calib', where Np is the number of pixels.
        
        Generates 4 different tables:

            - Prismatic High Definition (8 pixels/energy or 64 pixels/detector)

            - Prismatic Low Definition (3 pixels/energy or 24 pixels/detector)

            - Single (1 pixel/energy or 8 pixels/detector)

            - Number (integer)
        
        Args:

            - Vanadiumdatafile (string): String to single data file used for normalization, Vanadium Ei scan (required).

        Kwargs:

            - A4datafile (string): String to single data file used for normalization, AlO A4 scan (default False).

            - savelocation (string): String to save location folder (calibration)

            - tables (list): List of needed conversion tables (Default: ['Single','PrismaticLowDefinition','PrismaticHighDefinition'], increasing number of pixels).

            - plot (boolean): Set to True if pictures of all fit are to be stored in savelocation

            - mask (boolean): If True the lower 100 pixels are set to 0

        .. warning::
            At the moment, the active detector area is defined by NumberOfSigmas (currently 3) times the Gaussian width of Vanadium peaks.

        """
        self.initialize()
        
        with hdf.File(Vanadiumdatafile,'r') as VanFile:
            if not A4datafile == False: # pragma: no cover
                A4File = hdf.File(A4datafile,'r')
                A4FileInstrument = getInstrument(A4File)
                A4FileInstrumentType = A4FileInstrument.name.split('/')[-1]


            VanFileInstrument = getInstrument(VanFile)
            

            VanFileInstrumentType = VanFileInstrument.name.split('/')[-1]
            
            if not A4datafile == False: # pragma: no cover
                if VanFileInstrumentType == A4FileInstrumentType:
                    InstrumentType = VanFileInstrumentType
                else:
                    raise AttributeError('The provided Vanadium and Powder files does not have the same instrument type ({} and {} respectively).'.format(VanFileInstrumentType,A4FileInstrumentType))        
            InstrumentType = VanFileInstrumentType
            if InstrumentType=='CAMEA':

                if savelocation[-1]!='/':
                    savelocation+='/'
                
                Data = np.array(VanFileInstrument.get('detector/counts')).transpose(2,0,1).astype(float)
                
                if False: #Mask pixels where spurion is present
                    Data[:,:,:100]=0

                Ei = np.array(VanFileInstrument.get('monochromator/energy')).astype(float)
                analysers = 8
                pixels = self.wedges[0].detectors[0].pixels
                detectors = len(self.A4[0])*len(self.A4)
                detectorsorInWedge = len(self.A4[0])
                wedges = len(self.A4)
                if pixels!=Data.shape[2]:
                    raise ValueError('The number of pixels ({}) in the data file does not match instrument description ({})!'.format(pixels,Data.shape[2]))

                bins = []
                for table in tables:
                    if isinstance(table,int):
                        bins.append(table)
                    else:
                        raise AttributeError("Provided table attribute ({}) not recognized an integer.".format(table))
                if len(bins)==0:
                    raise AttributeError("No binning has been chosen for normalization routine.")
                # Initial finding of peaks
                peakPos = np.ones((detectors,analysers),dtype=float)*(-1)
                peakVal = np.zeros_like(peakPos,dtype=float)
                peakWidth = np.ones_like(peakPos,dtype=float)
                peakBackg = np.zeros_like(peakPos,dtype=float)

                # Looking only at pixel direction (integration over E)
                ESummedData = Data.sum(axis=1)
                dataSubtracted = np.array(ESummedData.copy(),dtype=float)

                

                if plot: # pragma: no cover
                    plt.ioff()
                    plt.figure(figsize=(16,11))
                    if not os.path.exists(savelocation+'Raw'):
                        os.makedirs(savelocation+'Raw')
                    for i in range(detectors):
                        plt.clf()
                        plt.scatter(np.arange(pixels),np.sum(Data[:][i],axis=0),s=5)
                        plt.ylim(0,np.max(np.sum(Data[i],axis=0))*1.1)
                        plt.xlabel('Pixel')
                        plt.ylabel('Intensity [arg]')
                        plt.title('Vanadium normalization detector '+str(i))
                        plt.tight_layout()
                        plt.savefig(savelocation+'Raw/detector'+str(i)+'.png',format='png', dpi=150)
                
                for j in range(analysers):
                    peakPos[:,j],peakVal[:,j] = findPeak(dataSubtracted) # Find a peak in data
                    for i in range(detectors):
                        guess = [peakVal[i,j],float(peakPos[i,j]),20,np.min(ESummedData[i])]
                        try:
                            res = scipy.optimize.curve_fit(Gaussian,np.arange(ESummedData.shape[1]),dataSubtracted[i,:],p0=[guess])
                        except RuntimeError:
                            raise RuntimeError('Fitting did not converge at detector {} analyser {}'.format(i,j))
                        peakPos[i,j] = res[0][1]
                        peakVal[i,j] = res[0][0]
                        peakWidth[i,j]= res[0][2]
                        if peakPos[i,j]>ESummedData.shape[1]:
                            raise ValueError('Peak found at {} for analyser {} and detector {}'.format(peakPos[i,j],j,i))
                        # Generate peak as the one fitted and subtract it from signal
                        x=np.arange(pixels)
                        y = Gaussian(x,peakVal[i,j],peakPos[i,j],peakWidth[i,j],peakBackg[i,j])
                        peak = y>peakVal[i,j]*0.05
                        dataSubtracted[i,peak]= 0
                
                if plot: # pragma: no cover
                    x = np.arange(pixels)
                    for k in range(wedges):
                        plt.clf()
                        plt.suptitle('Fits')
                        for i in range(detectorsorInWedge):
                            y=np.zeros_like(x,dtype=float)
                            plt.subplot(4, 4, i+1)
                            plt.scatter(np.arange(pixels),ESummedData[i+13*k],s=4)
                            for j in range(analysers):
                                y += Gaussian(x,peakVal[i+13*k,j],peakPos[i+13*k,j],peakWidth[i+13*k,j],peakBackg[i+13*k,j])
                                plt.plot([peakPos[i+13*k,j],peakPos[i+13*k,j]],[0,np.max(ESummedData[i+13*k])*1.1])
                            plt.plot(x,y,'k')
                            plt.xlabel('Pixel')
                            plt.ylabel('Intensity [arg]')
                            plt.title('Detector {}'.format(i+13*k))
                            plt.ylim(0,np.max(ESummedData[i+13*k])*1.1)

                        plt.tight_layout()
                        plt.savefig(savelocation+r'/Raw/Fit_wedge_'+str(k)+'.png',format='png', dpi=150)
                        print('Saving: {}'.format(savelocation+r'/Raw/Fit_wedge_'+str(k)+'.png'))

                
                ## Sort the positions such that peak 1 is the furthermost left peak and assert diff(pos)>100
                sortedPeakPosArg = np.argsort(peakPos,axis=1)
                sortedPeakPos = np.sort(peakPos,axis=1)
                sortedPeakPos[np.logical_or(sortedPeakPos>pixels,sortedPeakPos<0)]=5*pixels # High number

                sortedPeakPosArg2 = np.argsort(sortedPeakPos,axis=1)
                sortedPeakPos.sort(axis=1)

                #differences = np.diff(sortedPeakPos,axis=1)
                #outliers = np.zeros_like(peakPos,dtype=bool)
                #outliers[:,:-1]=differences<pixels/100
                #sortedPeakPos[outliers]=5*pixels
                sortedPeakPosArg3 = np.argsort(sortedPeakPos,axis=1)
                argSort = np.array([sortedPeakPosArg[i,sortedPeakPosArg2[i,sortedPeakPosArg3[i,:]]] for i in range(detectors)])
                sortedPeakPos = np.sort(sortedPeakPos,axis=1)
                peaks=np.sum(sortedPeakPos<7*pixels,axis=1) # Number of peaks found

                if np.any(peaks!=analysers):
                    raise ValueError('Wrong number of peaks, {} found in detector(s): {}\nIn total error in {} detector(s).'.format(peaks[peaks!=analysers],np.arange(peaks.shape[0])[peaks!=analysers],np.sum(peaks!=analysers)))

                pixelpos  = np.array([peakPos[i,argSort[i]] for i in range(detectors)])
                widths    = np.array([peakWidth[i,argSort[i]] for i in range(detectors)])

                ## Define the active detector area
                sigmas = NumberOfSigmas # Active area is all pixels inside of pm 3 sigmas

                lowerPixel = pixelpos-sigmas*widths
                upperPixel = pixelpos+sigmas*widths

                split = (lowerPixel[:,1:]-upperPixel[:,:-1])/2+upperPixel[:,:-1]

                extendedSplit=np.zeros((split.shape[0],split.shape[1]+2))
                extendedSplit[:,1:-1] = split
                extendedSplit[:,-1]=np.ones((split.shape[0]))*pixels

                x=np.arange(pixels)
                activePixels = np.zeros((detectors,analysers,pixels),dtype=bool)
                for i in range(detectors):
                    if plot: # pragma: no cover
                        plt.clf()
                        plt.title('Detector {} Active pixels'.format(i))
                        plt.scatter(x,ESummedData[i],s=4,color='black')
                    for j in range(analysers):
                        activePixels[i,j] = np.logical_and(x>lowerPixel[i,j],x<upperPixel[i,j])
                        if plot: plt.scatter(x[np.logical_and(x>lowerPixel[i,j],x<upperPixel[i,j])], # pragma: no cover
                            ESummedData[i,np.logical_and(x>lowerPixel[i,j],x<upperPixel[i,j])],s=4,color='red')
                    if plot: # pragma: no cover
                        plt.ylim(0,np.max(ESummedData[i])*1.1)
                        plt.xlabel('Pixel')
                        plt.ylabel('Intensity [arg]')
                        plt.savefig(savelocation+'/Raw/Active_'+str(i)+'.png',format='png', dpi=150)

                Eguess = np.zeros_like(peakPos,dtype=int)
                for i in range(Eguess.shape[0]):
                    for j in range(analysers):
                        Eguess[i,j]=np.argmax(Data[i,:,int(pixelpos[i,j])])
                
                fitParameters = []
                activePixelRanges = []
                for detpixels in bins:
                    if detpixels*analysers*3>len(Ei):
                        warnings.warn('Fitting might be unstable due to {} pixels being fitted using only {} energies ({} free parameters).'.format(detpixels,len(Ei),detpixels*analysers*3),category=RuntimeWarning,stacklevel=2)
                        
                    if plot: # pragma: no cover
                        EiX = np.linspace(Ei[0],Ei[-1],len(Ei))
                        if not os.path.exists(savelocation+'/{}_pixels'.format(detpixels)):
                            os.makedirs(savelocation+'/{}_pixels'.format(detpixels)) 
                        colors=np.zeros((3,detpixels))
                        if pixels==1:
                            colors[:,0]=[0.65,0.2,0.45]
                        else:
                            colors[0]=np.linspace(0.3,1.0,detpixels)
                            colors[1]=np.linspace(0.2,0.2,detpixels)
                            colors[2]=np.linspace(0.8,0.1,detpixels)
                        plt.suptitle('{} pixels'.format(detpixels))

                    fittedParameters=np.zeros((detectors,analysers,detpixels,4))
                    activePixelDetector=[]
                    for i in range(detectors):
                        activePixelAnalyser = []
                        if plot: # pragma: no cover
                            plt.clf()
                            plt.title('Detector {}, {} pixels'.format(i,detpixels))
                            x =np.linspace(0,detpixels,len(Ei))
                        for j in range(analysers):
                            center = int(round(sortedPeakPos[i,j]))
                            width = activePixels[i,j].sum()
                            pixelAreas = np.linspace(-width/2.0,width/2.0,detpixels+1,dtype=int)+center+1 #Add 1 such that the first pixel is included 20/10-17
                            for k in range(detpixels):
                                binPixelData = Data[i,:,pixelAreas[k]:pixelAreas[k+1]].sum(axis=1)
                                ECenter = Ei[np.argmax(binPixelData)]
                                ECutLow = ECenter-0.4
                                ECutHigh= ECenter+0.4
                                TopId = np.argmin(np.abs(Ei-ECutHigh))
                                BotId = np.argmin(np.abs(ECutLow-Ei))
                                if TopId<BotId:
                                    _ = TopId
                                    TopId = BotId
                                    BotId = _
                                binPixelData = binPixelData[BotId:TopId]
                                EiLocal = Ei[BotId:TopId]
                                Bg = np.min(binPixelData[[0,-1]])
                                guess = np.array([np.max(binPixelData), ECenter,0.005,Bg],dtype=float)
                                try:
                                    res = scipy.optimize.curve_fit(Gaussian,EiLocal,binPixelData.astype(float),p0=guess)
                                    
                                except: # pragma: no cover
                                    if not os.path.exists(savelocation+'/{}_pixels'.format(detpixels)):
                                        os.makedirs(savelocation+'/{}_pixels'.format(detpixels))
                                    if not plot:
                                        plt.ioff
                                    plt.figure()
                                    plt.scatter(EiLocal,binPixelData)
                                    plt.plot(Ei,Gaussian(Ei,*guess))
                                
                                    plt.savefig(savelocation+'/{}_pixels/Detector{}_{}.png'.format(detpixels,i,k),format='png',dpi=150)
                                    plt.close()

                                fittedParameters[i,j,k]=res[0]
                                if plot: # pragma: no cover
                                    plt.plot(EiX,Gaussian(EiX,*fittedParameters[i,j,k]),color='black')
                                    plt.scatter(EiLocal,binPixelData,color=colors[:,k])
                            activePixelAnalyser.append(np.linspace(-width/2.0,width/2.0,detpixels+1,dtype=int)+center+1)
                        activePixelDetector.append(activePixelAnalyser)
                        if plot: # pragma: no cover
                            plt.grid('on')
                            plt.xlabel('Ei [meV]')
                            plt.ylabel('Weight [arb]')
                            plt.tight_layout(rect=(0,0,1,0.95))
                            plt.savefig(savelocation+'/{}_pixels/Detector{}.png'.format(detpixels,i),format='png',dpi=150)
                            print('Saving: {}'.format(savelocation+'/{}_pixels/Detector{}.png'.format(detpixels,i)))

                    if not A4datafile is False: # pragma: no cover
                        # Perform A4 calibration
                        A4FileValue = np.array(A4FileInstrument.get('detector/polar_angle'))
                        EiFile = np.array(A4FileInstrument.get('monochromator/energy'))[0]
                        A4FileIntensity = np.array(A4FileInstrument.get('detector/data'))

                        factorsqrtEK = 0.694692
                        ki = np.sqrt(EiFile)*factorsqrtEK

                        Qvec = 1.8049 # Angstrom <----------------------CHANGE!

                        # q = 2 k sin(theta)
                        theta = np.arcsin(Qvec/(2*ki))
        
                        A4 = np.array(self.A4)
                        A4=A4.reshape(A4.shape[0]*A4.shape[1],A4.shape[2],order='C')
                        EPrDetector = len(self.wedges[0].detectors[0].split)+1

                        
                        pixelList = np.array(activePixelDetector).reshape(A4.shape[0],EPrDetector,detpixels+1).astype(int)
                        PixelEdge = np.array([[pixelList[:,:,i],pixelList[:,:,i+1]] for i in range(detpixels)]).transpose((2,3,0,1))
                        PixelEnergy = fittedParameters[:,:,:,1].reshape(A4.shape[0],EPrDetector*detpixels)

                        ## Find detector analyser combi corresponding to energy
                        SoftwarePixel = np.array([np.argmin(np.abs(x-EiFile)) for x in PixelEnergy])

                        MeanA4Instr = np.zeros((A4.shape[0],EPrDetector*detpixels))
                        MeanIntensity = np.zeros((len(A4FileValue),A4.shape[0],EPrDetector*detpixels))
                        for i in range(A4.shape[0]): # For each detector
                            for j in range(EPrDetector):
                                for k in range(detpixels):
                                    MeanIntensity[:,i,j*detpixels+k] = np.sum(A4FileIntensity[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)
                                    MeanA4Instr[i,j*detpixels+k] = np.mean(A4[i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]])
                                    
                        x = A4FileValue
                        A4FitValue = np.zeros((A4.shape[0]))

                        
                        if plot==True:
                            plt.clf()
                        for i in range(104):
                            y = MeanIntensity[:,i,SoftwarePixel[i]]
                            if plot==True:
                                plt.scatter(x,y)
                            
                            guess=[np.max(y),x[np.argmax(y)],3,0]
                            try:
                                fit = scipy.optimize.curve_fit(Gaussian,x,y,p0=[guess])
                            except:
                                A4FitValue[i]=guess[1]
                            else:
                                A4FitValue[i] = fit[0][1]

                        if plot==True: # pragma: no cover
                            if not os.path.exists(savelocation+'A4'):
                                os.makedirs(savelocation+'A4')
                            plt.savefig(savelocation+'A4'+'/A4_{}.png'.format(detpixels),format='png',dpi=150)

                        A4FitValue+=2*theta*180.0/np.pi # offset relative to expected from powder line

                        if plot==True: # pragma: no cover
                            plt.clf()
                            plt.scatter(range(A4.shape[0]),A4FitValue)
                            plt.scatter(range(A4.shape[0]),MeanA4Instr[:,int(np.round(np.mean(SoftwarePixel)))]*180.0/np.pi)
                            plt.legend(['File','Geometry'])
                            plt.savefig(savelocation+'A4'+'/Points_{}.png'.format(detpixels),format='png',dpi=150)

                            diff = A4FitValue-MeanA4Instr[:,int(np.round(np.mean(SoftwarePixel)))]*180.0/np.pi#+2*theta*180.0/np.pi
                            plt.clf()
                            plt.scatter(range(A4.shape[0]),diff)
                            plt.savefig(savelocation+'A4'+'/diff_{}.png'.format(detpixels),format='png',dpi=150)

                    else: # Use nominal A4 values from calculation
                        #A4FitValue = []
                        #for i in range(8):
                        #    for j in range(13):
                        #        A4FitValue.append(-(i*8+j*0.55))
                        #A4FitValue = np.array(A4FitValue)
                        A4 = np.array(self.A4).reshape(104,1024)
                        A4Pixel = []
                        for i in range(len(fittedParameters)):
                            for j in range(len(fittedParameters[i])):
                                for k in range(len(fittedParameters[i][j])):
                                    #print(activePixelDetector[i][j][k],activePixelDetector[i][j][k+1])
                                    A4Pixel.append(np.mean(A4[i,activePixelDetector[i][j][k]:activePixelDetector[i][j][k+1]]))
                        A4Pixel = np.array(A4Pixel).reshape(len(fittedParameters),len(fittedParameters[i]),len(fittedParameters[i][j]))
                        #print(A4Pixel.shape)
                        #print(len(activePixelDetector))
                        #print(len(activePixelDetector[0]))
                        #print(len(activePixelDetector[0][0]))
                        #print(len(fittedParameters),len(fittedParameters[i]),len(fittedParameters[i][j]))
                        #h=kk
                        A4FitValue = np.rad2deg(A4Pixel)

                    fitParameters.append(fittedParameters)
                    activePixelRanges.append(np.array(activePixelDetector))
                    tableString = 'Normalization for {} pixel(s) using VanData {} and A4Data{}\nPerformed {}\nDetector,Energy,Pixel,Amplitude,Center,Width,Background,lowerBin,upperBin,A4Offset\n'.format(detpixels,Vanadiumdatafile,A4datafile,datetime.datetime.now())
                    for i in range(len(fittedParameters)):
                        for j in range(len(fittedParameters[i])):
                            for k in range(len(fittedParameters[i][j])):
                                tableString+=str(i)+','+str(j)+','+str(k)+','+','.join([str(x) for x in fittedParameters[i][j][k]])
                                tableString+=','+str(activePixelRanges[-1][i][j][k])+','+str(activePixelRanges[-1][i][j][k+1])
                                tableString+=','+str(A4FitValue[i,j,k])+'\n'
                    tableName = 'Normalization_{}.calib'.format(detpixels)
                    print('Saving {} pixel data to {}'.format(detpixels,savelocation+tableName))
                    file = open(savelocation+tableName,mode='w')

                    file.write(tableString)
                    file.close()
        
            if not A4datafile is False: # pragma: no cover
                A4File.close()

def parseXML(Instr,fileName):
    import xml.etree.ElementTree as ET

    tree = ET.parse(fileName)
    instr_root = tree.getroot()

    
    
        
    for attrib in instr_root.keys():
        if attrib=='position':
            Instr.position = np.array(instr_root.attrib[attrib].split(','),dtype=float)
        Instr.settings[attrib]=instr_root.attrib[attrib]
    
    
    Instr._wedges=[]
    for wedge in list(instr_root):#.getchildren():
        
        if wedge.tag in dir(Wedge):
            Wedgeclass_ = getattr(Wedge, wedge.tag)
        else:
            raise ValueError("Element is supposed to be a Wedge, but got '{}'.".format(wedge.tag))
        wedgeSettings = {}
        
        for attrib in wedge.keys():
            if attrib=='concept':
                wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=str)
            else:        
                wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=float)
            
        temp_wedge = Wedgeclass_(**wedgeSettings)
        
        
        
        
        for item in list(wedge):#.getchildren():
            if item.tag in dir(Detector):
                class_ = getattr(Detector, item.tag)
            elif item.tag in dir(Analyser):
                class_ = getattr(Analyser,item.tag)
            else:
                raise ValueError("Item '{}' not recognized as MJOLNIR detector or analyser.".format(item.tag))
            
            itemSettings = {}
            for attrib in item.keys():
                attribVal = item.get(attrib).strip().split(',')
                if len(attribVal)==1:
                    if(attrib=='split'):
                        try:
                            val=float(attribVal[0])
                        except ValueError:
                            val=[]
                        itemSettings[attrib]=val
                    else:
                        itemSettings[attrib]=float(attribVal[0])
                else:
                    if(attrib=='split'):
                        #print(type(attribVal))
                        itemSettings[attrib]=attribVal
                    else:
                        itemSettings[attrib]=np.array(attribVal,dtype=float)    
            try:
                temp_item = class_(**itemSettings)
            except TypeError as e:
                print(e.args[0])
                raise ValueError('Item {} misses argument(s):{}'.format(class_,e.args[0].split(':')[0]))
            except AttributeError as e:
                raise AttributeError('Error in passing {} with attributes {}'.format(class_,itemSettings))
            except ValueError:
                raise ValueError('Item {} not initialized due to error.'.format(class_))
            #print(temp_item)
            temp_wedge.append(temp_item)
            #print()

        #print(str(temp_wedge))
        Instr.append(temp_wedge)
    #print(str(Instr))
    
   

def getNX_class(x,y,attribute):
    try:
        variableType = y.attrs['NX_class']
    except:
        variableType = ''
    if variableType==attribute:
        return x

def getInstrument(file):
    location = file.visititems(lambda x,y: getNX_class(x,y,b'NXinstrument'))
    return file.get(location)



def Gaussian(x,A,mu,sigma,b):
    return A*np.exp(-np.power(mu-x,2.0)*0.5*np.power(sigma,-2.0))+b


def findPeak(data):
    return [np.argmax(data,axis=1),np.max(data,axis=1)]

def convertToHDF(fileName,title,sample,fname,CalibrationFile=None,pixels=1024,cell=[5,5,5,90,90,90],factor=10000): # pragma: no cover
    """Convert McStas simulation to h5 format.
    
    Args:

        - fileName (str): File name of created file ('*.hdf')

        - title (str): Title of HdF file

        - sample (str): Name of sample

        - fname (str): Location folder of McStas Data (must end with '/')

    Kwargs:

        - CalibrationFile (str or list of str): Location of calibration file(s) wanted in HdF file (default None)

        - pixels (int): Number of pixels on detectors (default 1024)

        - cell (list): Cell parameters passed into the hdf file (default [5,5,5,90,90,90])

        - factor (float): Overall scale factor for intensity

    """
    def addMetaData(entry,title):
        dset = entry.create_dataset('start_time',(1,),dtype='<S70')
        dset[0] = b'2018-03-22T16:44:02+01:00'

        dset = entry.create_dataset('end_time',(1,),dtype='<S70')
        dset[0] = b"2018-03-22T18:44:02+01:00"

        dset = entry.create_dataset('experimental_identifier',(1,),dtype='<S70')
        dset[0] = b"UNKNOWN"

        dset = entry.create_dataset('instrument',(1,),dtype='<S70')
        dset[0] = b"CAMEA"

        dset = entry.create_dataset('comment',(1,),dtype='<S70')
        dset[0] = b"I feel uncommented"

        dset = entry.create_dataset('title',(1,),dtype='<S70')
        dset[0] = np.string_(title)

        dset = entry.create_dataset('proposal_id',(1,),dtype='<S70')
        dset[0] = b"2018-00666"

        dset = entry.create_dataset('proposal_title',(1,),dtype='<S70')
        dset[0] = b"I need my Title!"

        cont = entry.create_group('local_contact')
        cont.attrs['NX_class'] = np.string_('NXuser')
        dset = cont.create_dataset('name',(1,),dtype='S70')
        dset[0] = b"UNKNOWN"

        us = entry.create_group('proposal_user')
        us.attrs['NX_class'] = np.string_('NXuser')
        dset = us.create_dataset('name',(1,),dtype='S70')
        dset[0] = b"Jakob Lass"
        dset = us.create_dataset('email',(1,),dtype='S70')
        dset[0] = b"jakob@lass.dk"

        pus = entry.create_group('user')
        pus.attrs['NX_class'] = np.string_('NXuser')
        dset = pus.create_dataset('name',(1,),dtype='S70')
        dset[0] = b"Jakob Lass"
        dset = pus.create_dataset('email',(1,),dtype='S70')
        dset[0] = b"jakob@lass.dk"
        dset = pus.create_dataset('address',(1,),dtype='S70')
        dset[0] = b"UNKNOWN"
        dset = pus.create_dataset('affiliation',(1,),dtype='S70')
        dset[0] = b"UNKNOWN"

        

    def addMono(inst):
        mono = inst.create_group('monochromator')
        mono.attrs['NX_class'] = np.string_('NXmonochromator')

        dset = mono.create_dataset('type',(1,),dtype='S70')
        dset[0] = b"Pyrolithic Graphite"
        
        dset = mono.create_dataset('d_spacing',(1,),'float32')
        dset[0] = 3.354
        dset.attrs['units'] = 'anstrom'

        dset = mono.create_dataset('horizontal_curvature',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'meter'

        dset = mono.create_dataset('vertical_curvatur',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'meter'

        dset = mono.create_dataset('horizontal_curvature_zero',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'meter'

        dset = mono.create_dataset('vertical_curvature_zero',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'meter'

        dset = mono.create_dataset('gm',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'degree'

        dset = mono.create_dataset('gm_zero',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'degree'

        dset = mono.create_dataset('tlm',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'degree'

        dset = mono.create_dataset('tlm_zero',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'degree'

        dset = mono.create_dataset('tum',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'degree'

        dset = mono.create_dataset('tum_zero',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'degree'

        monoSlit = inst.create_group('monochromator_slit')
        monoSlit.attrs['NX_class'] = np.string_('NXmonochromatorslit')

        for x in ['bottom','left','right','top']:
            dset = monoSlit.create_dataset(x,(1,),'float32')
            dset[0] = 0.0
            dset.attrs['units'] = 'mm'

            dset = monoSlit.create_dataset(x+'_zero',(1,),'float32')
            dset[0] = 0.0
            dset.attrs['units'] = 'mm'
        for x in ['x_gab','y_gab']:
            dset = monoSlit.create_dataset(x,(1,),'float32')
            dset[0] = 0.0
            dset.attrs['units'] = 'mm'
    
    def addAna(inst):
        ana = inst.create_group('analyzer')
        ana.attrs['NX_class'] = np.string_('NXcrystal')

        dset = ana.create_dataset('type',(1,),dtype='S70')
        dset[0] = b"Pyrolithic Graphite"
        
        dset = ana.create_dataset('d_spacing',(1,),'float32')
        dset[0] = 3.354
        dset.attrs['units'] = 'anstrom'

        dset = ana.create_dataset('nominal_energy',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'mev'

        

    def addDetector(inst):
        det = inst.create_group('detector')
        det.attrs['NX_class'] = np.string_('NXdetector')

    def readDetSequence():
        detlist = []
        dir_path = os.path.dirname(os.path.realpath(__file__))
        fin = open(dir_path+'/'+'detsequence.dat','r')
        for line in fin:
            detlist.append(line.strip())
        fin.close()
        return detlist

    def readDetFile(fname,pixels=1024,factor=10000):
        detdata = np.zeros((pixels),dtype='int32')
        f = open(fname,'r')
        psddata = f.readlines()
        f.close()
        idx = 0
        for line in psddata:
            if line.find('EI=') > 0:
                l = line.split('=')
                ei = float(l[1])
            if line.find('A3=') > 0:
                l = line.split('=')
                a3 = float(l[1])
            if line.find('A4=') > 0:
                l = line.split('=')
                a4 = float(l[1])
            if line.find('variables:') > 0:
                idx = idx + 1
                break
            idx = idx + 1
        detind = 0
        for i in range(idx+1,pixels+idx-1):
            l = psddata[i].split()

            detdata[detind] = int(round(factor*float(l[1])))
            #if l[1]!='0':
            #    print(float(l[1])*factor)
            detind = detind + 1
        return detdata,a3,a4,ei

    def readScanPointData(dir,detlist,Numpoints,pixels=1024,factor=10000):
        frame = np.zeros((104,pixels),dtype='int32')
        i = 0
        for detfile in detlist:
            detdata, a3, a4, ei = readDetFile(dir +'/' + str(Numpoints) + '/' + detfile,factor=factor)
            frame[i] = detdata
            i = i + 1
        return frame,a3,a4,ei

    def readScanData(dir,Numpoints,pixels=1024,factor=10000):
        detlist = readDetSequence()
        data = np.zeros((Numpoints,104,pixels),dtype='int32')
        a3 = []
        a4 = []
        ei = []
        for n in range(Numpoints):
            frame, a3n, a4n, ein = readScanPointData(dir,detlist,n,factor=factor)
            a3.append(a3n)
            a4.append(a4n)
            ei.append(ein)
            data[n] = frame
        return data,a3,a4,ei
        
    def addSample(entry,name,cell):
        sam = entry.create_group('sample')
        sam.attrs['NX_class'] = np.string_('NXsample')
        dset = sam.create_dataset('name',(1,),dtype='S70')
        dset[0] = np.string_(name)

        ub = np.zeros((3,3,),dtype='float32')
        ub[0,0] = 1.
        ub[1,1] = 1.
        ub[2,2] = 1.
        dset = sam.create_dataset('orientation_matrix',data=ub)
        dset = sam.create_dataset('plane_vector_1',data=[1,0,0,0,0,0,0])
        dset = sam.create_dataset('plane_vector_2',data=[0,1,0,0,0,0,0])

        normal = np.zeros((3,),dtype='float32')
        normal[2] = 1.0
        dset = sam.create_dataset('plane_normal',data=normal)

        cell = np.array(cell,dtype='float32')
        dset = sam.create_dataset('unit_cell',data=cell)

        dset = sam.create_dataset('azimuthal_angle',data=0.0)
        dset = sam.create_dataset('x',data=0.0)
        dset = sam.create_dataset('y',data=0.0)

        for x in ['sgu','sgl']:
            dset = sam.create_dataset(x,data=0.0)
            dset = sam.create_dataset(x+'_zero',data=0.0)


        

    def isVaried(data):
        if len(data)>1 and data[0]!=data[1]:
            return True
        else:
            return False

    def makeTheta(ei):
        theta = []
        tth = []
        for e in ei:
            k = np.sqrt(float(e)/2.072)
            fd = np.pi/(k*3.354)
            th = np.degrees(np.arcsin(fd))
            theta.append(th)
            tth.append(2.*th)
        return theta,tth

    def storeScanData(entry,data,a3,a4,ei):
        nxdata = entry.create_group('data')
        nxdata.attrs['NX_class'] = np.string_('NXdata')
        
        det = entry['CAMEA/detector']
        dset = det.create_dataset('counts',data=data.swapaxes(1,2), compression="gzip", compression_opts=9)
        dset.attrs['target'] = np.string_('/entry/CAMEA/detector/counts')
        nxdata['counts'] = dset

        dset = det.create_dataset('summed_counts',data=np.sum(data,axis=(1,2)))
        dset.attrs['target'] = np.string_('/entry/CAMEA/detector/summed_counts')
        nxdata['summed_counts'] = dset

        sam = entry['sample']

        scanType='Unknown'
        scanvars = ''
        if isVaried(a3):
            dset = sam.create_dataset('rotation_angle',data=a3)
            dset_zero = sam.create_dataset('rotation_angle_zero',data=np.array([0.0]))
            nxdata['rotation_angle'] = dset
            nxdata['rotation_angle_zero'] = dset_zero
            scanType = 'cscan a3 {} da3 {} np {} mn 10000'.format(np.mean(a3),np.mean(np.diff(a3)),len(a3))
            scanvars+='a3'
        else:
            dset = sam.create_dataset('rotation_angle',(1,),dtype='float32')
            dset_zero = sam.create_dataset('rotation_angle_zero',data=np.array([0.0]))

        dset.attrs['units'] = np.string_('degrees')
        dset_zero.attrs['units'] = np.string_('degrees')

        if isVaried(a4):
            dset = sam.create_dataset('polar_angle',data=a4)
            nxdata['polar_angle'] = dset
            dset_zero = sam.create_dataset('polar_angle_zero',(1,),dtype='float32',data=0.0)
            nxdata['polar_angle_zero'] = dset_zero
            scanType = 'cscan a4 {} da4 {} np {} mn 10000'.format(np.mean(a4),np.mean(np.diff(a4)),len(a4))
            scanvars+='a4'
        else:
            dset = sam.create_dataset('polar_angle',(1,),dtype='float32',data=a4[0])
            dset_zero = sam.create_dataset('polar_angle_zero',(1,),dtype='float32',data=0.0)
        dset.attrs['units'] = np.string_('degrees')
        dset_zero.attrs['units'] = np.string_('degrees')
        

        mono = entry['CAMEA/monochromator']
        theta,tth = makeTheta(ei)

        if isVaried(ei):
            dset = mono.create_dataset('energy',data=ei)
            nxdata['incident_energy'] = dset
            mono.create_dataset('rotation_angle',data=theta);
            mono.create_dataset('polar_angle',data=tth)
            scanType = 'cscan ei {} dei {} np {} mn 10000'.format(np.mean(ei),np.mean(np.diff(ei)),len(ei))
            scanvars+='ei'
        else:
            dset = mono.create_dataset('energy',(1,),dtype='float32')
            dset[0] = ei[0]
            dset = mono.create_dataset('rotation_angle',(1,),dtype='float32')
            dset[0] = theta[0]
            dset = mono.create_dataset('polar_angle',(1,),dtype='float32')
            dset[0] = tth[0]
        dset = entry['CAMEA/monochromator/rotation_angle']    
        dset.attrs['units'] = np.string_('degrees')
        dset = entry['sample/polar_angle']    
        dset.attrs['units'] = np.string_('degrees')
        

        dset.attrs['units'] = np.string_('degrees')

        #dset = mono.create_dataset('summed_counts',data=np.sum(data,axis=(1,2)));
        #dset.attrs['units'] = np.string_('counts')
        
        
        makeMonitor(entry,Numpoints)
        entry.create_dataset('scancommand',data=scanType)
        entry.create_dataset('scanvars',data=scanvars)
    def makeMonitor(entry,Numpoints):
        control = entry.create_group('control')
        control.attrs['NX_class'] = np.string_('NXmonitor')
        mons = [10000]*Numpoints
        control.create_dataset('data',data=mons,dtype='int32')
        dset = control.create_dataset('preset',(1,),dtype='int32')
        dset[0] = 10000
        dset = control.create_dataset('mode',(1,),dtype='S70')
        dset[0] = b"monitor"
        time = [36.87]*Numpoints
        control.create_dataset('time',data=time,dtype='float32')

        time = [36.87*1e9]*Numpoints
        control.create_dataset('absolut_time',data=time,dtype='float32')
        
        pb = entry.create_group('proton_beam')
        pb.attrs['NX_class'] = np.string_('NXmonitor')
        vals = [0]*Numpoints
        dset = pb.create_dataset('data',data=vals,dtype='int32')

    f = hdf.File(fileName,'w')
    f.attrs['file_name'] = np.string_(fileName)
    f.attrs['file_time'] = np.string_(b'2018-03-22T16:44:02+01:00')
    
    entry = f.create_group('entry')
    entry.attrs['NX_class'] = np.string_('NXentry')

    addMetaData(entry,np.string_(title))
    
    #------------ Instrument
    inst = entry.create_group(b'CAMEA')
    inst.attrs['NX_class'] = np.string_('NXinstrument')
    
    if not CalibrationFile is None:
        #calib = inst.create_group(b'calibration')
        if not isinstance(CalibrationFile,list):
            CalibrationFile=[CalibrationFile]
        for i in range(len(CalibrationFile)):
            calibrationData = np.genfromtxt(CalibrationFile[i],skip_header=3,delimiter=',')
            binning = CalibrationFile[i].split('/')[-1].split('_')[-1].split('.')[0]
            pixelCalib = inst.create_group('calib{}'.format(binning))

            pixelCalib.create_dataset('final_energy'.format(binning),data=calibrationData[:,4],dtype='float32')
            pixelCalib.create_dataset('background'.format(binning),data=calibrationData[:,6],dtype='float32')
            pixelCalib.create_dataset('width'.format(binning),data=calibrationData[:,5],dtype='float32')
            pixelCalib.create_dataset('amplitude'.format(binning),data=calibrationData[:,3],dtype='float32')
            pixelCalib.create_dataset('boundaries'.format(binning),data=calibrationData[:,7:9],dtype='int')
            pixelCalib.create_dataset('a4offset'.format(binning),data=calibrationData[:,9],dtype='float32')
            
    
    addMono(inst)
    addAna(inst)
    
    addDetector(inst)
    
    addSample(entry,np.string_(sample),cell)
    import os
    Numpoints = sum([os.path.isdir(fname+'/'+i) for i in os.listdir(fname)])
    data,a3,a4,ei = readScanData(fname,Numpoints,factor=factor)
    storeScanData(entry,data,a3,a4,ei)
    f.close()


# _________________________TESTS____________________________________________



def test_Instrument_init():
    Instr = Instrument()

    assert(np.all(Instr.position==(0,0,0)))

    Det = Detector.Detector(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.Analyser(position=(0.5,0,0),direction=(1,0,1))
    
    wedge = Wedge.Wedge(detectors=[Det,Det],analysers=Ana)

    Instr.wedges=[wedge,wedge]

    assert(Instr.settings['Initialized']==False)



def test_Instrument_error():
    
    try:
        Instr = Instrument(fileName='wrongDummyFile.bin')
        assert False
    except ValueError:
        assert True

    Instr = Instrument()

    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    try:
        Instr.wedges=Ana
        assert False
    except AttributeError:
        assert True

    try:
        Instr.wedges=[Ana,Ana]
        assert False
    except AttributeError:
        assert True

    try:
        Instr.append("Wrong object type")
        assert False
    except AttributeError:
        assert True
    
    try:
        Instr.append(["List of",3.0,"wrong objects"])
        assert False
    except AttributeError:
        assert True

    try:
        Instr.settings = {'Name','New dictionary'}
        assert False
    except NotImplementedError:
        return True


def test_Instrument_warnings():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.wedges = wedge

    with warnings.catch_warnings(record=True) as w: # From https://docs.python.org/3.1/library/warnings.html
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        Instr.wedges = wedge
        # Verify some things
        assert len(w) == 1
        assert issubclass(w[0].category, UserWarning)
        assert 'The list of wedges is not empty! Appending new wedges(s)' in str(w[0].message)


def test_Instrument_append():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Instr.append([wedge,wedge])
    Instr.append(wedge)

    assert(len(Instr.wedges)==3)


def test_Instrument_plot():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append(wedge)
    plt.ioff()
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    Instr.plot(ax)

def test_Instrument_Setting(): 
    Instr = Instrument()
    Instr.settings['SettingVersion']=1.0
    assert(Instr.settings['SettingVersion']==1.0)


def test_Instrument_Initialization():
    Instr = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0),concept='ManyToMany')
    pixels=33
    split = [12]
    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0),pixels=pixels,split=split)
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    

    wedge.append([Det,Det,Ana,Ana,Ana])

    try:
        Instr.initialize()
        assert False
    except ValueError:
        assert True

    try:
        print(Instr.A4)
        assert False
    except RuntimeError:
        assert True

    try:
        print(Instr.Ef)
        assert False
    except RuntimeError:
        assert True

    Instr.append(wedge)
    try:
        Instr.initialize()
        assert False
    except ValueError:
        assert True
    Instr.wedges[0].detectors[0].split = [0,12,20,pixels]
    Instr.initialize()

    assert(len(Instr.A4)==1)
    assert(len(Instr.A4[0])==2)
    assert(len(Instr.A4[0][0])==pixels)
    assert(len(Instr.A4)==len(Instr.Ef))
    assert(len(Instr.A4[0])==len(Instr.Ef[0]))
    assert(len(Instr.A4[0][0])==len(Instr.Ef[0][0]))
    assert(Instr.settings['Initialized']==True)

    try:
        Instr.A4 = []
        assert False
    except NotImplementedError:
        assert True

    try:
        Instr.Ef = []
        assert False
    except NotImplementedError:
        assert True



def test_Instrument_saveload():
    import os
    Instr = Instrument(position=(0,1,0))
    Instr2 = Instrument()

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append(wedge)

    tempFile = 'temp.bin'
    Instr.save(tempFile)
    Instr2.load(tempFile)
    os.remove(tempFile)
    

    assert(Instr==Instr2)



def test_parseXML(): # Improve this test!

    tempFileName = '__temp__.xml'
        
    Instr = Instrument()
    Instr.settings['Author'] = 'Jakob Lass'

    wedge = Wedge.Wedge(position=(0.5,0,0))

    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

    wedge.append([Det,Ana])
    Instr.append([wedge,wedge])
    Instr.append(wedge)
    Instr.saveXML(tempFileName)
        
    InstrLoaded = Instrument(fileName=tempFileName)
    os.remove(tempFileName)

    assert(Instr==InstrLoaded) 


def test_XML_errors():

    fileString = ""
    fileString+="<?xml version='1.0'?>"
    fileString+="<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>"
    fileString+="<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>"
    fileString+="<FlatAnalyser direction='0.707,0.0,0.707' d_spacing='3.35' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>"
    fileString+="<TubeDetector1D position='1.198,0.0580,0.71' direction='0.998,0.04841,0.0' pixels='456' length='0.883' diameter='0.02' split='57, 114, 171, 228, 285, 342, 399'></TubeDetector1D>"
    fileString+="</Wedge>"
    fileString+="</Instrument>"

    temp_file = 'Tempfile.xml'
    f = open(temp_file,'w')
    f.write(fileString)
    f.close()

    try:
        Instr = Instrument(fileName=temp_file)
        del Instr
        assert False
    except ValueError:
        assert True

    fileString = ""
    fileString+="<?xml version='1.0'?>"
    fileString+="<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>"
    fileString+="<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>"
    fileString+="<FlatAnalyser position='0.0580,0.71' direction='0.707,0.0,0.707' d_spacing='3.35' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>"
    fileString+="<TubeDetector1D position='1.198,0.0580,0.71' direction='0.998,0.04841,0.0' pixels='456' length='0.883' diameter='0.02' split='57, 114, 171, 228, 285, 342, 399'></TubeDetector1D>"
    fileString+="</Wedge>"
    fileString+="</Instrument>"
    f = open(temp_file,'w')
    f.write(fileString)
    f.close()
    try:
        Instr = Instrument(fileName=temp_file)
        assert False
    except AttributeError:
        assert True

    fileString = ""
    fileString+="<?xml version='1.0'?>"
    fileString+="<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>"
    fileString+="<FlatAnalyser position='0.0,0.0,0.0' concept='ManyToMany'>"
    fileString+="<FlatAnalyser position='0.0580,0.71' direction='0.707,0.0,0.707' d_spacing='3.35' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>"
    fileString+="<TubeDetector1D position='1.198,0.0580,0.71' direction='0.998,0.04841,0.0' pixels='456' length='0.883' diameter='0.02' split='57, 114, 171, 228, 285, 342, 399'></TubeDetector1D>"
    fileString+="</FlatAnalyser>"
    fileString+="</Instrument>"
    f = open(temp_file,'w')
    f.write(fileString)
    f.close()
    try:
        Instr = Instrument(fileName=temp_file)
        assert False
    except ValueError:
        assert True
    os.remove(temp_file)

def test_instrument_string_dummy(): # Todo: Improve test!
    Instr = Instrument()

    string = str(Instr)
    del string
    assert True
    
def test_instrument_create_xml():

    Instr = Instrument()
    filename = 'temp'
    Instr.generateCAMEAXML(filename)

    Instr2 = Instrument(fileName=filename+'.xml')
    os.remove(filename+'.xml')
    assert(len(Instr2.wedges)==8)


@pytest.mark.unit
def test_Normalization_tables(quick):

    Instr = Instrument(fileName='Data/CAMEA_Updated.xml')
    Instr.initialize()

    NF = 'Data/camea2018n000038.hdf'
    #AF = 'TestData/1024/A4Normalization.h5'

    try:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation='Data/',plot=False,tables=[]) # No binning specified 
        assert False
    except AttributeError:
        assert True

    try:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation='Data/',plot=False,tables=['Nothing?']) # Wrong binning
        assert False
    except AttributeError:
        assert True

    if not quick==True:
        Instr.generateCalibration(Vanadiumdatafile=NF,  savelocation='Data/',plot=False,tables=[1,3,8]) 
    else:
        Instr.generateCalibration(Vanadiumdatafile=NF ,savelocation='Data/',plot=False,tables=[1]) 


