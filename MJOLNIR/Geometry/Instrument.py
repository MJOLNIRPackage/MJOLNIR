import sys,os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import numpy as np
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector,Wedge
from MJOLNIR import _tools
from MJOLNIR import TasUBlibDEG as TasUBlib
from MJOLNIR.TasUBlibDEG import cosd,sind
from MJOLNIR.Data import Sample,RLUAxes
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize
from scipy.stats import norm
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
        self._pixelPosition = []
        self._analyserPosition = []
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
            self._pixelPosition.append(detectorPixelPositions)
            self._analyserPosition.append(analyserPixelPositions)




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

        offset = 28.0# -4.835960288880082# offset such that last pixel of detector 0 is at 0


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
    def generateCalibration(self,Vanadiumdatafile,A4datafile=False,savelocation='calibration/', 
    tables=['Single','PrismaticLowDefinition','PrismaticHighDefinition'], plot=False, mask=True, adaptiveBinning=False, ignoreTubes=None,
    sample='V',sampleMass=None,sampleDebyeWallerFactor=1.0,formulaUnitsPerUnitCell=1.0,sampleIncoherent=5.08):
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

            - adaptiveBinning (boolean): If true pixel bins are asigned to give same Gaussian area of intensity for all pixels (default False)

            - ignoreTubes (list of ints): List containing tubes to be ignored in fitting (default [])

            - sample (string): Chemical composition of sample used for normalization (default 'V')

            - sampleMass (float): Mass of normalization sample (default None)

            - sampleDebyeWallerFactor (float): Correction of normalization data with Debye-Waller factor (default 1.0)

            - formulaUnitsPerUnitCell (float): Numer of formalunits in normalization sample per unit cell (default 1.0) 
            
            - sampleIncoherent (float): Incoherent strength of sample (default 5.08)

        .. warning::
            At the moment, the active detector area is defined by NumberOfSigmas (currently 3) times the Gaussian width of Vanadium peaks.

        When a sample mass is provided, the normalization tables are per default rescaled with Vanadium unless default parameters are overwritten.

        """
        self.initialize()

        if ignoreTubes is None:
            ignoreTubes = []
        
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
                monitor = np.mean(np.array(VanFile.get('entry/control/data'))) # Extract monitor count
                
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
                        if i in ignoreTubes:
                            peakPos[i,j] = -1
                            peakVal[i,j] = -1
                            peakWidth[i,j]= 0
                            continue

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
                
                if sampleMass != None: # A sample mass has been proivided

                    sampleMolarMass = _tools.calculateMolarMass(sample,formulaUnitsPerUnitCell)
                    
                    # TODO: Check placement of sampleDebyeWallerFactor!
                    vanadiumFactor = sampleDebyeWallerFactor*sampleMolarMass/(sampleMass*sampleIncoherent*monitor)
                else:
                    vanadiumFactor = 1.0

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
                        if i in ignoreTubes:
                            for j in range(analysers):
                                for k in range(detpixels):
                                    fittedParameters[i,j,k,0]
                                activePixelAnalyser.append(np.linspace(-0/2.0,0/2.0,detpixels+1,dtype=int)+0+1)
                        else:
                            if plot: # pragma: no cover
                                plt.clf()
                                plt.title('Detector {}, {} pixels'.format(i,detpixels))
                                x =np.linspace(0,detpixels,len(Ei))
                            for j in range(analysers):
                                center = int(round(sortedPeakPos[i,j]))
                                width = activePixels[i,j].sum()
                                
                                if adaptiveBinning: # Adaptive binning with equal Gaussian area for each pixel piece

                                    binStart,binEnd = [center-width/2.0,center+width/2.0]
                                    binMid = center#(binStart+binEnd)/2.0

                                    binWidth = (binEnd-binStart)/(2.0*NumberOfSigmas)
                                    totalArea = norm.cdf(NumberOfSigmas)-norm.cdf(-NumberOfSigmas)
                                    areaOutside = norm.cdf(-NumberOfSigmas)

                                    areaPerBin = totalArea/detpixels

                                    areaLeftOfBin = [areaOutside+n*areaPerBin for n in range(detpixels+1)]

                                    pixelAreas = np.round(norm.ppf(areaLeftOfBin)*binWidth+binMid).astype(int)
                                else:
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
                                    fittedParameters[i,j,k,0] *= np.sqrt(2*np.pi)*fittedParameters[i,j,k,2] #np.sum(binPixelData) # Use integrated intensity as amplitude 29-07-2020 JL
                                    fittedParameters[i,j,k,0] *= vanadiumFactor # Multiply vanadium factor to perform absolut normalization 22-02-2021 JL
                                    if plot: # pragma: no cover
                                        plt.plot(EiX,GaussianNormalized(EiX,*fittedParameters[i,j,k]),color='black')
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
                    tableString = 'Normalization for {} pixel(s) using VanData {} and A4Data{}\nPerformed {}\nDetector,Energy,Pixel,IntegratedIntensity,Center,Width,Background,lowerBin,upperBin,A4Offset\n'.format(detpixels,Vanadiumdatafile,A4datafile,datetime.datetime.now())
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

def GaussianNormalized(x,A,mu,sigma,b):
    return A/(np.sqrt(2*np.pi)*sigma)*np.exp(-np.power(mu-x,2.0)*0.5*np.power(sigma,-2.0))+b


def findPeak(data):
    return [np.argmax(data,axis=1),np.max(data,axis=1)]

# TODO: Make the reader create a true HDF file with all attributes set
def convertToHDF(fileName,title,sample,fname,CalibrationFile=None,pixels=1024,cell=[5,5,5,90,90,90],factor=10000,detectors=104,\
    ub=None,plane_vector_1=None,plane_vector_2=None,plane_normal=None,rotation_angle_zero=0.0,polar_angle_offset=0.0): # pragma: no cover
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

        dset = entry.create_dataset('experiment_identifier',(1,),dtype='<S70')
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

        dset = ana.create_dataset('analyzer_selection',(1,),'float32')
        dset[0] = 0
        dset.attrs['units'] = 'anstrom'

        dset = ana.create_dataset('nominal_energy',(1,),'float32')
        dset[0] = 0.0
        dset.attrs['units'] = 'mev'

        

    def addDetector(inst):
        det = inst.create_group('detector')
        det.attrs['NX_class'] = np.string_('NXdetector')

        detsel = det.create_dataset('detector_selection',(1,),'float32')
        detsel[0] = 0

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
        a3 = 0
        a4 = 0
        ei = 0
        for line in psddata:
            if line.find('EI=') >0 or line.find(' SourceE=') > 0:
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

    def readScanPointData(dir,detlist,Numpoints,pixels=1024,factor=10000,detectors=104):
        frame = np.zeros((detectors,pixels),dtype='int32')
        i = 0
        for detfile in detlist:
            detdata, a3, a4, ei = readDetFile(dir +'/' + str(Numpoints) + '/' + detfile,factor=factor,pixels=pixels)
            frame[i] = detdata
            i = i + 1
        return frame,a3,a4,ei

    def readScanData(dir,Numpoints,pixels=1024,factor=10000,detectors=104):
        detlist = readDetSequence()[:detectors]
        data = np.zeros((Numpoints,detectors,pixels),dtype='int32')
        a3 = []
        a4 = []
        ei = []
        for n in range(Numpoints):
            frame, a3n, a4n, ein = readScanPointData(dir,detlist,n,factor=factor,pixels=pixels,detectors=detectors)
            a3.append(a3n)
            a4.append(a4n)
            ei.append(ein)
            data[n] = frame
        return data,a3,a4,ei
        
    def addSample(entry,name,cell,ub=None,plane_vector_1=None,plane_vector_2=None,plane_normal=None):
        sam = entry.create_group('sample')
        sam.attrs['NX_class'] = np.string_('NXsample')
        dset = sam.create_dataset('name',(1,),dtype='S70')
        dset[0] = np.string_(name)
        if ub is None:
            ub = np.zeros((3,3,),dtype='float32')
            ub[0,0] = 1.
            ub[1,1] = 1.
            ub[2,2] = 1.
        dset = sam.create_dataset('orientation_matrix',data=ub)
        if plane_vector_1 is None:
            plane_vector_1=[1,0,0,0,0,0,0]
        dset = sam.create_dataset('plane_vector_1',data=plane_vector_1)
        if plane_vector_2 is None:
            plane_vector_2=[0,1,0,0,0,0,0]
        dset = sam.create_dataset('plane_vector_2',data=plane_vector_2)

        if plane_normal is None:
            normal = np.zeros((3,),dtype='float32')
            normal[2] = 1.0

        dset = sam.create_dataset('plane_normal',data=plane_normal)

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

    def storeScanData(entry,data,a3,a4,ei,rotation_angle_zero=0.0,polar_angle_offset=0.0):
        nxdata = entry.create_group('data')
        nxdata.attrs['NX_class'] = np.string_('NXdata')
        
        det = entry['CAMEA/detector']
        dset = det.create_dataset('counts',data=data.swapaxes(1,2), compression="gzip", compression_opts=9)
        dset.attrs['target'] = np.string_('/entry/CAMEA/detector/counts')
        dset.attrs['units'] = np.string_('counts')
        nxdata['counts'] = dset

        dset = det.create_dataset('summed_counts',data=np.sum(data,axis=(1,2)))
        dset.attrs['target'] = np.string_('/entry/CAMEA/detector/summed_counts')
        dset.attrs['units'] = np.string_('counts')
        nxdata['summed_counts'] = dset

        sam = entry['sample']

        scanType='Unknown'
        scanvars = ''
        if isVaried(a3):
            dset = sam.create_dataset('rotation_angle',data=np.array(a3))
            dset_zero = sam.create_dataset('rotation_angle_zero',data=np.array([rotation_angle_zero]))
            dset.attrs['target'] = np.string_('/entry/sample/rotation_angle')
            nxdata['rotation_angle'] = dset
            scanType = 'cscan a3 {} da3 {} np {} mn 10000'.format(np.mean(a3),np.mean(np.diff(a3)),len(a3))
            scanvars+='a3'
        else:
            dset = sam.create_dataset('rotation_angle',(1,),dtype='float32',data=a3[0])
            dset_zero = sam.create_dataset('rotation_angle_zero',data=np.array([rotation_angle_zero]))

        dset.attrs['units'] = np.string_('degrees')
        dset_zero.attrs['units'] = np.string_('degrees')

        if isVaried(a4):
            dset = sam.create_dataset('polar_angle',data=a4)
            dset.attrs['target'] = np.string_('/entry/CAMEA/analyzer/polar_angle')
            nxdata['polar_angle'] = dset
            scanType = 'cscan a4 {} da4 {} np {} mn 10000'.format(np.mean(a4),np.mean(np.diff(a4)),len(a4))
            scanvars+='a4'
        else:
            #dset = sam.create_dataset('polar_angle',(1,),dtype='float32',data=-a4[0])
            #dset_zero = sam.create_dataset('polar_angle_offset',(1,),dtype='float32',data=polar_angle_offset)

            dset = entry['CAMEA/analyzer'].create_dataset('polar_angle',(1,),dtype='float32',data=a4[0])
        dset_zero = entry['CAMEA/analyzer'].create_dataset('polar_angle_offset',(1,),dtype='float32',data=polar_angle_offset)

        dset.attrs['units'] = np.string_('degrees')
        dset_zero.attrs['units'] = np.string_('degrees')
        

        mono = entry['CAMEA/monochromator']
        theta,tth = makeTheta(ei)

        if isVaried(ei):
            dset = mono.create_dataset('energy',data=ei)
            dset.attrs['units'] = np.string_('meV')
            dset.attrs['target'] = np.string_('/entry/CAMEA/monochromator/energy')
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

        #dset = mono.create_dataset('summed_counts',data=np.sum(data,axis=(1,2)));
        #dset.attrs['units'] = np.string_('counts')
        
        
        makeMonitor(entry,Numpoints)
        entry.create_dataset('scancommand',data=[np.string_(scanType)])
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

    with hdf.File(fileName,'w') as f:
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
        
        addSample(entry,np.string_(sample),cell,ub,plane_vector_1,plane_vector_2,plane_normal)
        import os
        Numpoints = sum([os.path.isdir(fname+'/'+i) for i in os.listdir(fname)])
        data,a3,a4,ei = readScanData(fname,Numpoints,factor=factor,pixels=pixels,detectors=detectors)
        storeScanData(entry,data,a3,a4,ei,rotation_angle_zero=rotation_angle_zero,polar_angle_offset=polar_angle_offset)
        

def prediction(A3Start,A3Stop,A3Steps,A4Positions,Ei,Cell,r1,r2,points=False):
    s = Sample.Sample(*Cell,projectionVector1=r1,projectionVector2=r2)
 
    class simpleDataSet():
        def __init__(self,sample):
            self.sample = [sample]

   
    def converterToQxQy(A3,A4,Ei,Ef):

        ki = np.sqrt(Ei)*_tools.factorsqrtEK
        kf = np.sqrt(Ef)*_tools.factorsqrtEK

        r = [0,0,0,A3,A4,0.0,0.0,Ei,Ef]
        QV = TasUBlib.calcTasUVectorFromAngles(r)
        q = np.sqrt(ki**2 +kf**2-
            2. *ki *kf * np.cos(np.deg2rad(A4)))

        
        Qx,Qy = (QV*q.reshape(-1,1))[:,:2].T
        return (Qx,Qy)

    def converterToA3A4(Qx,Qy, Ei,Ef,A3Off=0.0):
        Qx = np.asarray(Qx)
        Qy = np.asarray(Qy)

        QC = np.array([Qx,Qy])
        q = np.linalg.norm(QC)

        U1V = np.array([Qx.flatten(),Qy.flatten(),0.0],dtype=float)

        U1V/=np.linalg.norm(U1V)
        U2V = np.array([0.0,0.0,1.0],dtype=float)
        
        
        TV = TasUBlib.buildTVMatrix(U1V, U2V)
        R = np.linalg.inv(TV)
        
        ss = 1.0
        
        cossgl = np.sqrt(R[0,0]*R[0,0]+R[1,0]*R[1,0])
        om = TasUBlib.arctan2d(R[1,0]/cossgl, R[0,0]/cossgl)
     
        ki = np.sqrt(Ei)*_tools.factorsqrtEK
        kf = np.sqrt(Ef)*_tools.factorsqrtEK
        
        cos2t =(ki**2 + kf**2 - q**2) / (2. * np.abs(ki) * np.abs(kf))
        
        A4 = ss*TasUBlib.arccosd(cos2t)
        theta = TasUBlib.calcTheta(ki, kf, A4)
        A3 = -om + -ss*theta + A3Off
        return A3,-A4
        

    t = simpleDataSet(s)
    fig = plt.figure(figsize=(26,13))
    Ax = [RLUAxes.createRLUAxes(t,figure=fig,ids=(3,3,i+1)) for i in range(9)]
    
    MJOLNIRDir = os.path.dirname(_tools.__file__)
    calibPath = os.path.join(MJOLNIRDir,'Normalization_1.calib')

    calib = np.loadtxt(calibPath,delimiter=',',skiprows=3)

    A4Instrument = calib[:,-1].reshape(104,8)
    EfInstrument = calib[:,4].reshape(104,8)
    EfInstrument[np.isclose(EfInstrument,0.0)]=np.nan

    A3PointValues = np.linspace(A3Start,A3Stop,A3Steps)-s.theta
    steps = A3Steps
    if steps>30:
        steps = 30
    A3Values = np.linspace(A3Start,A3Stop,steps)-s.theta

    endColor = np.array([0.12156863, 0.46666667, 0.70588235, 0.5])
    startColor = np.array([0.83921569, 0.15294118, 0.15686275, 0.5])
    diff = (endColor-startColor)/7.0

    def format_axes_func(self,x,y,Ei,Ef,A3Off=s.theta):
        A3,A4 = converterToA3A4(x,y,Ei,Ef,-A3Off)
        return ax.sample.format_coord(x,y)+', A3 = {:.2f} deg, A4 = {:.2f} deg, Ef = {:.2f}'.format(A3,A4,Ef)

    for i,(ax,Ef,A4) in enumerate(zip(Ax,EfInstrument.reshape(104,8).T,A4Instrument.reshape(104,8).T)):
        ax.Ef = np.nanmean(Ef)
        for A4Position in A4Positions:
            if points:
                for A3 in A3PointValues:
                    H,L = converterToQxQy(A3=A3, A4=A4+A4Position,Ei = Ei, Ef = Ef)
                    ax.scatter(H,L,color=startColor[:-1]+diff[:-1]*i)

            A4Edge = np.array([np.concatenate([a4,[a4[-1]]*steps,a4[::-1],[a4[0]]*steps]) for a4 in A4.reshape(8,13)])
            EfEdge = np.array([np.concatenate([ef,[ef[-1]]*steps,ef[::-1],[ef[0]]*steps]) for ef in Ef.reshape(8,13)])
            A3Edge = np.array(list(np.concatenate([[A3Values[0]]*13,A3Values,[A3Values[-1]]*13,A3Values[::-1]]))*8).reshape(8,-1)
            Hedge,Ledge = np.array([converterToQxQy(A3=a3, A4=a4+A4Position,Ei = Ei, Ef = ef) for a3,a4,ef in zip(A3Edge,A4Edge,EfEdge)]).transpose([1,0,2])
            [a.plot(hedge,ledge,color=startColor+diff*i) for hedge,ledge in zip(Hedge,Ledge) for a in [Ax[-1],ax]]

            
        ax.set_xlabel(_tools.generateLabelDirection(s.projectionVector1) + ' RLU')
        ax.set_ylabel(_tools.generateLabelDirection(s.projectionVector2) + ' RLU')
        ax.set_title(r'$\Delta $E {:.2f}'.format(Ei-np.nanmean(Ef))+' meV, E$_f$ {:.2f} meV'.format(np.nanmean(Ef)))
        

    Ax[-1].set_xlabel(_tools.generateLabelDirection(s.projectionVector1) + ' RLU')
    Ax[-1].set_ylabel(_tools.generateLabelDirection(s.projectionVector2) + ' RLU')
    Ax[-1].set_title('All Planes')

    fig.tight_layout()
    Ax[0].format_coord = lambda x,y: format_axes_func(Ax[0],x,y,Ei,np.nanmean(EfInstrument,axis=0)[0],-s.theta)
    Ax[1].format_coord = lambda x,y: format_axes_func(Ax[1],x,y,Ei,np.nanmean(EfInstrument,axis=0)[1],-s.theta)
    Ax[2].format_coord = lambda x,y: format_axes_func(Ax[2],x,y,Ei,np.nanmean(EfInstrument,axis=0)[2],-s.theta)
    Ax[3].format_coord = lambda x,y: format_axes_func(Ax[3],x,y,Ei,np.nanmean(EfInstrument,axis=0)[3],-s.theta)
    Ax[4].format_coord = lambda x,y: format_axes_func(Ax[4],x,y,Ei,np.nanmean(EfInstrument,axis=0)[4],-s.theta)
    Ax[5].format_coord = lambda x,y: format_axes_func(Ax[5],x,y,Ei,np.nanmean(EfInstrument,axis=0)[5],-s.theta)
    Ax[6].format_coord = lambda x,y: format_axes_func(Ax[6],x,y,Ei,np.nanmean(EfInstrument,axis=0)[6],-s.theta)
    Ax[7].format_coord = lambda x,y: format_axes_func(Ax[7],x,y,Ei,np.nanmean(EfInstrument,axis=0)[7],-s.theta)
    fig.tight_layout()

    return Ax