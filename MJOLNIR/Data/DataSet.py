import sys, os
sys.path.append('.')
sys.path.append('..')
sys.path.append('../..')
import scipy
from scipy.ndimage import filters
import matplotlib.pyplot as plt
import numpy as np
import pickle as pickle
from MJOLNIR.Geometry import GeometryConcept,Analyser,Detector,Wedge,Instrument
import h5py as hdf
import scipy.optimize
import datetime
import warnings

dataLocation = 'entry/data/data'#'entry/Detectors/Detectors'
EiLocation = 'entry/data/incident_energy' # 'entry/Ei'
monLocation = 'entry/control/data'#'entry/Monitor'


NumberOfSigmas= 3 # Defining the active area of a peak on a detector as \pm n*sigma



class DataSet(object):
    def __init__(self, instrument=None,datafiles=None,normalizationfiles=None, calibrationfiles=None, convertedfiles=None, **kwargs):
        """DataSet object to hold all informations about data.
        
        Kwargs:
            
            - instrument (Instrument): Instrument object describing the data (default None).

            - datafiles (list of strings): List of datafiles to be used in conversion.

            - normalizationfiles (string or list of strings): Location of Vanadium normalization file(s).

            - calibrationfiles (string or list of strings): Location of calibration normalization file(s).

            - convertedfiles (string or list of strings): Location of converted data files.

        Raises:

            - ValueError
            
            - NotImplementedError
        
        
        """
        
        self._instrument = None
        self._datafiles = []
        self._normalizationfiles = []
        self._convertedfiles = []
        self._calibrationfiles = []
        if instrument is not None:
            self.instrument = instrument


        if datafiles is not None:
            self.datafiles = datafiles

        if normalizationfiles is not None:
            self.normalizationfiles = normalizationfiles
        
        if convertedfiles is not None:
            self.convertedfiles = convertedfiles

        if calibrationfiles is not None:
            self.calibrationfiles = calibrationfiles


        self._settings = {}
            
        
        
        # Add all other kwargs to settings
        for key in kwargs:
            self.settings[key]=kwargs[key]
        

    @property
    def instrument(self):
        return self._instrument

    @instrument.getter
    def instrument(self):
        return self._instrument

    @instrument.setter
    def instrument(self,instrument):
        if not issubclass(type(instrument),Instrument.Instrument):
            raise AttributeError('Instrument provided is not of instrument type!')
        else:
            self._instrument = instrument


    @property
    def datafiles(self):
        return self._datafiles

    @datafiles.getter
    def datafiles(self):
        return self._datafiles

    @datafiles.setter
    def datafiles(self,datafiles):
        try:
            self._datafiles = IsListOfStrings(datafiles)
        except Exception as e:
            #print('Error {} while parsing input.!'.format(e))
            raise(e)


    @property
    def normalizationfiles(self):
        return self._normalizationfiles

    @normalizationfiles.getter
    def normalizationfiles(self):
        return self._normalizationfiles

    @normalizationfiles.setter
    def normalizationfiles(self,normalizationfiles):
        try:
            self._normalizationfiles = IsListOfStrings(normalizationfiles)
        except Exception as e:
            #print('Error {} while parsing input.!'.format(e))
            raise(e)


    @property
    def convertedfiles(self):
        return self._convertedfiles

    @convertedfiles.getter
    def convertedfiles(self):
        return self._convertedfiles

    @convertedfiles.setter
    def convertedfiles(self,convertedfiles):
        try:
            self._convertedfiles = IsListOfStrings(convertedfiles)
        except Exception as e:
            #print('Error {} while parsing input.!'.format(e))
            raise(e)


    @property
    def calibrationfiles(self):
        return self._calibrationfiles

    @calibrationfiles.getter
    def calibrationfiles(self):
        return self._calibrationfiles

    @calibrationfiles.setter
    def calibrationfiles(self,calibrationfiles):
        try:
            self._calibrationfiles = IsListOfStrings(calibrationfiles)
        except Exception as e:
            #print('Error {} while parsing input.!'.format(e))
            raise(e)

    @property
    def settings(self):
        return self._settings

    @settings.getter
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self,*args,**kwargs):
        raise NotImplementedError('Settings cannot be overwritten.')    

    def load(self,filename):
        """Method to load an object from a pickled file."""
        try:                                # Opening the given file with an error catch
            fileObject = open(filename, 'rb')
        except IOError as e:                        # Catch all IO-errors
            print("Error in opening file:\n{}".format(e))
        else:
            tmp_dict = pickle.load(fileObject)
            
            fileObject.close()
            # TODO: Make checks that the object loaded is of correct format?
            self=tmp_dict

    def save(self, filename):
        try:                                # Opening the given file with an error catch
            fileObject = open(filename, 'wb')
        except IOError as e:                        # Catch all IO-errors
            print("Error in opening file:\n{}".format(e))
        else:
                pickle.dump(self, fileObject)
                fileObject.close()  

    def initialize(self):
        ## Needed things for initialization: Instrument(?!), datafiles, normalizationfiles(?)
        if self.instrument is None:
            raise RuntimeError("DataSet object does not contain an instrument!")
        if len(self.datafiles) == 0:
            raise RuntimeError("DataSet object does not contain data files!")
        if len(self.normalizationfiles) ==0:
            raise RuntimeError("DataSet object does not contain normalization files!")

        if self.instrument.settings['Initialized']==False:
            self.instrument.initialize()

    def __eq__(self, other): 
        return np.logical_and(set(self.__dict__.keys()) == set(other.__dict__.keys()),self.__class__ == other.__class__)

    def __str__(self):
        string = '{} with settings:\n'.format(self.__class__)
        for attrib in self.settings:
            string+='{}:\t{}\n'.format(attrib,self.settings[attrib])
        for attrib in self.__dict__:
            if attrib=='_settings':
                continue
            string+='{}:\t{}\n'.format(attrib,self.__dict__[attrib])
        string+='\n'    
        return string


    def EnergyCalibration(self,datafile=None,savelocation='energycalibration/',tables=['Single','PrismaticLowDefinition','PrismaticHighDefinition'],InstrumentType='CAMEA',plot=False):
        """Method to generate look-up tables for normalization. Saves calibration file(s) as 'EnergyCalibration_Np.calib', where Np is the number of pixels.
        
        Generates 4 different tables:

            - Prismatic High Definition (8 pixels/energy or 64 pixels/detector)

            - Prismatic Low Definition (3 pixels/energy or 24 pixels/detector)

            - Single (1 pixel/energy or 8 pixels/detector)

            - Number (integer)

        Kwargs:

            - datafile (string): String to data single data file used for normalization (required).

            - savelocation (string): String to save location folder (energycalibration)

            - tables (list): List of needed conversion tables (Default: ['Single','PrismaticLowDefinition','PrismaticHighDefinition','Unbinned'], increasing number of pixels).

            - InstrumentType (string): Type of instrument (default CAMEA).

            - plot (boolean): Set to True if pictures of all fit are to be stored in savelocation


        .. note::
            Assumes that file contains a scan over  Ei. If one is in doubt whether the Vanadium normalization data has enough data points to allow a full 8 pixel binning, start with 1 pixel and 3, as the 8 pixel binning may end up raising an error.
        

        .. warning::
            At the moment, the active detector area is defined by NumberOfSigmas (currently 3) times the Guassian width of Vanadium peaks.

        """
        if self.instrument.settings['Initialized']==False:
            self.initialize()

        if datafile is None:
            datafile = self.normalizationfiles[0]
            print('Using {} for normalization table.'.format(datafile))

        if InstrumentType=='CAMEA':

            NormFile = hdf.File(datafile)

            if savelocation[-1]!='/':
                savelocation+='/'
            
            Data = np.array(NormFile.get(dataLocation)).transpose(1,0,2)
            Ei = np.array(NormFile.get(EiLocation))
            analysers = 8
            pixels = len(self.instrument.A4[0][0]) # <------------------- Change!
            detectors = len(self.instrument.A4[0])*len(self.instrument.A4)
            detectorsorInWedge = len(self.instrument.A4[0])
            wedges = len(self.instrument.A4)

            if pixels!=Data.shape[2]:
                raise ValueError('The number of pixels ({}) in the data file does not match instrument description ({})!'.format(pixels,Data.shape[2]))

            bins = []
            for table in tables:
                if table=='Unbinned':
                    bins.append(pixels)
                elif table=='PrismaticHighDefinition':
                    bins.append(8)
                elif table=='PrismaticLowDefinition':
                    bins.append(3)
                elif table=='Single':
                    bins.append(1)
                elif isinstance(table,int):
                    bins.append(table)
                else:
                    raise AttributeError("Provided table attribute ({}) not recognized, should be 'Unbinned','PrismaticHighDefinition','PrismaticLowDefinition','Single', and/or integer.".format(table))
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

            

            if plot:
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
                    guess = [peakVal[i,j],float(peakPos[i,j]),pixels/100.0,np.min(ESummedData[i])]
                    #if i ==7:
                    #    print(guess)
                    res = scipy.optimize.curve_fit(Gaussian,np.arange(ESummedData.shape[1]),dataSubtracted[i,:],p0=[guess])
                    #if np.abs(res[0][1]-guess[1])>10:
                    #    print('Peak in detector {} with analyser {}, guess {}, diff {}'.format(i,j,guess[1],np.abs(res[0][1]-guess[1])))
                    peakPos[i,j] = res[0][1]
                    peakVal[i,j] = res[0][0]
                    peakWidth[i,j]= res[0][2]

                    # Generate peak as the one fitted and subtract it from signal
                    x=np.arange(pixels)
                    y = Gaussian(x,peakVal[i,j],peakPos[i,j],peakWidth[i,j],peakBackg[i,j])
                    #dataSubtracted[i]-=y
                    peak = y>peakVal[i,j]*0.05
                    dataSubtracted[i,peak]= 0

            if plot:
                plt.clf()
                #figman = plt.get_current_fig_manager()
                #figman.window.setGeometry(1366, 24, 1920, 1176)
                plt.suptitle('Fits')
                x = np.arange(pixels)
                for k in range(wedges):
                    for i in range(detectorsorInWedge):
                        y=np.zeros_like(x,dtype=float)
                        plt.subplot(4, 4, i+1)
                        plt.scatter(np.arange(pixels),ESummedData[i+13*k],s=4)
                        #    plt.scatter(peakPos[i],peakVal[i])
                        for j in range(analysers):
                            y += Gaussian(x,peakVal[i+13*k,j],peakPos[i+13*k,j],peakWidth[i+13*k,j],peakBackg[i+13*k,j])
                            plt.plot([peakPos[i+13*k,j],peakPos[i+13*k,j]],[0,np.max(ESummedData[i+13*k])*1.1])
                        plt.plot(x,y,'k')
                        plt.xlabel('Pixel')
                        plt.ylabel('Intensity [arg]')
                        plt.title('Detector {}'.format(i))
                        #plt.legend(['Fit','Data'])
                        plt.ylim(0,np.max(ESummedData[i+13*k])*1.1)

                    plt.tight_layout()
                    plt.savefig(savelocation+'/Raw/Fit_wedge_'+str(k)+'.png',format='png', dpi=150)

            
            ## Sort the positions such that peak 1 is the furthermost left peak and assert diff(pos)>100
            sortedPeakPosArg = np.argsort(peakPos,axis=1)
            sortedPeakPos = np.sort(peakPos,axis=1)
            sortedPeakPos[np.logical_or(sortedPeakPos>pixels,sortedPeakPos<0)]=5*pixels # High number

            sortedPeakPosArg2 = np.argsort(sortedPeakPos,axis=1)
            sortedPeakPos.sort(axis=1)

            differences = np.diff(sortedPeakPos,axis=1)
            outliers = np.zeros_like(peakPos,dtype=bool)
            outliers[:,:-1]=differences<pixels/10
            sortedPeakPos[outliers]=5*pixels
            sortedPeakPosArg3 = np.argsort(sortedPeakPos,axis=1)
            argSort = np.array([sortedPeakPosArg[i,sortedPeakPosArg2[i,sortedPeakPosArg3[i,:]]] for i in range(detectors)])
            sortedPeakPos = np.sort(sortedPeakPos,axis=1)
            peaks=np.sum(sortedPeakPos<5*pixels,axis=1) # Number of peaks found


            if np.any(peaks!=analysers):
                raise ValueError('Wrong number of peaks, {} found in detector(s): {}\nIn total error in {} detector(s).'.format(peaks[peaks!=analysers],np.arange(peaks.shape[0])[peaks!=analysers],np.sum(peaks[peaks!=analysers])))

            pixelpos  = np.array([peakPos[i,argSort[i]] for i in range(detectors)])
            #amplitudes= np.array([peakVal[i,argSort[i]] for i in range(detectors)])
            widths    = np.array([peakWidth[i,argSort[i]] for i in range(detectors)])
            #backgrounds=np.array([peakBackg[i,argSort[i]] for i in range(detectors)])
            
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
                if plot:
                    plt.clf()
                    plt.title('Detector {} Active pixels'.format(i))
                    plt.scatter(x,ESummedData[i],s=4,color='black')
                for j in range(analysers):
                    activePixels[i,j] = np.logical_and(x>lowerPixel[i,j],x<upperPixel[i,j])
                    if plot: plt.scatter(x[np.logical_and(x>lowerPixel[i,j],x<upperPixel[i,j])],
                        ESummedData[i,np.logical_and(x>lowerPixel[i,j],x<upperPixel[i,j])],s=4,color='red')
                if plot:
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
                    warnings.warn('Fitting might be unstable due to {} pixels being fitted using only {} energies ({} free parameters).'.format(detpixels,len(Ei),detpixels*analysers*3))
                    
                if plot:
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
                    if plot:
                        plt.clf()
                        plt.title('Detector {}, {} pixels'.format(i,detpixels))
                        x=np.linspace(0,detpixels,len(Ei))
                    for j in range(analysers):

                        center = int(round(sortedPeakPos[i,j]))
                        width = activePixels[i,j].sum()
                        pixelAreas = np.linspace(-width/2.0,width/2.0,detpixels+1,dtype=int)+center+1 #Add 1 such that the first pixel is included 20/10-17
                        
                        for k in range(detpixels):
                            binPixelData = Data[i,:,pixelAreas[k]:pixelAreas[k+1]].sum(axis=1)
                            #print('{},{}'.format(pixelAreas[k],pixelAreas[k+1]))
                            #[print(x) for x in binPixelData]
                            guess = [np.max(binPixelData), Ei[np.argmax(binPixelData)],0.1,0]
                            try:
                                res = scipy.optimize.curve_fit(Gaussian,Ei,binPixelData,p0=guess)
                            except:
                                if not os.path.exists(savelocation+'/{}_pixels'.format(detpixels)):
                                    os.makedirs(savelocation+'/{}_pixels'.format(detpixels)) 
                                plt.figure()
                                plt.scatter(Ei,binPixelData)
                                plt.plot(Ei,Gaussian(Ei,guess[0],guess[1],guess[2],guess[3]))
                                plt.savefig(savelocation+'/{}_pixels/Detector{}_{}.png'.format(detpixels,i,k),format='png',dpi=300)
                                plt.close()
                                #raise ValueError('Fitting not converged, probably due to too few points.')



                            fittedParameters[i,j,k]=res[0]
                            if plot:
                                plt.plot(EiX,Gaussian(EiX,*fittedParameters[i,j,k]),color='black')
                                plt.scatter(Ei,binPixelData,color=colors[:,k])
                        activePixelAnalyser.append(np.linspace(-width/2.0,width/2.0,detpixels+1,dtype=int)+center+1)
                    activePixelDetector.append(activePixelAnalyser)
                    if plot:
                        plt.grid('on')
                        plt.xlabel('Ei [meV]')
                        plt.ylabel('Weight [arb]')
                        plt.tight_layout(rect=(0,0,1,0.95))
                        plt.savefig(savelocation+'/{}_pixels/Detector{}.png'.format(detpixels,i),format='png',dpi=300)


                fitParameters.append(fittedParameters)
                activePixelRanges.append(np.array(activePixelDetector))
                tableString = 'Normalization for {} pixel(s) using data {}\nPerformed {}\nDetector,Energy,Pixel,Amplitude,Center,Width,Background,lowerBin,upperBin\n'.format(detpixels,datafile,datetime.datetime.now())
                for i in range(len(fittedParameters)):
                    for j in range(len(fittedParameters[i])):
                        for k in range(len(fittedParameters[i][j])):
                            tableString+=str(i)+','+str(j)+','+str(k)+','+','.join([str(x) for x in fittedParameters[i][j][k]])
                            tableString+=','+str(activePixelRanges[-1][i][j][k])+','+str(activePixelRanges[-1][i][j][k+1])+'\n'
                
                tableName = 'EnergyNormalization_{}.calib'.format(detpixels)
                print('Saving {} pixel data to {}'.format(detpixels,savelocation+tableName))
                file = open(savelocation+tableName,mode='w')

                file.write(tableString)
                file.close()
                self.calibrationfiles.append(savelocation+tableName)




    def ConvertDatafile(self,datafiles=None,calibrationfile=None,savelocation=None):
        """Conversion method for converting scan file(s) to hkl file. Converts the given h5 file into NXqom format and saves in a file with same name, but of type .nxs.
        Copies all of the old data file into the new to ensure complete reduncency. Determins the binning wanted from the file name of normalization file.

        Kwargs:

            - datafiles (string or list of): File path(s), file must be of hdf format (default self.datafiles).

            - calibrationfile (string): File path to normalization file (default self.calibrationfiles[-1]).

            - savelocation (string): File path to save location of data file(s) (defaults to same as raw file).

        Raises:

            - IOError

            - AttributeError
            
        """

        if calibrationfile is None:
            binnings = []
            for nfile in self.calibrationfiles:
                binnings.append(int(nfile.split('_')[-1].split('.')[0]))
            if len(binnings)==0:
                raise AttributeError('No normalization file provided either through input of in the DataSet object.')
            else:
                maxarg = np.argmax(binnings)
                calibrationfile = self.calibrationfiles[maxarg]
                print('Using normalization file {}'.format(calibrationfile))

        binning = int(calibrationfile.split('_')[-1].split('.')[0]) # Find binning from normalization file name
        

        if datafiles is None:
            if len(self.datafiles)==0:
                raise AttributeError('No data files file provided either through input of in the DataSet object.')
            datafiles = self.datafiles


        if not isinstance(datafiles,list):
            datafiles=[datafiles]
        for datafile in datafiles:
            normalization = np.genfromtxt(calibrationfile,delimiter=',',skip_header=3)
            EPrDetector = len(self.instrument.wedges[0].detectors[0].split)+1
            
            
            file = hdf.File(datafile)
        
            A4 = np.array(self.instrument.A4)
            A4=A4.reshape(A4.shape[0]*A4.shape[1],A4.shape[2],order='C')
            Ef = np.array(self.instrument.Ef)
            Ef=Ef.reshape(Ef.shape[0]*Ef.shape[1],Ef.shape[2],order='C')

            PixelEdge = normalization[:,7:].reshape(A4.shape[0],EPrDetector,binning,2).astype(int)
            
            Data = np.array(file.get('/entry/data/data'))
            
            #ScanPoints = Data.shape[0]
            
            if not Data.shape[1:]==A4.shape:
                raise AttributeError('The shape of the data{} does not match instrument{}!'.format(Data.shape[1:],A4.shape))
            
            
            factorsqrtEK = 0.694692
            
            A4File = np.array(file.get('entry/CAMEA/detector/polar_angle'))
            A4Shape = A4.shape
            A4Total = -A4.reshape((1,A4Shape[0],A4Shape[1]))-np.deg2rad(A4File).reshape((A4File.shape[0],1,1))#-np.deg2rad(2.3173119802914783)
            #PixelEdgeA4Shaped = PixelEdge.reshape((1,PixelEdge.shape[0],EPrDetector,binning,2))
            
            A4Mean = np.zeros((A4Total.shape[0],A4Total.shape[1],EPrDetector*binning))
            #EfMean = np.zeros((Ef.shape[0],EPrDetector*binning))
            DataMean=np.zeros((Data.shape[0],Data.shape[1],EPrDetector*binning),dtype=int)
            for i in range(A4Total.shape[1]): # for each detector
                for j in range(EPrDetector):
                    for k in range(binning):
                        A4Mean[:,i,j*binning+k] = np.mean(A4Total[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)
                        DataMean[:,i,j*binning+k] = np.sum(Data[:,i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=1)
                        #EfMean[i,j*binning+k] = np.mean(Ef[i,PixelEdge[i,j,k,0]:PixelEdge[i,j,k,1]],axis=0)
            
            EfMean = normalization[:,4].reshape(A4.shape[0],EPrDetector*binning)
            Normalization = (normalization[:,3]*np.sqrt(2*np.pi)*normalization[:,5]).reshape(A4.shape[0],EPrDetector*binning)

            kf = factorsqrtEK*np.sqrt(EfMean)
            Ei = np.array(file.get('/entry/CAMEA/monochromator/energy'))
            
            ki = factorsqrtEK*np.sqrt(Ei)
            
            A3 = np.deg2rad(np.array(file.get('/entry/sample/rotation_angle/')))
            
            # Qx = ki-kf*cos(A4), Qy = -kf*sin(A4)

            #Qx = ki.reshape((*ki.shape,1,1,1))-(kf.reshape((1,1,*Ef.shape))*np.cos(A4Total)).reshape((1,*A4Total.shape))
            Qx = ki.reshape((ki.shape[0],1,1,1))-(kf.reshape((1,1,EfMean.shape[0],EfMean.shape[1]))*np.cos(A4Mean)).reshape((1,A4Mean.shape[0],A4Mean.shape[1],A4Mean.shape[2]))
            #Qy = np.zeros((*ki.shape,1,1,1))-kf.reshape((1,1,*Ef.shape))*np.sin(A4Total.reshape((1,*A4Total.shape)))
            Qy = np.zeros((ki.shape[0],1,1,1))-kf.reshape((1,1,EfMean.shape[0],EfMean.shape[1]))*np.sin(A4Mean.reshape((1,A4Mean.shape[0],A4Mean.shape[1],A4Mean.shape[2])))
            
            QX = Qx.reshape((1,Qx.shape[0],Qx.shape[1],Qx.shape[2],Qx.shape[3]))*np.cos(A3.reshape((A3.shape[0],1,1,1,1)))-Qy.reshape((1,Qy.shape[0],Qy.shape[1],Qy.shape[2],Qy.shape[3]))*np.sin(A3.reshape((A3.shape[0],1,1,1,1)))
            QY = Qx.reshape((1,Qx.shape[0],Qx.shape[1],Qx.shape[2],Qx.shape[3]))*np.sin(A3.reshape((A3.shape[0],1,1,1,1)))+Qy.reshape((1,Qy.shape[0],Qy.shape[1],Qy.shape[2],Qy.shape[3]))*np.cos(A3.reshape((A3.shape[0],1,1,1,1)))
            if QX.shape.count(1)!=2:
                raise ValueError('At least two parameters changed simulatneously!')
            
            #EnergyShape = (1,len(Ei),1,*Ef.shape)
            EnergyShape = (1,len(Ei),1,EfMean.shape[0],EfMean.shape[1])
            #DeltaE = (Ei.reshape((*Ei.shape,1,1))-Ef.reshape((1,*Ef.shape))).reshape(EnergyShape)
            DeltaE = (Ei.reshape((Ei.shape[0],1,1))-EfMean.reshape((1,EfMean.shape[0],EfMean.shape[1]))).reshape(EnergyShape)
            
            
            
            #Intensity = Data.reshape(*QX.shape)
            Intensity = DataMean.reshape((QX.shape[0],QX.shape[1],QX.shape[2],QX.shape[3],QX.shape[4]))
        
            DeltaE=DeltaE.repeat(QX.shape[0],axis=0)
            DeltaE=DeltaE.repeat(QX.shape[2],axis=2)
            
            Monitor = np.array(file.get('/entry/control/data'),dtype=int)
            Monitor.shape = (len(Monitor),1,1)
            Monitor = np.repeat(Monitor,Intensity.shape[3],axis=1)
            Monitor = np.repeat(Monitor,Intensity.shape[4],axis=2)
            Monitor.shape = Intensity.shape

            Normalization.shape = (1,1,1,Normalization.shape[0],Normalization.shape[1])
            Normalization = np.repeat(Normalization,Intensity.shape[0],axis=0)
            Normalization = np.repeat(Normalization,Intensity.shape[1],axis=1)
            Normalization = np.repeat(Normalization,Intensity.shape[2],axis=2)
            ## TODO: Don't let all things vary at the same time!!
            
            if not savelocation is None:
                saveloc = savelocation+datafile.replace('.h5','.nxs').split('/')[-1]
            else:
                saveloc = datafile.replace('.h5','.nxs')
            saveNXsqom(datafile,file,saveloc,Intensity,Monitor,QX,QY,DeltaE,calibrationfile,Normalization)
            
            file.close()
            self.convertedfiles.append(saveloc)

    def binData3D(self,dx,dy,dz,datafiles=None):
        """Bin a converted data file into voxels with sizes dx*dy*dz. Wrapper for the binData3D functionality.

        Args:

            - dx (float): step sizes along the x direction (required).

            - dy (float): step sizes along the y direction (required).

            - dz (float): step sizes along the z direction (required).

        Kwargs:

            - datafile (string or list of strings): Location(s) of data file to be binned (default converted file in DataSet).

        Raises:

            - AttributeError

        Returns:

            - Datalist: List of converted data files having 4 sub arrays: Intensity(counts), Monitor, Normalization, Normalization count
            - bins: 3 arrays containing edge positions in x, y, and z directions.
        """
        
        if datafiles is None:
            if len(self.convertedfiles)==0:
                raise AttributeError('No data file to be binned provided in either input or DataSet object.')
            else:
                datafiles = self.convertedfiles
        elif not isinstance(datafiles,list):
            datafiles = [datafiles]
        else:
            raise AttributeError('datafiles attribute not understood.')

        returnData = []

        for data in datafiles:
            
            file = hdf.File(data,'r')

            I = np.array(file.get('entry/data/data'))
            posx = np.array(file.get('entry/data/qx'))
            posy = np.array(file.get('entry/data/qy'))
            energy = np.array(file.get('entry/data/en'))
            Norm = np.array(file.get('entry/data/normalization'))
            Monitor = np.array(file.get('entry/data/monitor'))
            file.close()

            pos = [posx,posy,energy]
            _tempData,bins = binData3D(dx,dy,dz,pos,I,norm=Norm,mon=Monitor)
            returnData.append(_tempData)

        return returnData,bins
        





def Gaussian(x,A,mu,sigma,b):
    #A,mu,sigma,b = args
    return A*np.exp(-np.power(mu-x,2.0)*0.5*np.power(sigma,-2.0))+b


def findPeak(data):
    return [np.argmax(data,axis=1),np.max(data,axis=1)]



        
def IsListOfStrings(object):
    if isinstance(object, list):
        isListOfStr = True
        for item in object:
            if not isinstance(item, str):
                isListOfStr=False
                break
        if isListOfStr:
            return object
        else:
            raise AttributeError('Data files provided are not a list of strings or string!')
    elif isinstance(object,str):
        return [object]
    else:
        raise AttributeError('Data files provided are not a list of strings or string!')
    


def saveNXsqom(datafile,fs,savefilename,Intensity,Monitor,QX,QY,DeltaE,normalizationfile,Normalization):
    
    fd = hdf.File(savefilename,'w')
    group_path = fs['/entry'].parent.name
    
    group_id = fd.require_group(group_path)
    
    
    fs.copy('/entry', group_id, name="/entry")
    
    definition = fd.create_dataset('entry/definition',(1,),dtype='S70',data=np.string_('NXsqom'))
    definition.attrs['NX_class'] = 'NX_CHAR'
    
    process = fd.create_group('entry/reduction')
    process.attrs['NX_class']=b'NXprocess'
    proc = process.create_group('MJOLNIR_algorithm_convert')
    proc.attrs['NX_class']=b'NXprocess'
    author= proc.create_dataset('author',shape=(1,),dtype='S70',data=np.string_('Jakob Lass'))
    author.attrs['NX_class']=b'NX_CHAR'
    
    date= proc.create_dataset('date',shape=(1,),dtype='S70',data=np.string_(datetime.datetime.now()))
    date.attrs['NX_class']=b'NX_CHAR'
    
    description = proc.create_dataset('description',shape=(1,),dtype='S70',data=np.string_('Conversion from pixel to Qx,Qy,E in reference system of instrument.'))
    description.attrs['NX_class']=b'NX_CHAR'
    
    rawdata = proc.create_dataset('rawdata',shape=(1,),dtype='S200',data=np.string_(datafile))
    rawdata.attrs['NX_class']=b'NX_CHAR'
    
    normalizationString = proc.create_dataset('normalization table',shape=(1,),dtype='S200',data=np.string_(normalizationfile))
    normalizationString.attrs['NX_class']=b'NX_CHAR'
    
    data = fd.get('entry/data')
    data['rawdata']=data['data']
    del data['data']
    
    
    fileLength = Intensity.size
    
    Int = data.create_dataset('data',shape=(fileLength,),dtype='int32',data=Intensity.flatten())
    Int.attrs['NX_class']='NX_INT'
    
    monitor = data.create_dataset('monitor',shape=(fileLength,),dtype='int32',data=Monitor.flatten())
    monitor.attrs['NX_class']=b'NX_INT'

    normalization = data.create_dataset('normalization',shape=(fileLength,),dtype='float32',data=Normalization.flatten())
    normalization.attrs['NX_class']=b'NX_FLOAT'
    
    qx = data.create_dataset('qx',shape=(fileLength,),dtype='float32',data=QX.flatten())
    qx.attrs['NX_class']=b'NX_FLOAT'
    
    qy = data.create_dataset('qy',shape=(fileLength,),dtype='float32',data=QY.flatten())
    qy.attrs['NX_class']=b'NX_FLOAT'
    
    qz = data.create_dataset('qz',shape=(fileLength,),dtype='float32',data=np.zeros((fileLength,)))
    qz.attrs['NX_class']=b'NX_FLOAT'
    
    en = data.create_dataset('en',shape=(fileLength,),dtype='float32',data=DeltaE.flatten())
    en.attrs['NX_class']=b'NX_FLOAT'

    fd.close()


def calculateGrid3D(X,Y,Z):
    """Generate 3D grid with centers given by X,Y, and Z.
     Args:
        
        X (3D array): 3D array of x values generated by np.meshgrid.
                
        Y (3D array): 3D array of y values generated by np.meshgrid.
                
        Z (3D array): 3D array of z values generated by np.meshgrid.
        
    Example:

    >>> x = np.linspace(-1.5,1.5,20)
    >>> y = np.linspace(0,1.5,10)
    >>> z = np.linspace(-1.0,5.5,66)
    >>> X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
    >>> XX,YY,ZZ = calculateGrid3D(X,Y,Z)

    Now XX is a 21x11x67 array containing all x coordinates of the edges exactly midway bewteen the points. Same goes for YY and ZZ with y and z coordinates respectively.
    """

    xshape = X.shape
    
    XT = np.zeros((xshape[0]+1,xshape[1]+1,xshape[2]+1))
    YT = np.zeros_like(XT)
    ZT = np.zeros_like(XT)
    
    
    
    dx0 = np.diff(X,axis=0)
    dx1 = np.diff(X,axis=1)
    dx2 = np.diff(X,axis=2)
    dy0 = np.diff(Y,axis=0)
    dy1 = np.diff(Y,axis=1)
    dy2 = np.diff(Y,axis=2)
    dz0 = np.diff(Z,axis=0)
    dz1 = np.diff(Z,axis=1)
    dz2 = np.diff(Z,axis=2)
    
    
    XX = X.copy()
    XX[:-1]-=0.5*dx0
    XX[-1]-=0.5*dx0[-1]
    XX[:,:-1]-=0.5*dx1
    XX[:,-1]-=0.5*dx1[:,-1]
    XX[:,:,:-1]-=0.5*dx2
    XX[:,:,-1]-=0.5*dx2[:,:,-1]
    
    YY = Y.copy()
    YY[:-1]-=0.5*dy0
    YY[-1]-=0.5*dy0[-1]
    YY[:,:-1]-=0.5*dy1
    YY[:,-1]-=0.5*dy1[:,-1]
    YY[:,:,:-1]-=0.5*dy2
    YY[:,:,-1]-=0.5*dy2[:,:,-1]
    
    ZZ = Z.copy()
    ZZ[:-1]-=0.5*dz0
    ZZ[-1]-=0.5*dz0[-1]
    ZZ[:,:-1]-=0.5*dz1
    ZZ[:,-1]-=0.5*dz1[:,-1]
    ZZ[:,:,:-1]-=0.5*dz2
    ZZ[:,:,-1]-=0.5*dz2[:,:,-1]
    
    XT[:-1,:-1,:-1]=XX.copy()
    YT[:-1,:-1,:-1]=YY.copy()
    ZT[:-1,:-1,:-1]=ZZ.copy()
    
    
    XT[-1,:-1,:-1]=XT[-2,:-1,:-1]+dx0[-1]
    XT[:-1,-1,:-1]=XT[:-1,-2,:-1]+dx1[:,-1,:]
    XT[:-1,:-1,-1]=XT[:-1,:-1,-2]+dx2[:,:,-1]
    XT[:-1,-1,-1]=0.5*(XT[:-1,-1,-2]+dx2[:,-1,-1]+XT[:-1,-2,-1]+dx1[:,-1,-1])
    XT[-1,:-1,-1]=0.5*(XT[-1,:-1,-2]+dx2[-1,:,-1]+XT[-2,:-1,-1]+dx0[-1,:,-1])
    XT[-1,-1,:-1]=0.5*(XT[-1,-2,:-1]+dx1[-1,-1,:]+XT[-2,-1,:-1]+dx0[-1,-1,:])
    XT[-1,-1,-1]=(XT[-1,-2,-1]+dx1[-1,-1,-1]+XT[-2,-1,-1]+dx0[-1,-1,-1]+XT[-1,-1,-2]+dx2[-1,-1,-1])/3
    
    YT[-1,:-1,:-1]=YT[-2,:-1,:-1]+dy0[-1]
    YT[:-1,-1,:-1]=YT[:-1,-2,:-1]+dy1[:,-1,:]
    YT[:-1,:-1,-1]=YT[:-1,:-1,-2]+dy2[:,:,-1]
    YT[:-1,-1,-1]=0.5*(YT[:-1,-1,-2]+dy2[:,-1,-1]+YT[:-1,-2,-1]+dy1[:,-1,-1])
    YT[-1,:-1,-1]=0.5*(YT[-1,:-1,-2]+dy2[-1,:,-1]+YT[-2,:-1,-1]+dy0[-1,:,-1])
    YT[-1,-1,:-1]=0.5*(YT[-1,-2,:-1]+dy1[-1,-1,:]+YT[-2,-1,:-1]+dy0[-1,-1,:])
    YT[-1,-1,-1]=(YT[-1,-2,-1]+dy1[-1,-1,-1]+YT[-2,-1,-1]+dy0[-1,-1,-1]+YT[-1,-1,-2]+dy2[-1,-1,-1])/3
    
    ZT[-1,:-1,:-1]=ZT[-2,:-1,:-1]+dz0[-1]
    ZT[:-1,-1,:-1]=ZT[:-1,-2,:-1]+dz1[:,-1,:]
    ZT[:-1,:-1,-1]=ZT[:-1,:-1,-2]+dz2[:,:,-1]
    ZT[:-1,-1,-1]=0.5*(ZT[:-1,-1,-2]+dz2[:,-1,-1]+ZT[:-1,-2,-1]+dz1[:,-1,-1])
    ZT[-1,:-1,-1]=0.5*(ZT[-1,:-1,-2]+dz2[-1,:,-1]+ZT[-2,:-1,-1]+dz0[-1,:,-1])
    ZT[-1,-1,:-1]=0.5*(ZT[-1,-2,:-1]+dz1[-1,-1,:]+ZT[-2,-1,:-1]+dz0[-1,-1,:])
    ZT[-1,-1,-1]=(ZT[-1,-2,-1]+dz1[-1,-1,-1]+ZT[-2,-1,-1]+dz0[-1,-1,-1]+ZT[-1,-1,-2]+dz2[-1,-1,-1])/3
    
    
    return XT,YT,ZT




def binData3D(dx,dy,dz,pos,data,norm=None,mon=None):
    """ 3D binning of data.

    Args:

        - dx (float): Step size in x (required).

        - dy (float): Step size in x (required).

        - dz (float): Step size in x (required).

        - pos (2D array): Position of data points as flattened lists (X,Y,Z) (required).

        - data (array): Flattened data array (required).

    Kwargs:

        - norm (array): Flattened normalization array.

        - mon (array): Flattened monitor array.

    returns:

        Rebinned intensity (and if provided Normalization, Monitor, and Normalization Count) and X, Y, and Z bins in 3 1D arrays.


    Example:

    >>> pos = [Qx,Qy,E]
    >>> Data,bins = DataSet.binData3D(0.05,0.05,0.2,pos,I,norm=Norm,mon=Monitor)

    """

    

    
    diffx = np.abs(np.max(pos[0])-np.min(pos[0]))
    diffy = np.abs(np.max(pos[1])-np.min(pos[1]))
    diffz = np.abs(np.max(pos[2])-np.min(pos[2]))
    
    xbins = np.round(diffx/dx).astype(int)+1
    ybins = np.round(diffy/dy).astype(int)+1
    zbins = np.round(diffz/dz).astype(int)+1
    
    _X = np.linspace(np.min(pos[0]),np.max(pos[0]),xbins)
    _Y = np.linspace(np.min(pos[1]),np.max(pos[1]),ybins)
    _Z = np.linspace(np.min(pos[2]),np.max(pos[2]),zbins)
    
    X,Y,Z = np.meshgrid(_X,_Y,_Z,indexing='ij')
    
    XX,YY,ZZ = calculateGrid3D(X,Y,Z)
    
    bins=[XX[:,0,0],YY[0,:,0],ZZ[0,0,:]]
    
    

    intensity =    np.histogramdd(np.array(pos).T,bins=bins,weights=data.flatten())[0].astype(data.dtype)

    returndata = [intensity]
    if mon is not None:
        MonitorCount=  np.histogramdd(np.array(pos).T,bins=bins,weights=mon.flatten())[0].astype(mon.dtype)
        returndata.append(MonitorCount)
    if norm is not None:
        Normalization= np.histogramdd(np.array(pos).T,bins=bins,weights=norm.flatten())[0].astype(norm.dtype)
        NormCount =    np.histogramdd(np.array(pos).T,bins=bins,weights=np.ones_like(data).flatten())[0].astype(int)
        returndata.append(Normalization)
        returndata.append(NormCount)

    return returndata,bins

#__________________________TESTS_______________________

def test_DataSet_Creation():
    Instr = Instrument.Instrument()

    dataset = DataSet(instrument=Instr,OhterSetting=10.0)
    
    if(dataset.settings['OhterSetting']!=10.0):
        assert False

    try:
        dataset = DataSet(instrument=np.array([1.0]))
        assert False
    except:
        assert True

def test_Dataset_Initialization():
    Instr = Instrument.Instrument()
    emptyDataset = DataSet()
    try:
        emptyDataset.initialize()
        assert False
    except RuntimeError:
        assert True
    
    dataset = DataSet(instrument=Instr,OhterSetting=10.0,datafiles='SomeFile',normalizationfiles=['VanData','SecondFile'],convertedfiles='Converted.nxs')
    assert(dataset.datafiles[0]=='SomeFile')
    assert(dataset.normalizationfiles[0]=='VanData')
    assert(dataset.normalizationfiles[1]=='SecondFile')
    assert(len(dataset.normalizationfiles)==2)
    assert(dataset.instrument==Instr)
    assert(len(dataset.calibrationfiles)==0)
    assert(dataset.convertedfiles[0]=='Converted.nxs')


    try:
        dataset.initialize()
        assert False
    except:
        assert True

    wedge = Wedge.Wedge(position=(0.5,0,0),concept='ManyToMany')
    pixels=33
    split = [12,20]
    Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0),pixels=pixels,split=split)
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    wedge.append([Det,Det,Ana,Ana,Ana])
    Instr.append(wedge)

    dataset.initialize()

def test_DataSet_Error():
    try:
        ds = DataSet(instrument='Nope')
        assert False
    except AttributeError:
        assert True
    Instr = Instrument.Instrument()

    ds = DataSet()
    
    try:
        ds.datafiles = 100
        assert False
    except AttributeError:
        assert True
    
    try:
        ds.normalizationfiles = 10
        assert False
    except AttributeError:
        assert True

    try:
        ds.settings={}
        assert False
    except NotImplementedError:
        assert True

    try:
        ds.convertedfiles = 10
        assert False
    except AttributeError:
        assert True

    try:
        ds.calibrationfiles = 10
        assert False
    except AttributeError:
        assert True

    try:
        print(ds.instrument)
        ds.initialize()
        assert False
    except RuntimeError:
        assert True

    ds.instrument=Instr
    print(ds.datafiles)
    print(ds.normalizationfiles)
    try:
        ds.initialize()
        assert False
    except RuntimeError:
        assert True

    ds.datafiles = 'Here.h5'
    try:
        ds.initialize()
        assert False
    except RuntimeError:
        assert True

            

def test_DataSet_Equality():
    D1 = DataSet(datafiles='Here',normalizationfiles=['Van/1.nxs','Van/2.nxs'])
    assert(D1==D1)

def test_DataSet_SaveLoad():
    Instr = Instrument.Instrument()
    D1 = DataSet(datafiles='Here',normalizationfiles = 'Van.nxs',instrument=Instr)

    temp = 'temporary.bin'

    D1.save(temp)
    D2 = DataSet()
    D2.load(temp)
    os.remove(temp)
    assert(D1==D2) 

def test_DataSet_str():
    Instr = Instrument.Instrument()
    D1 = DataSet(datafiles='Here',normalizationfiles = 'Van.nxs',instrument=Instr)
    string = str(D1)
    print(string)

def test_Normalization_tables():
    Instr = Instrument.Instrument(filename='TestData/CAMEA_Full.xml')
    Instr.initialize()

    NF = 'TestData/VanNormalization.h5'


    dataset = DataSet(instrument=Instr,normalizationfiles=NF)

    try:
        dataset.EnergyCalibration(NF,'TestData/',plot=False,tables=[]) # No binning specified 
        assert False
    except AttributeError:
        assert True

    try:
        dataset.EnergyCalibration(NF,'TestData/',plot=False,tables=['Nothing?']) # Wrong binning
        assert False
    except AttributeError:
        assert True


    dataset.EnergyCalibration(NF,'TestData/',plot=True,tables=['Single']) 
    dataset.EnergyCalibration(savelocation='TestData',plot=False,tables=['PrismaticHighDefinition','PrismaticLowDefinition',2]) 
    assert(dataset.calibrationfiles[-1]=='TestData/EnergyNormalization_2.calib')
    


def test_DataSet_Convert_Data():
    Instr = Instrument.Instrument(filename='TestData/CAMEA_Full.xml')
    Instr.initialize()

    NF = 'TestData/VanNormalization.h5'
    DataFiles = 'TestData/VanNormalization.h5'
    dataset = DataSet(instrument=Instr,normalizationfiles=NF,datafiles=DataFiles)

    calibrationfiles = 'TestData/EnergyNormalization_8.calib'
    

    try:
        dataset.ConvertDatafile(datafiles=DataFiles)
        assert False
    except AttributeError: # Cant find normalization table
        assert True

    if not os.path.exists(calibrationfiles):
        dataset.EnergyCalibration(savelocation='TestData/',tables=['PrismaticHighDefinition'])
    else:
         dataset.calibrationfiles.append(calibrationfiles)

    dataset.ConvertDatafile(datafiles=DataFiles,calibrationfile=calibrationfiles)
    dataset.ConvertDatafile(savelocation='TestData')


def test_DataSet_3DMesh():
    
    x = np.linspace(0,1,2)
    y = np.linspace(0,1,5)
    z = np.linspace(1,2,5)

    X,Y,Z = np.meshgrid(x,y,z,indexing='ij')
    XT1,YT1,ZT1 = calculateGrid3D(X,Y,Z)

    assert(XT1.shape==(3,6,6))
    assert(np.all(XT1[:,0,0]==np.array([-0.5,0.5,1.5])))
    assert(np.all(YT1[0,:,0]==np.array([-0.125,0.125,0.375,0.625,0.875,1.125])))
    assert(np.all(YT1[0,:,0]==ZT1[0,0,:]-1.0))



def test_DataSet_BinData():
    I = np.random.randint(0,100,(10,20,30))
    Norm = np.random.rand(10,20,30)
    Posx = np.linspace(0,1,10)
    Posy = np.linspace(0,1,20)
    Posz = np.linspace(1,2,30)
    PosX,PosY,PosZ = np.meshgrid(Posx,Posy,Posz,indexing='ij')



    pos = [PosX.flatten(),PosY.flatten(),PosZ.flatten()]
    Data,bins = binData3D(0.5,0.25,0.25,pos,I,norm=Norm)

    ReBinnedI = Data[0]
    RebinnedNorm = Data[1]
    RebinnedNormCount = Data[2]


    assert(ReBinnedI.shape==(3,5,5))
    assert(np.all(bins[0]==np.linspace(-0.25,1.25,4)))
    assert(np.all(bins[1]==np.linspace(-0.125,1.125,6)))
    assert(np.all(bins[2]==np.linspace(1-0.125,2.125,6)))

    assert(RebinnedNorm.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.shape==ReBinnedI.shape)
    assert(RebinnedNormCount.dtype==int)
    assert(RebinnedNorm.dtype==Norm.dtype)
    assert(ReBinnedI.dtype==I.dtype)


def test_DataSet_full_test():
    import MJOLNIR.Data.Viewer3D
    import warnings
    plt.ioff()
    Instr = Instrument.Instrument(filename='TestData/CAMEA_Full_2.xml')
    Instr.initialize()

    NF = 'TestData/VanNormalization.h5'
    DataFile='TestData/cameasim2018n000005.h5'


    dataset = DataSet(instrument=Instr,normalizationfiles=NF,datafiles=DataFile)
    dataset.EnergyCalibration(tables=[8],savelocation='TestData')
    dataset.ConvertDatafile(savelocation='TestData')

    Data,bins = dataset.binData3D(0.05,0.05,0.2)

    BinnedData = Data[0]

    warnings.simplefilter("ignore")
    Intensity = np.divide(BinnedData[0]*BinnedData[3],BinnedData[1]*BinnedData[2])
    warnings.simplefilter('once')

    viewer = MJOLNIR.Data.Viewer3D.Viewer3D(Intensity,bins)
    del viewer