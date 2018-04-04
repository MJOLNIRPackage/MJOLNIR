Generate normalization table from data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Before real data can be converted from pixel position and counts into S(q,omega), one needs decide the binning used for the data as well as generate the normalization tables used. This is done using a Vanadium scan file containing a suitable number of energy steps. Three different binnings for each energy at all of the detectors are default for the CAMEA backend:

 - 8 pixels ('PrismaticHighDefinition')
 - 3 pixels ('PrismaticLowDefinition')
 - 1 pixel  ('Single')
 - n pixels (n is integer)

Having chosen binning(s) one creates the tables either with or without creating fit plots at the same time. Creating these does indeed increase runtime a lot but is needed when one wants to inspect the fitting performed. A error will be raised if the number of peaks found in the data file does not match the number of analyser the detectors are exposed to. 

+-------------------------------------------------------------------+------------------------------------------------------------------+
| .. literalinclude:: ../../Tutorials/Generate_normalization.py     |  .. RawData:                                                     |
|     :lines: 4-                                                    |                                                                  |
|     :language: python                                             |  .. figure:: ../../TestData/Raw/Fit_wedge_4.png                  |
|     :linenos:                                                     |    :width: 90%                                                   |
|                                                                   |                                                                  |
|                                                                   | Plot of fit to data integrated in the energy direction for wedge |
|                                                                   | 4.                                                               |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
|  .. SimpleInstrumentA4:                                           |  .. SimpleInstrumentEf:                                          |
|                                                                   |                                                                  |
|  .. figure:: ../../TestData/Raw/Active_52.png                     |  .. figure:: ../../TestData/8_pixels/Detector52.png              |
|    :width: 90%                                                    |    :width: 90%                                                   |
|                                                                   |                                                                  |
| Active area of detector 52 as defined by 3 sigmas away from center| Fit of peaks in vanadium data for detector 52 when using a       |
| pixel, where red denotes active and black inactive.               | a binning of 8 pixels per analyser.                              |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+

In the end, it is the data in the normalization file, in the above case denoted EnergyNormalization_8.calib and located in the TestData, folder that is essensial. It contains the normalization and energy location of all peaks on all detectors in the format:
 - Detector
 - Energy
 - Pixel
 - Amplitude
 - Center
 - Width
 - Background
 - lowerBin
 - upperBin

on each line starting with detector 0, analyser 0, pixel 0 increasing index of pixel, then analyser and lastly detector. 

