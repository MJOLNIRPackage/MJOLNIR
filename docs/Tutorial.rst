Tutorials
=========

Contains a list of different tutorials in order to familiarize users with the objects and methods used in MJOLNIR.


Build a simple instrument
^^^^^^^^^^^^^^^^^^^^^^^^^
In order to build a complete instrument using the MJOLNIR :ref:`GeometryModule` module, a lot of different objects need to come together. This tutorial sets out to introduce these step by step with each example increases complexity. The following is an example of how to build a simple instrument consisting of one :ref:`Detectors<TubeDetector1D>` and one :ref:`Analysers<FlatAnalyser>` grouped togehter by a :ref:`Wedge<Wedge>` object:


+---------------------------------------------------------------+--------------------------------------------------------+
| .. literalinclude:: ../Tutorials/Build_a_simple_instrument.py | .. _Build_a_simple_instrument_fig:                     |
|     :lines: 4-27,29-                                          |                                                        |
|     :language: python                                         | .. figure:: ../Tutorials/Build_a_simple_instrument.png |
|     :linenos:                                                 |   :width: 75%                                          |
|                                                               |                                                        |
|                                                               |   The figure produced by the current code example.     |
+---------------------------------------------------------------+--------------------------------------------------------+


Calcualte A4 and Ef
^^^^^^^^^^^^^^^^^^^
This is an example of initialization of an instrument with two detectors and two analysers. The detectors each have 10 pixels, where the first 5 are looking at the first analyser and the last 5 are looking at the second analyser. This gives rise to different energies. As the second detector is moved away from the straight line through the analyser its A4 values also changes.

+---------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+
| .. literalinclude:: ../Tutorials/Calcualte_A4_and_Ef.py       |[[array([ 1.57079633,  1.57079633,  1.57079633,  1.57079633,  1.57079633,                                                                |
|     :lines: 4-                                                |1.57079633,  1.57079633,  1.57079633,  1.57079633,  1.57079633]), array([ 1.55420013,  1.55430836,  1.55451765,  1.55481497,  1.55518368,|
|     :language: python                                         |1.52099283,  1.52123672,  1.5217106 ,  1.52238902,  1.52323891])]]                                                                       |
|     :linenos:                                                 |                                                                                                                                         |
|                                                               |[[array([ 3.83622353,  4.27947022,  4.81164763,  5.44262522,  6.18119629,                                                                |
|                                                               |3.83622353,  4.27947022,  4.81164763,  5.44262522,  6.18119629]), array([ 3.8414743 ,  4.28381055,  4.81454296,  5.44345335,  6.17926234,|
|                                                               |3.84580202,  4.29004983,  4.82331398,  5.45545541,  6.19525547])]]                                                                       |
+---------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------+



Load instrument from XML file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to not having to create an instrument from scratch each time a data treament is performed, one way is to load the instrument from a XML file formated as below. What is important is the structure of the XML file, where the outer object is the instrument with all of it settings; middle part is the wedge(s) and inner part all of the detectors and analysers. All objects in the instrument has their attributes defined in the opening bracket of the XML object and nothing between it and the closing bracket.

+-------------------------------------------------------------------+------------------------------------------------------------------+
|      .. literalinclude:: ../Tutorials/SimpleInstrument.xml                                                                           |
|         :language: xml                                                                                                               |
|         :linenos:                                                                                                                    |
|                                                                                                                                      |
|                                                                                                                                      |
|                                                                                                                                      |
+-------------------------------------------------------------------+------------------------------------------------------------------+
| .. literalinclude:: ../Tutorials/Load_Instrument_From_XML_file.py |  .. SimpleInstrument:                                            |
|     :lines: 4-21,23-33,35-44,46                                   |                                                                  |
|     :language: python                                             |  .. figure:: ../Tutorials/SimpleInstrument.png                   |
|     :linenos:                                                     |    :width: 90%                                                   |
|                                                                   |                                                                  |
|                                                                   | Plot of instrument loaded from the XML file. Remember, that the  |
|                                                                   | sample is located at the origin (0,0,0)                          |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
|  .. SimpleInstrumentA4:                                           |  .. SimpleInstrumentEf:                                          |
|                                                                   |                                                                  |
|  .. figure:: ../Tutorials/SimpleInstrument_A4.png                 |  .. figure:: ../Tutorials/SimpleInstrument_Ef.png                |
|    :width: 90%                                                    |    :width: 90%                                                   |
|                                                                   |                                                                  |
| Scattering angles of the different dectors, where individual      | Energies of the different dectors, where individual pixels hit   |
| pixels hit  different detectors.                                  | different detectors.                                             |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+




Generate normalization table from data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Before real data can be converted from pixel position and counts into S(q,omega), one needs decide the binning used for the data as well as generate the normalization tables used. This is done using a Vanadium scan file containing a suitable number of energy steps. Three different binnings for each energy at all of the detectors are default for the CAMEA backend:

 - 8 pixels ('PrismaticHighDefinition')
 - 3 pixels ('PrismaticLowDefinition')
 - 1 pixel  ('Single')
 - n pixels (n is integer)

Having chosen binning(s) one creates the tables either with or without creating fit plots at the same time. Creating these does indeed increase runtime a lot but is needed when one wants to inspect the fitting performed. A error will be raised if the number of peaks found in the data file does not match the number of analyser the detectors are exposed to. 

+-------------------------------------------------------------------+------------------------------------------------------------------+
| .. literalinclude:: ../Tutorials/Generate_normalization.py        |  .. RawData:                                                     |
|     :lines: 4-                                                    |                                                                  |
|     :language: python                                             |  .. figure:: ../TestData/Raw/Fit_wedge_4.png                     |
|     :linenos:                                                     |    :width: 90%                                                   |
|                                                                   |                                                                  |
|                                                                   | Plot of fit to data integrated in the energy direction for wedge |
|                                                                   | 4.                                                               |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
|  .. SimpleInstrumentA4:                                           |  .. SimpleInstrumentEf:                                          |
|                                                                   |                                                                  |
|  .. figure:: ../TestData/Raw/Active_52.png                        |  .. figure:: ../TestData/8_pixels/Detector52.png                 |
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

Convert data to Q and Omega
^^^^^^^^^^^^^^^^^^^^^^^^^^^
With the above normalization table created, on can easily convert raw data files using the method in the :ref:`DataSet<DataSet>` called ConvertDatafile as 

.. literalinclude:: ../Tutorials/Convert_Data.py
    :lines: 4-
    :language: python 
    :linenos:  

The code then converts the scan file, be it either Ei, A4, or A3 scan, into a new HDF file following the Nexus NXsqom data convention. This is followed in order to facilitate easy interfacing with other software later used. The location of the converted file is the same as the original but the new file has the ending .nxs. Furthermore, in order to store all information and ensure that no data is lost, all information already present in the raw data file is copied into the new. This also include the original raw data counts.


Bin data and visualize
^^^^^^^^^^^^^^^^^^^^^^
Having converted the data into the Nexus NXsqom format, one wants to both rebin the data and visualize it. As different detector pixels covers different positions in reciprocal space, one needs to rebin the data in order to avoid large areas of no data. This is done using the method in the :ref:`DataSet<DataSet>` called binData3D. As input one needs to provide the step size in the x, y, and z directions, the 3D position and intensity. Furthermore, normalization and monitor count can be specified in order to also bin these. Returned is the rebinned data together with normalization and monitor count, if applicable, and the bins used.

+-------------------------------------------------------------------+------------------------------------------------------------------+
|      .. literalinclude:: ../Tutorials/Binning_data.py             |    .. Visualization_E_546:                                       |
|         :lines: 4-22,24,27,30-                                    |                                                                  |
|         :language: python                                         |    .. figure:: ../Tutorials/Visualization_E_546.png              |
|         :linenos:                                                 |      :width: 90%                                                 |
|                                                                   |                                                                  |
| Binning of the converted data into Qx and Qy bins of size 0.05AA  |  Cut through data along the energy direction showing Qx and Qy   |
| and energy in 0.2 meV. Intensity is calculated and with the bins  |  for a phonon scan at the energy 5.46 meV.                       |
| passed to the visualizer.                                         |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+
|  .. Visualization_Qx_119:                                         |  .. SimpleInstrumentEf:                                          |
|                                                                   |                                                                  |
|  .. figure:: ../Tutorials/Visualization_Qx_119.png                |  .. figure:: ../Tutorials/Visualization_Qy_089.png               |
|    :width: 90%                                                    |    :width: 90%                                                   |
|                                                                   |                                                                  |
| Cut of data along the Qx direction.                               | Cut of data along the Qx direction.                              |
|                                                                   |                                                                  |
|                                                                   |                                                                  |
+-------------------------------------------------------------------+------------------------------------------------------------------+

The bins and calculated intensity is then passed on to the Viewer3D object, that generates a matplotlib figure. This plot is made interactive through the slider in the bottom, that shows the current position along the axis as well as the value of this. By pressing the up/down arrows (or +/- buttons as well as scrolling) one can change this value and thus investigate the third dimension. By default the energy direction is chosen but it can be changed by pressing the 0, 1, or 2 numerical buttons. This is to be understood as which direction is to be sliced. Further, by clicking on the plot, the x and y value corresponding to the point is printed to the terminal. This is then intented to be used in further data treatment when cuts are needed. 

