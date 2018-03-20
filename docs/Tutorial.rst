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
|     :lines: 4-21,23-34,36-45,47                                   |                                                                  |
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









