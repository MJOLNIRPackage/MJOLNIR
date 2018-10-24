Load instrument from XML file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to not having to create an instrument from scratch each time a data treament is performed, one way is to load the instrument from a XML file formated as below. What is important is the structure of the XML file, where the outer object is the instrument with all of it settings; middle part is the wedge(s) and inner part all of the detectors and analysers. All objects in the instrument has their attributes defined in the opening bracket of the XML object and nothing between it and the closing bracket.

+----------------------------------------------------------------------+------------------------------------------------------------------+
|      .. literalinclude:: ../../Tutorials/SimpleInstrument.xml                                                                           |
|         :language: xml                                                                                                                  |
|         :linenos:                                                                                                                       |
|                                                                                                                                         |
+----------------------------------------------------------------------+------------------------------------------------------------------+
| .. literalinclude:: ../../Tutorials/Load_Instrument_From_XML_file.py |  .. SimpleInstrument:                                            |
|     :lines: 4-                                                       |                                                                  |
|     :language: python                                                |  .. figure:: ../../Tutorials/SimpleInstrument.png                |
|     :linenos:                                                        |    :width: 90%                                                   |
|                                                                      |                                                                  |
|                                                                      | Plot of instrument loaded from the XML file. Remember, that the  |
|                                                                      | sample is located at the origin (0,0,0)                          |
|                                                                      |                                                                  |
+----------------------------------------------------------------------+------------------------------------------------------------------+
|  .. SimpleInstrumentA4:                                              |  .. SimpleInstrumentEf:                                          |
|                                                                      |                                                                  |
|  .. figure:: ../../Tutorials/SimpleInstrument_A4.png                 |  .. figure:: ../../Tutorials/SimpleInstrument_Ef.png             |
|    :width: 90%                                                       |    :width: 90%                                                   |
|                                                                      |                                                                  |
| Scattering angles of the different dectors, where individual         | Energies of the different dectors, where individual pixels hit   |
| pixels hit  different detectors.                                     | different detectors.                                             |
|                                                                      |                                                                  |
+----------------------------------------------------------------------+------------------------------------------------------------------+





