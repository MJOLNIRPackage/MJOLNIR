.. _GeometryModule:

Geometry Module
===============


The geometry module is created in order to build a virtual copy of the instrument in the MJOLNIR. This is done using the following classes

.. currentmodule:: Geometry



.. autosummary::
   :nosignatures:


    
   GeometryConcept.GeometryConcept
   GeometryConcept.GeometryObject
   Detector.Detector
   Detector.TubeDetector1D
   Analyser.Analyser
   Analyser.FlatAnalyser
   Wedge.Wedge
   Instrument.Instrument

Below is an extended description of the different classes and their methods.



.. automodule:: Geometry

Geometry base classes
^^^^^^^^^^^^^^^^^^^^^

General object from which other geometry objects inherits.


.. automodule:: GeometryConcept

.. _GeometryConcept:

.. autoclass:: GeometryConcept
    :members:


-----------------

.. _GeometryObject:

.. autoclass:: GeometryObject
   :members: 



Detectors
^^^^^^^^^

.. automodule:: Detector

.. _Detector:

.. autoclass:: Detector
   :members: 

------------------

.. _TubeDetector1D:

.. autoclass:: TubeDetector1D
   :members:



Analysers
^^^^^^^^^

.. automodule:: Analyser

.. _Analyser:

.. autoclass:: Analyser
   :members: 


------------------

.. _FlatAnalyser:

.. autoclass:: FlatAnalyser
   :members:


Wedge
^^^^^

.. automodule:: Wedge

.. _Wedge:

.. autoclass:: Wedge
    :members:


Instrument
^^^^^^^^^^

.. automodule:: Instrument

.. _Instrument:

.. autoclass:: Instrument
    :members:
