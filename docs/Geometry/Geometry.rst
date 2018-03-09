Geometry Module
===============


The geometry module is created in order to build a virtual copy of the instrument in the MJOLNIR. This is done using the following classes

.. currentmodule:: Geometry

.. autosummary::

   GeometryConcept.GeometryConcept
   GeometryConcept.GeometryObject
   Detector.Detector
   Detector.TubeDetector1D
   Analyser.Analyser
   Analyser.FlatAnalyser
   Wedge.Wedge


Below is an extended description of the different classes and their methods.



.. automodule:: Geometry

Geometry base class
^^^^^^^^^^^^^^^^^^^

General object from which other geometry objects inherits.

.. automodule:: GeometryConcept

.. autoclass:: GeometryConcept
    :members:


-----------------

.. autoclass:: GeometryObject
   :members: 


Detectors
^^^^^^^^^

.. automodule:: Detector

.. autoclass:: Detector
   :members: 

------------------

.. autoclass:: TubeDetector1D
   :members:



Analysers
^^^^^^^^^

.. automodule:: Analyser

.. autoclass:: Analyser
   :members: 


------------------

.. autoclass:: FlatAnalyser
   :members:


Wedge
^^^^^

.. automodule:: Wedge

.. autoclass:: Wedge
    :members:
