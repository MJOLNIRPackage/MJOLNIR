.. MJOLNIR documentation master file, created by
   sphinx-quickstart on Wed Feb 28 17:31:13 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MJOLNIR's documentation!
===================================
Welcome to the documentation page of the MJOLNIR software package developed for data treatment, conversion, and visualization for multiplexing neutron scattering 
spectrometers. With the increased point density measured during a single experiment, flexible and versatile software is needed to deal with this new data. This 
software package tries to both introduce an easy-to-use scripting interface to visualization and treatment while still offering customizable plotting features.

Table of contents
-----------------
Below is an overview of the provided documentation containing everything from physical/mathematical explanations of different concets through code specific references to 
tutorials and guides on how to make use of MJOLNIR.

.. toctree::
   :maxdepth: 1

   Introduction
   faq
   Tutorials/Tutorial
   MJOLNIR
   InDepthDocumentation/Documentation
   Optimizations/Optimization


.. toctree::
   :maxdepth: 1

   Commissioning/Commissioning


| 

As the software is to be used by a broad user group, it is of great importance to support a wide range of python versions. Below is the current list of supported versions.

|travis| |coverall| |python27| |python34| * |python35| |python36| |python37| *

.. note::
    Python 3.4 is believed to be compatible but is not tested due to updates in testing framework.
    Python 3.7 is supported by the package but some of the dependencies might not be supported on all operating systems.

.. |travis| image:: https://travis-ci.org/Jakob-Lass/MJOLNIR.svg?branch=1.1.0'
    :target: https://travis-ci.org/Jakob-Lass/MJOLNIR

.. |coverall| image:: https://coveralls.io/repos/github/Jakob-Lass/MJOLNIR/badge.svg?branch=1.1.0'
    :target: https://coveralls.io/github/Jakob-Lass/MJOLNIR?branch=1.1.0'

.. |python27| image:: https://img.shields.io/badge/python-2.7-brightgreen.svg 
    :target: https://travis-ci.org/Jakob-Lass/MJOLNIR
 
.. |python34| image:: https://img.shields.io/badge/python-3.4-yellow.svg
    :target: https://travis-ci.org/Jakob-Lass/MJOLNIR   

.. |python35| image:: https://img.shields.io/badge/python-3.5-brightgreen.svg
    :target: https://travis-ci.org/Jakob-Lass/MJOLNIR 

.. |python36| image:: https://img.shields.io/badge/python-3.6-brightgreen.svg
    :target: https://travis-ci.org/Jakob-Lass/MJOLNIR

.. |python37| image:: https://img.shields.io/badge/python-3.7-green.svg
    :target: https://travis-ci.org/Jakob-Lass/MJOLNIR

