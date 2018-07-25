Introduction
============

This is the introductonary page for the MOLJNIR software package. The main purpose of this document is to give an overview of different features of the software and how you can contribute to it.

The software is currently developed solely by Jakob Lass, PhD student at both the Niels Bohr Institute, Copenhagen - Denmark, and the Paul Scherrer Institute, Villigen - Switzerland, and is not developed by a professional team. The software is intended to be used for data treatment and visualization for the CAMEA upgrad at the RITA II instrument at SINQ, PSI Villigen - Switzerland. 

The software is found at GitHub_ and is intended to be used together with Python versions 2.7, (3.4,) 3.5, and 3.6. This compability is ensured by the use of automated unit test through the Travis project (Travis_). Python 3.4 is no longer tested due to updates in the Travis testing frameworl. Further than just testing, as to ensure a thorough testing the coverage of these are monitored using Coverals (Coveralls_). However, certain algorithms and methods are not suited to be tested through simple tests. This includes graphical methods where one for example uses a plotting routine to generate a specific output. Though the visual inspection is far outside of the testing scope for this software, some of the methods are still tested by simple run through test. That is, if they can be run and generate a plot without crashing and throwing an error, it is believed that they work as intended. This is where acutal user testing is needed. 

.. Module documentation
.. ^^^^^^^^^^^^^^^^^^^^
.. Each module is supposed to be independent from the rest of this software suit. That is, it is supposed
.. to be working on its one without the need of other peices or moduels. However,
.. some possitive synagy is possible....


Software Structure
^^^^^^^^^^^^^^^^^^

The software is devided into individual modules being Instrument, DataSet, and Statistics. With this division it is intended that each part of the software suit is to be fully independent of the others but may be used together. The same goes for the tutorials that are intended to cover all of the methods and workflows a user would come into contact with while using the software.

Installation
^^^^^^^^^^^^
The inteded way for the software to be installed on a computer is currently to navigate to the github page GitHub_ and download the latest release, usually found in the main branch. For an up to date version with new features, the 1.0.0 branch can be used, but is not to be assumed complete. It is recommended to create a new Anaconda environment or similar for the installation, but this is not required. By utilizing the pip package manager:

.. code-block:: bash

   cd MJOLNIR
   pip install .


Alternatively one can install it directly by running 

.. code-block:: bash

   cd MJOLNIR
   python install.py install


Both of the above methods installs the software to the specified environment. If only a temporary local build is wanted, the latter aproach can be used with 'install' replaced  by 'build.



License
^^^^^^^
The software package is released under the software lincense Mozilla Public License Version 2.0 to allow for redistribution and usage of the code. If this linces file has not been shipped with your distribution of the code, please find it here: licence_.




.. _Licence: https://choosealicense.com/licenses/

.. _GitHub: https://github.com/Jakob-Lass/MJOLNIR/

.. _Coveralls: https://coveralls.io/github/Jakob-Lass/MJOLNIR/

.. _Travis: https://travis-ci.org/Jakob-Lass/MJOLNIR/

.. Contribution
.. ^^^^^^^^^^^^
.. include Contribution/Contribution.rst


Bug Report
^^^^^^^^^^
If an un error or unexpected behaviour of the software is observed, of if a feature is needed, you are more than welcome to create an issue or feature request at the GitHub page (Issues_). Dealing and fixing the reported bugs will be most easily done if both operation system, software version, a minimal working example, and other relevant informations are provided. Further as time goes by, it is hoped that this page will also contain explanations and answers to the most frequently asked question of the software. 

Currently there are the following open issues |open| and closed |closed|

.. |open| image:: https://img.shields.io/github/issues/Jakob-Lass/MJOLNIR.svg?style=plastic 
    :alt: Open GitHub issues
    :target: https://github.com/Jakob-Lass/MJOLNIR/

.. |closed| image:: https://img.shields.io/github/issues-closed/Jakob-Lass/MJOLNIR.svg?style=plastic
   :alt: Closed GitHub issues
   :target: https://github.com/Jakob-Lass/MJOLNIR/




.. _Issues: https://github.com/Jakob-Lass/MJOLNIR/issues
