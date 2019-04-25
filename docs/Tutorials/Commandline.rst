.. _Command-Line-Tutorials:

Command Line Tutorials
======================

Despite the full usability and customizability of using a scripting interface to the visualization software, during an experiment or in order to quickly get an overview of data one is more interested in just inspecting the data using a standardized set of plotting parameters. It could also be that Python is not a language the user masters and thus a command line interface might be easier when simple visualization is needed. 

In any case, the following tutorials seek to introduce and explain the possibilities of using the command line scripts to quickly plot different parts of the data. However, running the scripts from the command line is operation system dependent. That is, if you are running either Linux or Mac, chances are that you can simply run:

.. code-block:: bash

    MJOLNIRScriptFile *args

where the shebang-command in the top of the scripting file is run on linux or Mac. 

When installing the MJOLNIR package through PiP it is sought for that the relevant script files are place in the correct position for the OS. All of the scripts 
starts with the prefix MJOLNIR in order to avoid any conflict with already available command line calls.

Below is an overview of currently implemented scripts and their arguments.
 
.. toctree::
   :maxdepth: 1

   Commandline/MJOLNIRCalibrationInspector
   Commandline/MJOLNIRConvert
   Commandline/MJOLNIR3DView