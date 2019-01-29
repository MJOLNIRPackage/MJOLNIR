.. _Command-Line-Tutorials:

Command Line Tutorials
======================

Despite the full usability and customizability of using a scripting interface to the visualization software, during an experiment or in order to quickly get an overview of data one is more interested in just inspecting the data using a standerdized set of plotting parameters. It could also be that Python is not a language the user masters and thus a command line interface might be easier when simple visualization is needed. 

In any case, the following tutorials seek to introduce and explain the possibilities of using the command line scripts to quickly plot different parts of the data. However, running the scripts from the command line is operation system dependent. That is, if you are running either Linux or Mac, chances are that you can simply run::

    python ScriptFile.py *args

or if python has been installed as an environment variable::

    ./ScriptFile.py *args

where in the latter case the shebang-command in the top of the scripting file is run.
 
.. toctree::
   :maxdepth: 1

   NormalizationTable

