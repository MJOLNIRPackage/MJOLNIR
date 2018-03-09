Introduction
============

This is the introductonary page for the MOLJNIR software package.
The main purpose of this document is to give an overview of different features of the software and how you can contribute to it.

Module documentation
--------------------
Each module is supposed to be independent from the rest of this software suit. That is, it is supposed
to be working on its one without the need of other peices or moduels. However,
some possitive synagy is possible....



Software Structure
------------------

The software is devided into modules with submodules if applicable.

License
-------
The software package is released under the software lincense Mozilla Public License Version 2.0 to allow for redistribution and usage of the code. If this linces file has not been shipped with your distribution of the code, please find it here: licence_.




.. _Licence: https://choosealicense.com/licenses/



Contribtion
------------

After the initial upstart phase you are more than welcome to contribute to the software. This is best done by:

* First create an issue on the GitHub page describing the scope of the contribution
   
   * Title: *Contribution: Title of contribution*.
   * Short describtion of features.
   * List of package dependencies.

* After discussion of feature scope you are welcome to start the programming
* Suggested changes are submitted through a pull request
* Each contribution needs to include:

    * souce code in a suitable subfolder (see Software structure)
    * Documentation of code located in the docs-folder having identical structure to the modules added
    * Suitable tests for the new functionality added to the .travis.yml-file
    * Needed packages added to the reqirements.txt file


Contribution Example:
_____________________

    Title: Extension of fitting rountine

    Describtion: An extension of the fitting module is needed to allow
    users to incorporate Bayesian update of parameters fitted in 3D
    with given priors. This is to be done by adding a subroutine to 
    the fit object.

    Initial thoughts: It is believed that building upon the article [CIT2002]_, 
    where this feature was incorporated in a C++ script.


.. [CIT2002] Citation of relevant article or other document.
