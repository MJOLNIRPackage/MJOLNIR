After the initial upstart phase you are more than welcome to contribute to the software. This is best done by:

* First create an issue on the GitHub page describing the scope of the contribution
   
   * Title: *Contribution: Title of contribution*.
   * Short description of features.
   * List of package dependencies.

* After discussion of feature scope you are welcome to start the programming
* Suggested changes are submitted through a pull request
* Each contribution needs to include:

    * source code in a suitable sub-folder (see Software structure)
    * Documentation of code located in the docs-folder having identical structure to the modules added
    * Suitable tests for the new functionality added to the .travis.yml-file
    * Needed packages added to the requirements.txt file


Contribution Example:
_____________________

    Title: Extension of fitting routine

    Description: An extension of the fitting module is needed to allow
    users to incorporate Bayesian update of parameters fitted in 3D
    with given priors. This is to be done by adding a subroutine to 
    the fit object.

    Initial thoughts: It is believed that building upon the XXX package 
    where this feature was already created.
