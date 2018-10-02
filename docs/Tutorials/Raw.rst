.. _Raw-plotting-and-fitting:

Raw plotting and fitting
^^^^^^^^^^^^^^^^^^^^^^^^
For specific details on the objevt c.f. :ref:`Viewer1D<Viewer1D>`. 


The visualizer is inteded to take a data-set, plot a single cut from it and perform fitting on it with simple functions (currently the Gaussian and Lorentzian are available). Durin initial initialization this is of great importance as to refine the UB matrix, find the A3 offset or set up the goniometers. As it is only intented to deal with one 1D set at a time, no functionality is currently developed to display more than 1 at a time. However, one can cycle through the provided data both along the different scan parameters and data values.

Below is a table of the shortcuts available in the different states of the program:

+------------------------------+--------+
| All States                   | Key    |
+==============================+========+
| Quit                         | q      |
+------------------------------+--------+
| Copy current values          | ctrl+c |
+------------------------------+--------+
| Cycle up in data             | up     |
+------------------------------+--------+
| Cycle down in data           | down   |
+------------------------------+--------+
| Cycle up in scan parameter   |        |
|                              | right  |
+------------------------------+--------+
| Cycle down in scan parameter | left   |
+------------------------------+--------+

The initial window shown consists of a text part at the top and a plot of the current data below. This initial window allows the following key presses:

+-------------------+-------------+
| Initial State     | Key         |
+===================+=============+
| Initialize fittng | i or ctrl+i |
+-------------------+-------------+

By pressing 'i' or 'ctrl+i' one starts the fitting. Here one can change between different fit functions. These are currently limited to a Gaussian and a Lorentzian function. Having choosen the fitting function one moves the mouse onto the canvas and left clicks corresponding to the guess one have for the parameter(s) shown in bold font (multiple parameters choosen depending on the fit function).


+----------------------------------+----------------+
| Fitting State                    | Key            |
+==================================+================+
| Choose Gaussian function         | 0              |
+----------------------------------+----------------+
| Choose Lorentzian function       | 1              |
+----------------------------------+----------------+
| Execute fit                      | e              |
+----------------------------------+----------------+
| Print current fit to terminal    | i or ctrl+i    |
+----------------------------------+----------------+
| Choose guess for bold parameters | leftmouseclick |
+----------------------------------+----------------+
| Cycle backwards in parameters    | r              |
+----------------------------------+----------------+

When the initial guess is satisfactory one presses 'e' to execute the fit. The x and y data is then fitted with the shown parameters as initial guesses and the errorbars as absolute standard diviations by the scipy.optimize.curve\_fit function. 

.. note::
    Errorbars being zero is reset to 1 and fit is performed with the least squares algorithm.

After a fit execution the found parameters are written in place of the guess parameters together with the square-root of the correlation matrix. (Wrongly) Assuming no correlation between parameters, the fitting variances
 are given by the diagonal elements. By pressing 'ctrl+c' the fitted parameters are copied to the clipboard. If another fit is wanted, one can press 'i' or 'ctrl+i' to initialize another fit.

+------------------------+-----+
| Executed State         | Key |
+========================+=====+
| Initialize another fit | i   |
+------------------------+-----+



Further shortcuts corresponds to the standard keys of `Matplotlib shortcuts`__. 




.. _MatplotlibShortCuts: https://matplotlib.org/users/navigation_toolbar.html#navigation-keyboard-shortcuts

__ MatplotlibShortCuts_
