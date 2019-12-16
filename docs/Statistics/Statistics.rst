Statistics Module
=================

This is the documentation part intended for the explanation of the statistics module. Due to the philosophy of this software suite where data visualization and actual data fitting and feature extraction are decoupled, this section is dedicated to explanation of the needed methods and workflows. 
The intended methods could include:

- Fit using either Gaussian or Poisson statistics for:
    - 1D extraction of data from two data points for constant energy
    - 1D extraction of data from 1 Q position and all energies
    - other..

- Advanced fitting routines for multidimensional data (2, 3, or n(?) dimensions) using:
    - Modified costfunctions to ensure continuity and breaks allowed for elastic peaks
    - Sequential fitting with blurring and rebinning of data to allow for initial fitting and guesses later used on the full data set
    - Constant signal to noise binning using voronoi tesselation to perform importance sampling

- Using swarm or other algorithms with all data fitting tool



.. automodule:: Statistics     

