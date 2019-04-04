=====================
Visualization methods
=====================

This section is dedicated to give both an overview of the available visualization methods and also an 
in depth explanation of the workings of them.


Discussed visualization tools for 'Live view'
.............................................

As discussed in the technical report of F. Groitel and S. Tóth, the is following non-exhaustive list 
features and visualizations needed by the user during an experiment at the CAMEA instrument

 - 1D line plots

    - Intensity as function of :math:`\Delta E` - So-called Dimer plot

    - Constant energy between two :math:`Q`-points

    - Constant :math:`Q`-point against :math:`\Delta E`

 - 2D Colour plots

    - Powder average, intensity against :math:`\Delta E` and :math:`\left|Q\right|`

    - Constant energy planes, intensity against :math:`Q_x` and  :math:`Q_x` or scattering plane 
    vectors for given energy width

    - Intensity against detector index and un-binned pixel

    - Intensity against detector index and :math:`\Delta E`

    - Intensity against energy and :math:`Q`-points, for a list of :math:`Q`-points plot the 
    intensity as function of :math:`\Delta E` and :math:`\vec{Q}`

The last bullet under the 2D plots has been added since the last review of the document (medio February 2017).

Some of the above visualizations have already been met by the developed software, and will 
below be explained in some detail as to elucidate the underlying algorithms and their limitations. 
Software specifications can be found in the 'Data Module'_


For more details, see <DataSet>ShortAnchor_.
