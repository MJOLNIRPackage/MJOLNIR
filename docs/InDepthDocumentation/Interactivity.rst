.. :Interactivity:

Interactivity
^^^^^^^^^^^^^

Working with the 3D data from CAMEA can be a challenge. Thus, a set of tools have been created...


For the QPlane, QE axis, and the QELine axis a set of interactive modes have been created. These include the resolution and cutting modes

.. table:: Overview of possible interactive modes for Axis in MJOLNIR
    :align: center

    +------------------+-----------+------------+
    | Interaction Mode | Activate  | Deactivate |
    +==================+===========+============+
    | Resolution       | r         | b          |
    +------------------+-----------+------------+
    | Cutting          | c         | b          |
    +------------------+-----------+------------+
    | Manage Mode [#]_ | m         | m/n/c/b    |
    +------------------+-----------+------------+
    | View 3D [#]_     | N/A       | N/A        |
    +------------------+-----------+------------+


.. [#] Management Mode is only available from the cutting mode as it is connected to clearing/moving cuts.

.. [#] View3D is activated by creating a view3D through the data set. For more details see `Interactivity <../Tutorials/Quick/Viewer3D.html>`_.



Resolution
----------

With the resolution mode activate, the resolution ellipsoid is calculated for the current mouse position on the data. It utilizes an updated version of the method presented in 
[eck14] G. Eckold and O. Sobolev, NIM A 752, pp. 54-64 (2014) and has been implemented `Tobias Weber <tweber@ill.fr>`_ in the software `TARKIN <https://doi.org/10.5281/zenodo.4117437>`_.

The ellipsoid is inherently 4D (dimensions are Qx, Qy, Qz, and energy), but the plots upon which the resolution is shown is 2D. By integration over the two coordinates not 
used in the figure, one can find the resolution matrix in the relevant plane. As this matrix (in this approximation) only dependent on the instrument, the sample rotation is to be taken into
account as well as the possible rescaling of the plotting axis. 

When the left mouse button is pressed on a specific point, the resolution matrix together with the in-plane eigenvectors and widths are printed and the ellipse(s) remain
on the plot. A resolution is calculated for each instrumental setting corresponding to the point of interest. If two different incoming energies result in an overlap at 
the point, two ellipses are plotted and information is printed for both.

.. table:: Interaction possibilities for the resolution mode
    :align: center

    +----------------------+------------+
    | Method               | Button     |
    +======================+============+
    | Calcualte Resolution | Left Click |
    +----------------------+------------+
    | Clear Ellipses       | c          |
    +----------------------+------------+
    | Return               | b          |
    +----------------------+------------+


Cutting
-------

On an axis, one can start the cutting mode by pressing 'c' and exit it again by pressing 'b'. In this mode 1D cuts can be created, managed, and plotted. Depending on the axis 
in question a number of cuts are possible, including a rectangular cut between two reciprocal space points at  constant energy, or an energy cut at a specific Q point. 
By pressing 'n' one cycles through the possible cuts with the type of cut printed in the coordinate format text (i.e. the text showing the current 
mouse position) and after all cutting types have been cycled through you will end up returning to the initial inactive starting position. To initiate the cut click with
the left mouse button. This point signifies, depending on the cut type, the start or the centre of the cut. After providing all positions necessary, the cut will be shown in red
and a 1D cut can be created by pressing 'c' and in some circumstances a 2D plot by 'Ctrl + c'. If a second cut is created the first will turn green to signify that it is now inactive. To change which cut is active and to move or delete
cuts the manage mode has been created. Activating it by pressing 'm' one can activate a cut by left clicking on it or delete it by right-clicking it. By click-and-drag the cut can be translated.
Cuts cannot be rescaled or rotated.

.. table:: Interaction possibilities for the cutting mode
    :align: center

    +-------------------+------------+
    | Method            | Button     |
    +===================+============+
    | New/Cycle Cut Type| n          |
    +-------------------+------------+
    | Define Cut        | Left Click |
    +-------------------+------------+
    | Enter Manage Mode | m          |
    +-------------------+------------+
    | Perform 1D Cut    | c          |
    +-------------------+------------+
    | Perform 2D Cut*   | Ctrl + c   |
    +-------------------+------------+
    | Snapping Cycle    | **.**      |
    +-------------------+------------+
    | Return            | b          |
    +-------------------+------------+

* Only in some circumstances, 2D plots are possible.

When in cutting mode, it might be advantageous to be able to snap to the underlying grid of the figure. A snapping option has been created where one can snap to the grid or a specific 
segment size, 1/n, for n = 1, 2, 3, 4, and 5. This mode is activated by pressing '**.**' while in cutting mode. The cycle goes as follows: None, 1, 2, 3, 4, 5 and then back to None.