The Reciprocal Lattice Units Axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This axis is created to easily plot data in RLU for all crystals in the scattering plane. This can at times prove difficult when dealing with non-orthogonal unit cells and/or if the scattering plane is non-trivial. All of this is taken care of when using this axis to plot. For detailed technical information about the implementation see section below or the source code.

.. code-block:: python
   :linenos:

   from MJOLNIR.Data import DataSet
   from MJOLNIR import _tools # Usefull tools useful across MJOLNIR 
   import numpy as np
   import matplotlib.pyplot as plt
   
   
   numbers = '137' # String of data numbers
   fileList = _tools.fileListGenerator(numbers,'/Path/To/Data/',2018) # Create file list from 2018 in specified folder
   
   ds = DataSet.DataSet(fileList)
   ds.convertDataFile(saveFile=False)
   mask = np.zeros_like(ds.I.data) # Define mask, see FAQ for explanation
   mask[:,:,:3]=True
   ds.mask = mask
   
   # Create RLU axis
   ax = ds.createRLUAxes()
   # Get the figure corresponding to the returned axis
   fig = ax.get_figure()
   
   # The axis should contain v1 and v2
   v1 = [-1.5,-1.5,0]
   v2 = [1.5,1.5,0]
   
   ax.set_axis(v1,v2)
   
   # Generate some points to plot
   points = np.array([[-1,-1],[-1,1],[1,1],[1,-1],[-1,-1]]).T
   
   # Use points as distance along projections
   projection = np.array(ax.sample.inv_tr(points[0],points[1]))
   # Convert into HKL (corresponds to the points above)
   HKL = np.array([points[0],points[1],np.zeros_like(points[0])])
   
   # Calculate actual position of HKL points
   P1P2 = np.array(ax.sample.calculateHKLtoProjection(HKL[0],HKL[1],HKL[2]))
   QxQy = np.array(ax.sample.inv_tr(P1P2[0],P1P2[1]))
   
   # Plot all 3 lines
   ax.plot(points[0],points[1],label='QxQy')
   ax.plot(projection[0],projection[1],label='Projection')
   ax.plot(QxQy[0],QxQy[1],'g',label='HK',marker='o', linestyle='dashed')
   
   
   ax.legend()
   
   ax.grid(True)
   
   fig.savefig('figure0.png',format='png')
   
   ############### Second example ###############
   numbers = '422' # String of data numbers
   fileList = _tools.fileListGenerator(numbers,'/Path/To/Data/',2018) # Create file list from 2018 in specified folder
   
   ds = DataSet.DataSet(fileList)
   ds.convertDataFile(saveFile=False)
   
   ax = ds.createRLUAxes()
   fig = ax.get_figure()
   
   # The axis should contain v1 and v2
   
   v1 = [-2,1,1]
   v2 = [2,-1,-1]
   
   ax.set_axis(v1,v2)
   
   ax.grid(True)
   
   fig.savefig('figure1.png',format='png')
   

The first example in the code above takes a data set measured on a crystal with unit cell [6.11  6.11  11.35  90.    90.   120.], that is, the crystal is hexagonal and is placed in the hk0 scattering plane. Behind the scene, data plotted in the axis is assumed to be in the :math:`Q_x`, :math:`Q_y` coordinate system rotated such that the first projection vector is along the x-axis (in this case (h,0,0). The data is then first converted into reciprocal space (HKL) and then into the projected in the scattering plane. This might seem like swatting flies with a sledge hammer, however, this allows for enough abstraction as to allow all cases to be plotted.

As the axis object makes use of the non-standard "GridHelperCurveLinear" many features regarding the tick marks and setting limits are not working, but see section below for dicussion. A custom function has been generated to set the limits of the shown graph called "set_limits". It takes at least two vectors in either HKL or projection along principal directions and assures that these are visible. Multiple vectors can be provided and it is then assured that all vector positions are visible. An example is shown in line 22 where it is assured that [-1,-1.5,0] and [2.0,1.0,0] are both visible. In this case, providing the vectors [-1,-1.5] and [2,1] is sufficient as the projections are simply [1,0,0] and [0,1,0]. 

.. figure:: RLUAxis.png
  :width: 50%
  :align: center



As to the data plotted in the axis. Plotting directly into the axis corresponds to plotting the :math:`Q_x` :math:`Q_y` system which in this case produces a square box around (0,0). However, wanting to plot things in terms of reciprocal lattice units one has two options: Plot corresponding to projection or calculate :math:`Q_x` :math:`Q_y` from HKL points. Both of these methods are shown in the first example resulting in the two parallellograms plotted on top of each other. Notice that the blue box does not have the same height as the two others due to the length of the H00 being 1.187 /AA and not unity. 

.. figure:: RLUAxis2.png
  :width: 50%
  :align: center



The second example has a crystal with unit cell [9.843  9.843  9.843 90.    90.    90.], i.e. simple cubic, but the scattering plane is (hkk). Setting the limits to include all points from -2 to 2 in H and from -1 to 1 in K is shown in the code by providing the position vectors. The same result is obtained by simply giving [-2,-1] and [2,1] to the method.

Technical details
-----------------

It was chosen to make use of the "GridHelperCurveLinear" despite the difficulties arising from it as this allows for plotting in RLU coordinates without having to skew data. That is, it is possible to keep data as measured by the instrument (sort of circular) while still providing all information about the reciprocal space to the user. As mentioned, some calculations happen behind the scenes when dealing with this object; most of them in the sample object itself. Mathematically what happens is as follows:

The general relationship between measured points from the instrument, denoted :math:`(Q_x,Q_y,Q_z)` for the two in-plane components along x and y, and the one out of plane, and the reciprocal lattice units HKL is given by the UB matrix

.. math::

    \begin{pmatrix}Q_x\\Q_y\\Q_z\end{pmatrix} = UB \cdot \begin{pmatrix}H\\K\\L\end{pmatrix}

From the geometrical constraints of the CAMEA backend, all scattering is performed in plane. In other words, the :math:`Q_z` is always 0. Thus one can use a simple projection matrix from 2D to 3D:

.. math::

    \begin{pmatrix}Q_x\\Q_y\end{pmatrix} = \underbrace{\begin{pmatrix}1 & 0 & 0\\0 & 1 & 0\end{pmatrix}}_{P_{23}} \cdot UB \cdot \begin{pmatrix}H\\K\\L\end{pmatrix}

As the measured space from the instrument side is 2D the points in terms of RLU are also to lay in a plane. One can make use of this by finding the two most simple projection vectors spanning this plane. Assuming they are found the projection along these are denoted  :math:`P_0` and :math:`P_1`. One can then  project the H, K, L points long these vectors as:

.. math::

    \begin{pmatrix}H\\K\\L\end{pmatrix} = P_M \cdot \begin{pmatrix}P_0\\P_1\end{pmatrix}

Now remains finding the projection matrix :math:`P_M`, which is given as the 3x2 column matrix of them divided by the square of their lengths. Putting it all together results in:

.. math::

   \begin{pmatrix}Q_x\\Q_y\end{pmatrix}=P_{23}\cdot UB \cdot P_M \cdot \begin{pmatrix}P_0\\P_1\end{pmatrix}

One further detail is that due to the way that the instrument positions, :math:`Q_x` and :math:`Q_y` are calculated, one needs to rotate this system with an angle corresponding to the "mis-alignment" of the orientation of the crystal. In reality this angle corresponds to the difference between the A3 zero offset and the one for which the first projection vector is along the x-axis.

In the code, the matrix coupling :math:`Q_x` and :math:`Q_y` to :math:`P_0` and :math:`P_1` is called convert, while the inverse is denoted convertinv. To use these two matrices in the plotting through the "GridHelperCurveLinear" axis, two functions are defined "tr" and "inv_tr" taking projection values to :math:`Q_x` and :math:`Q_y` and reverse respectively. These are provided to the axis and does the calculation for the plotting. For enabling the hover-over tool tip a function calculating :math:`Q_x` and :math:`Q_y` into HKL is created, being simply the matrix multiplication with inverse UB matrix. 

Examples
--------

Using the crystal of example 1 above, the unit cell parameters are  6.11   6.11  11.35  90.    90.   120. resulting in the reciprocal lattice vector matrix:

.. math::

    \mathrm{RLU\ matrix} = \begin{pmatrix}1.187 & 0.594 & -0.000\\0.000 & -1.028 & -0.000\\-0.000 & -0.000 & -0.554\end{pmatrix}

For the specific experiment, the UB matrix is found to be:

.. math::

    UB = \begin{pmatrix}1.097 & 0.155 & 0.000\\-0.454 & -1.177 & -0.000\\0.000 & 0.000 & -0.554\end{pmatrix}

The two projection vectors for the scattering plane HK0 is simply :math:`(1,0,0)^T` and :math:`(0,1,0)^T`, resulting in the convert matrix:

.. math::

    \mathrm{convert} = \begin{pmatrix}1.187 & 0.594\\0.000 & -1.028\end{pmatrix}

As seen, the transformation is non-orthogonal and thus results in the axis shown above. For the "mis-alignment" the rotation angle to correct is found to be -22.5 :math:`^{\mathrm{o}}`. That is, all of the data is rotated by -22.5 degrees before being plotted in the RLU axis.

For the second example shown above, with the cartesian unit cell but the non-trivial scattering plane, the matrices are:

.. math::

    UB &= \begin{pmatrix}-0.418 & 0.341 & 0.341\\0.482 & 0.296 & 0.296\\0.000 & 0.451 & -0.451\end{pmatrix}\\
    \mathrm{RLU\ matrix} &= \begin{pmatrix}0.638 & 0.000 & -0.000\\0.000 & -0.638 & -0.000\\-0.000 & -0.000 & -0.638\end{pmatrix}\\
    \mathrm{convert} &= \begin{pmatrix}0.638 & -0.000\\-0.000 & -0.903\end{pmatrix}

Here it is clear that the convert matrix is not necessarily simple to find from the RLU matrix.

Tick marks
----------

The discussion of the location of the tick marks is quite long and deserves a section of its own. It all revolves about the usage of the curve linear axis, which is provided as an experimental/development feature in the *matplotlib* package. This in change then introduces some difficulties with the compability of the code with both python 2 and 3 as some calls are only supported with python 3. Thus to keep the compability of MJOLNIR to both versions, the customization of tick locations and number is only used for python versions 3. 

 Codewise, what is done is to create two subclasses of the *mticker.MultipleLocator* and the *mticker.MaxLocator* classes designed to calculate the positions of tick marks using two different methods. The *MaxLocator* expects a integer input signifying the maximal number of ticks along the given axis, which are chosen at 'suitable' positions and will update with zooming and panning. 
 
 The *MultipleLocator* is designed to place tick marks at multiples of a given value. That is, if the base is 0.1, tick marks are located at all integer multiples of 0.1 inside of the viewable area. For static or panning, this solution is suitable and sufficient, but allowing the user to zoom requires the class to be more complex. In order to find a suitable base value, the class finds the window limits and given a wanted number of tick marks it finds the average wanted distance between the marks. As to be scale invariant this number is mapped to a value between 0.1 and 1, where it is compared to a predefined list of allowed fractions (1/1,1/2,1/3,1/4,1/5) where the closest is used.

 With the new tick marks found, a hook is created to overview the panning and zooming of the axis and force an update of the drawing of tick marks. 