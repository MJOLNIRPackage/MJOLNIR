Voronoi Tesselation and plotting functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The method plotA3A4 makes some assumptions about the physics in the measured data. That is, it 
assumes that all A4 values are the same for all of the pixel from one detector. This is however 
only true for the middle detector in the lower layer as it is positioned radially with the flat 
analysers aligned to it. The other detectors have slightly different A4 values for different 
pixels but the change is in the range of 0.5 degree.

Voronoi Tesselation
---------------------
A voronoi diagram is in e.g. 2D a plot showing the area around all of the provided points that 
are closer to one initial point than any other. In higher dimensions this corresponds to 
volumes and hyper-volumes closest to the given points. 

If one performed a voronoi tesselation on the reciprocal lattice points one would get exactly
 the Brillouin zone as these zones are defined to contain all points closer to its center 
 than any other center.

Walk-through of method
---------------------
Given a set of data files the program first asserts that the files are compatible. This is 
defined to be when the files have the same sample, and approximately same physical measurement 
conditions. That is, the incoming energy, temperature and magnetic and electric fields are to 
be equal within a tolerance. If all of this is satisfied, the files can be added together.

The measured A3 and A4 points for each file is calculated from the normalization tables and 
added to a common list of vertices needed for the voronoi algorithm. However, before the 
tesselation is performed, the boundary of each file is found from expanding the measurement 
points in both direction with a step size equal to half the size between the edge points and 
second rows. This will seek to make the edge pixels the same size as the second row pixels. 
Furthermore, as to stop the edge pixel to extend to infinity 9 additional points are 
generated. They are positioned around the data with one above the average x position, one 
below, one to each side, and also diagonally from the data. Thus the measurement points are 
fenced by the extra points and an edge pixel extending to infinity should be impossible. Then 
the voronoi is generated as described above resulting in n+10 regions, where n is the number 
of measured positions. For each pixel it is checked whether it intersects the boundary of the 
initial points; if not it is an internal point and a polygon patch can be created; else the 
overlap between the pixel and boundary polygon is kept and made into a patch.

One should now have n pixels with centroids at approximately the same position as the initial 
points. Next the intensity or colour is to be found. This could naively be done by just looping 
through the individual points and finding the corresponding pixel that contains it. However, 
this is simple too time-consuming and a trick is to be used instead. The trick is to sort 
intensities corresponding to the position of the measurement point such that smallest x and y 
values are first, with largest x and y values last. The same is done with the patches. Thus, 
one can then colour the pixels directly from the sorted intensity without any further difficulty. 
The only downside of this trick is that there has to be a one-to-one correspondance between each 
measurement point and patch. This might seem trivial but if a position is measured twice, that 
is, the same A3 and A4 are covered precisely by two different points the algorithm fails.

Dual of Delaunay triangulation
Given a metric it divides the space into regions where each point is closer to the given vertex than any other vertex.


Pictures:
Initial points plus boundary
Voronoi to infinity
Addition of surrounding points
Cutting of outer regions with boundary
Final picture

