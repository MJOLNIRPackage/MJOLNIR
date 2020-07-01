Optimizing of the plotA3A4 routine
----------------------------------
In order to make the method *plotA3A4* useful in a practical manner, 
it needs to be somewhat fast in its computation. Otherwise, one will 
not use it and instead use binning, or even worse another program! 
And as it is already mentioned in the documentation, the computation 
time for generating all pixels individually for the methods plotA3A4 
and plotQPatches is too long for practical use. Thus, an optimization 
of the underlying algorithm made sense to pursue. Before changing any 
code, the end test was set up: By first running the old method and 
then the new, two results are generated and they are checked to be 
identical. Is so, and the with a speed-up, the goal is reached.


Thus, before headlessly trying to perform optimization, a test of the 
current speed is needed. The following is a dump of times for the 
un-optimized function using four files and one plane with 8 
sub-planes (T0Phonon10meV.nxs, T0Phonon10meV90A4InterlaceA3.nxs, 
T0Phonon10meV93_5A4.nxs, and T0Phonon10meV93_5A4InterlaceA3.nxs, 
planes 8 through 15):

+-------------------------+----------------+------------------+
|Function                 | Time [s]       | Uncertainty [s]  |
+=========================+================+==================+
|testFiles                |  0.0003143     | 2.2458e-05       |  
+-------------------------+----------------+------------------+
|getA3A4                  |  2.8371810e-05 | 1.0929e-05       |  
+-------------------------+----------------+------------------+
|getData                  |  0.0561441     | 0.0118           |  
+-------------------------+-----+----------+------------------+
|concatINormMon           |  0.0405629     | 0.0068           |  
+-------------------------+----------------+------------------+
|A4Instr                  |  0.0001758     | 2.9835-05        |  
+-------------------------+----------------+------------------+
|genPointsAndBoundary     |  0.2238357     | 0.0151           |  
+-------------------------+----------------+------------------+
|voronoiTessellation      |  14.070171     | 1.6175           |  
+-------------------------+----------------+------------------+
|calcCentroids            |  0.4847603     | 0.0366           |  
+-------------------------+----------------+------------------+
|sortPoints               |  25.009635     | 0.7448           |  
+-------------------------+----------------+------------------+
|calculateQ               |  3.4248633     | 0.0778           |  
+-------------------------+----------------+------------------+
|binPlanes                |  0.0295210     | 0.0007           |  
+-------------------------+----------------+------------------+
|genPatchesAndCollection  |  8.7695337     | 0.0920           |  
+-------------------------+----------------+------------------+
|plotter                  |  0.1842771     | 0.0093           |
+-------------------------+----------------+------------------+


From this it is clear that three functions are to be looked at: *voronoiTessellation*, *sortPoints*, and *genPatchesAndCollection*.




Voronoi Tessellation subroutine
-------------------------------
Individual functions are timed by the using *my_timing_N* decorator. This allows 
for a partitioning of the program into smaller pieces that are easier to be 
optimized individually. Also when changing code in these sub-functions, the test 
was of course that the output is identical to the old one. This method is not the 
optimal for speeding up the code as it is most probable that changing larger coding 
structures might allow for a quicker speed-up. 

Below is a table of the resulting computational times for the tests used showing a clear speed-up of around 1.7 for real data. 

+-------------------------------+-----------------------+-----------------------+-----------------------+
| Tests (5 runs)                | Original [s]          | Optimized [s]         |  Gain                 |
+===============================+=======================+=======================+=======================+
| 100000 random points          | 10.27 :math:`\pm` 0.15| 4.45 :math:`\pm` 0.08 | 2.310 :math:`\pm` 0.05|
| between 0 and 1               |                       |                       |                       |
+-------------------------------+-----------------------+-----------------------+-----------------------+
| 1000000 random points         | 106.7 :math:`\pm` 0.8 | 47.75 :math:`\pm` 0.7 | 2.235 :math:`\pm` 0.04|
| between 0 and 1               |                       |                       |                       |
+-------------------------------+-----------------------+-----------------------+-----------------------+
|4 Data files and 8 planes ``*``| 10.93 :math:`\pm` 0.08| 6.3  :math:`\pm` 0.2  | 1.72 :math:`\pm` 0.06 |
+-------------------------------+-----------------------+-----------------------+-----------------------+

``*``: set-up from above.

The speed-up comes mainly from two changes; when generating the return data the 
intersections points of all of the polygons where recalculated while only the 
polygons cut by the boundary needed to be calculated. Secondly, when testing if 
all data points are within the boundary a list comprehension is changed into the 
vectorized *contains* function from the *shapely.vectorized* sub-library.
