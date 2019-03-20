FAQ
===

A collection of questions and answers ment to guide the user(s) of the MJOLNIR package through some of the most common pitfalls.


My RLU-axis has its tick marks positions stupidly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This 'feature' of the RLU axis is unfortunately linked to the current stage of the development of the `GridHelperCurveLinear <https://matplotlib.org/api/_as_gen/mpl_toolkits.axisartist.grid_helper_curvelinear.GridHelperCurveLinear.html#mpl_toolkits.axisartist.grid_helper_curvelinear.GridHelperCurveLinear>`_.
Because of the current lack of customization for this in-progress work some attributes cannot be changed. However, when generating an new RLUaxis 
one can specify the number of tick marks on the x and y axis. This number stays constant  throughout the lifetime of that axis, but a new can be created with a 
different number if needed. Further, specifying only the number of y-axis tick marks the program calculates a suitable number along the x-axis.


Everytime I run my script, I am bombarded with 'The file xx exists alread. Old file will be renamed to xx.'-warnings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The warning originates from the fact that a data file with the same name has been converted priviously in the sampe directory, 
and the program then tells you that the old conversion will be moved to a new file with the new name. It only serves as to try 
to prevent you from saving a conversion on top of another involuntarily. To  silence this annoyance, convert the data files only 
once or use the 'saveFile=False' flag in the DataSet.convertDataFile() method to not save the converted data.


