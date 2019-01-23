Template for tutorial generation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a short and simple template for creating a tutorial for the MJOLNIR package. What is needed for the developer to define is a title, an intro text, an outro text, and then the actual code to be formated. All of this is done by filling in the template and running the bash script "tutorials" found in the same folder. This initiates testing of the code written and a restructured text file is saved in the indicated folder. As an example, the code below can be used:

>>> #Simple function written to show how things work
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> 
>>> x = np.linspace(0,10,101)
>>> y = np.sin(x)
>>> 
>>> fig,ax = plt.subplots()
>>> ax.plot(x,y)
>>> 
>>> fig.savefig('figure0.png',format='png')
>>>  

This is just a simple example code generating a figure denoted "figure0.png". However, in the actual code this figure is called "MJOLNIR/docs/Tutorials/Tools/Here.png" and can be shown with the use of rst commands
 .. figure:: Here.png
  :width: 30%
  :align: center
