import sys 
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from Tutorial_Class import Tutorial


def Tester():
    #Simple function written to show how things work
    import matplotlib.pyplot as plt
    import numpy as np
    
    x = np.linspace(0,10,101)
    y = np.sin(x)
    
    fig,ax = plt.subplots()
    ax.plot(x,y)
    
    fig.savefig('/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials2/DataSet/Here.png',format='png')
    
title = 'Template for tutorial generation'

introText = 'This is a short and simple template for creating a tutorial for the MJOLNIR package. What is needed for the developer to define' \
    +' is a title, an intro text, an outro text, and then the actual code to be formated. All of this is done by filling in the template and '\
    +'running the bash script "tutorials" found in the same folder. This initiates testing of the code written and a restructured text file '\
    +'is saved in the indicated folder. As an example, the code below can be used:'

outroText = 'This is just a simple example code generating a figure denoted "figure0.png". However, in the actual code this figure is called "MJOLNIR/docs/Tutorials2/DataSet/Here.png" ' \
+'and can be shown with the use of rst commands\n .. figure:: Here.png\n  :width: 30%\n  :align: center\n'\
+'\n All of this is to work, I hope.'

introText = title+'\n'+'^'*len(title)+'\n'+introText


    
Example = Tutorial('template',introText,outroText,Tester,fileLocation = '/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials2/DataSet')

def test_template():
    Example.test()

#if __name__ == '__main__':
Example.generateTutorial()
