import sys 
sys.path.append('/home/lass/Dropbox/PhD/Software/MJOLNIR/')

from Tutorial_Class import Tutorial


def Tester():
    from MJOLNIR.Geometry import Instrument,Detector,Analyser,Wedge
    # Create instrument
    Instr = Instrument.Instrument()

    # Create two detectors and their corresponding analysers
    Det = Detector.TubeDetector1D(position=(1,0,1),direction=(1,0,0),pixels=10,split=[2,5,8],length=1.0)
    Det2 = Detector.TubeDetector1D(position=(1,0.1,1),direction=(1,0,0),pixels=10,split=[2,5,8],length=1.0)
    Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))
    Ana2 = Analyser.FlatAnalyser(position=(1,0,0),direction=(1,0,1))

    # Collect everything in the wedge object
    wedge = Wedge.Wedge(position=(0,0,0),detectors=[Det,Det2],analysers=[Ana,Ana2])

    Instr.append(wedge)

    Instr.initialize()
    print(Instr.A4)
    print(Instr.Ef)
    
title = 'Calculate A4 and Ef'

introText = 'With a virtual representation of the instrument, as explained in the `Instrument <instrument.rst#Build-simple-instrument>`__ tutorial, the first thing is to calculate the corresponding A4 and Ef '\
+'values. This is done by simply accessing the corresponding attributes of the instrument object. The instrument needs to have it initialize method called such that it calculates the A4 and Ef values when a wedge has been appended.'

outroText = 'The values printed are given below. '\
+'[[array([-1.57079633, -1.57079633, -1.57079633, -1.57079633, -1.57079633,\n'\
+'       -1.57079633]), array([-1.50556269, -1.5067601 , -1.50824438, -1.52086906, -1.52111537,\n'\
+'       -1.52159382])]]\n'\
+'\n'\
+'[[array([4.81164763, 5.44262522, 6.18119629, 3.83622353, 4.27947022,\n'\
+'       4.81164763]), array([4.83218148, 5.46500461, 6.20544497, 3.84580212, 4.2900502 ,\n'\
+'       4.82331491])]]\n\n'\
+'Notice that the number of pixels and how these are splitted are given to the detector. What happens is that the detector splits '\
+"it's pixels into two lists; pixels before and after 2 and 8 respectively are discarted as non-active, while pixels 2-4 and 5-7 corresponds to the two areas 'seeing' the two detectors. "\
+'That is, the "ManyToMany" concept has been given as default. Another concept exists, the "OneToOne", as tabulated below.\n\n'\
+'+------------+---------------------------------------------------+--------------------------------------------------------------------+\n'\
+'| Concept    | Effect                                            | Requirement                                                        |\n'\
+'+============+===================================================+====================================================================+\n'\
+'| OneToOne   | Each detector "sees" only one analyser.           | The number of detectors and analysers in the wedge is to be equal. |\n'\
+'+------------+---------------------------------------------------+--------------------------------------------------------------------+\n'\
+'| ManyToMany | Each detector "sees" all analysers in the wedge.  | The number of splits is to be equal to the number of analysers +2, |\n'\
+'|            |                                                   | where the extra 2 comes from discarting non-active ends.           |\n'\
+'+------------+---------------------------------------------------+--------------------------------------------------------------------+'


introText = title+'\n'+'^'*len(title)+'\n'+introText


    
Example = Tutorial('A4Ef',introText,outroText,Tester,fileLocation = '/home/lass/Dropbox/PhD/Software/MJOLNIR/docs/Tutorials/Instrument')

def test_A4Ef():
    Example.test()

#if __name__ == '__main__':
Example.generateTutorial()
