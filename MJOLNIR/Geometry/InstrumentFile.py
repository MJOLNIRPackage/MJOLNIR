from MJOLNIR.Geometry import Detector,Analyser,Wedge
import xml.etree.ElementTree as ET

#delimiterChar = ':'
commentChar = '#'
specialChar = '_'
PropertyCharStart = '('
PropertyCharEnd = ')'

def InstrumentFile(fileContent):
    settings = []

    #file = open(fileName, 'r')
    #fileContent = file.readlines()
    #file.close()

    for i in range(len(fileContent)):
        strippedLine=fileContent[i].strip()
        print(strippedLine)
        print('**************')
        if not strippedLine or strippedLine[0]==commentChar:    # Check if string is empty
            continue
        elif strippedLine[0]==specialChar:
            #try:
            [pos,objectType] = GetWedge(strippedLine[1:]) # Get name and value without special character
            
            #except ValueError:
            #    raise ValueError("Line {} does not contain a readable format!".format(i))

            if(objectType in dir(Detector) or objectType in dir(Analyser)):
                print(pos,objectType)
            else:
                raise ValueError('Object type \'{}\' in line {} of instrument file not recognized as correct object.'.format(objectType,i+1))

def GetWedge(string):
    """Returns position and type of object for a wedge declaration"""
    #try:
    #position = string.index(delimiterChar)
    #except ValueError:
    #    raise ValueError("Line does not contain the delimiter character '{}'".format(delimiterChar))
    #else:
    #variableName = string[:position].strip()    # Variable name is assumed to be in front of delimiterChar
    try:
        commentPosition = string.index(commentChar)
    except ValueError:
        commentPosition = len(string)

    
    variableContent = string[:commentPosition].strip()
    try:
        startPos = variableContent.index(PropertyCharStart)
        endPos = variableContent.index(PropertyCharEnd)
    except:
        pos = (0.0,0.0,0.0)
    else:
        pos = [float(x) for x in variableContent[startPos+1:endPos].split(',')]
    return [pos,variableContent[:startPos]]

def test_InstrumentFileReader():
    text = ['','# Test   ','_Analyser(0,0,0)',' \n_']

    print(InstrumentFile(text))