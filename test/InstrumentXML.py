from MJOLNIR.Geometry.InstrumentXML import createXMLString, parseXML


def test_parseXML(): # Improve this test!
	from MJOLNIR.Geometry import Wedge,Analyser,Detector,Instrument
	import os
	tempFileName = '__temp__'
		
	Instr = Instrument.Instrument()
	Instr.settings['Author'] = 'Jakob Lass'

	wedge = Wedge.Wedge(position=(0.5,0,0))

	Det = Detector.TubeDetector1D(position=(1.0,1,0),direction=(1,0,0))
	Ana = Analyser.FlatAnalyser(position=(0.5,0,0),direction=(1,0,1))

	wedge.append([Det,Ana])
	Instr.append([wedge,wedge])
	Instr.append(wedge)
		
	f = open(tempFileName,'w')

	f.write(createXMLString(Instr))
	f.close()
		
		
	InstrLoaded = parseXML(tempFileName)
	os.remove(tempFileName)
		
	assert(Instr==InstrLoaded)


def test_parseWrongXML():
	import os
	tempFileName = '__temp__'
	with open(tempFileName,'w') as f:
		f.write("<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>n"\
+"	<WrongWedge position='0.0,0.0,0.0' concept='ManyToMany'>\n"\
+"		<FlatAnalyser position='-0.0,1.1195,0.0' direction='0.7071067811865476,-0.0,0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>\n"\
+"		<FlatAnalyser position='-0.0,1.1827,0.0' direction='0.7071067811865475,-0.0,0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>\n"\
+"		<TubeDetector1D position='-0.01151899615488017,1.1999447123628586,0.71' direction='-0.01151899615488017,1.1999447123628586,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"		<TubeDetector1D position='-0.0,1.2,0.7' direction='-0.0,1.2,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"		<TubeDetector1D position='0.01151899615488017,1.1999447123628586,0.71' direction='0.01151899615488017,1.1999447123628586,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"	</WrongWedge>\n"\
+"</Instrument>")

	try:
		InstrLoaded = parseXML(tempFileName)
	except ValueError: # Wrong wedge type
		assert True

	os.remove(tempFileName)

	with open(tempFileName,'w') as f:
		f.write("<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>n"\
+"	<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>\n"\
+"		<WrongFlatAnalyser position='-0.0,1.1195,0.0' direction='0.7071067811865476,-0.0,0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></WrongFlatAnalyser>\n"\
+"		<FlatAnalyser position='-0.0,1.1827,0.0' direction='0.7071067811865475,-0.0,0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>\n"\
+"		<TubeDetector1D position='-0.01151899615488017,1.1999447123628586,0.71' direction='-0.01151899615488017,1.1999447123628586,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"		<TubeDetector1D position='-0.0,1.2,0.7' direction='-0.0,1.2,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"		<TubeDetector1D position='0.01151899615488017,1.1999447123628586,0.71' direction='0.01151899615488017,1.1999447123628586,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"	</Wedge>\n"\
+"</Instrument>")
	
	try:
		InstrLoaded = parseXML(tempFileName)
	except ValueError: # Wrong Analyser type
		assert True

	os.remove(tempFileName)

	with open(tempFileName,'w') as f:
		f.write("<Instrument Initialized='False' Author='Jakob Lass' Date ='16/03/18' position='0.0,0.0,0.0'>n"\
+"	<Wedge position='0.0,0.0,0.0' concept='ManyToMany'>\n"\
+"		<FlatAnalyser direction='0.7071067811865476,-0.0,0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>\n"\
+"		<FlatAnalyser position='-0.0,1.1827,0.0' direction='0.7071067811865475,-0.0,0.0' d_spacing='3.354' mosaicity='60' width='0.05' height='0.1'></FlatAnalyser>\n"\
+"		<TubeDetector1D position='-0.01151899615488017,1.1999447123628586,0.71' direction='-0.01151899615488017,1.1999447123628586,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"		<TubeDetector1D position='-0.0,1.2,0.7' direction='-0.0,1.2,0.0' pixels='1024' length='1' diameter='0.02' split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"		<TubeDetector1D split='0,189, 298, 404, 510, 618, 726, 837,1024'></TubeDetector1D>\n"\
+"	</Wedge>\n"\
+"</Instrument>")
	
	try:
		InstrLoaded = parseXML(tempFileName)
	except ValueError: # No position in flat analyser
		assert True

	os.remove(tempFileName)