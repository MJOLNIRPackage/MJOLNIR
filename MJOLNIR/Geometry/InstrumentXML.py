#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 16:43:47 2018

@author: lass
"""
import numpy as np
def parseXML(filename):
		
	from MJOLNIR.Geometry import Detector,Analyser,Wedge,Instrument
	import xml.etree.ElementTree as ET
	import numpy as np

	tree = ET.parse(filename)
	instr_root = tree.getroot()

	instrSettings = {}
		
	for attrib in instr_root.keys():
		instrSettings[attrib]=instr_root.attrib[attrib]

	Instr = Instrument.Instrument(**instrSettings)

	for wedge in list(instr_root):#.getchildren():
		
		if wedge.tag in dir(Wedge):
			Wedgeclass_ = getattr(Wedge, wedge.tag)
		else:
			raise ValueError("Element is supposed to be a Wedge, but got '{}'.".format(wedge.tag))
		wedgeSettings = {}
		
		for attrib in wedge.keys():
			if attrib=='concept':
				wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=str)
			else:		
				wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=float)
			
		temp_wedge = Wedgeclass_(**wedgeSettings)
		#print(temp_wedge)
		
		
		
		for item in list(wedge):#.getchildren():
			if item.tag in dir(Detector):
				class_ = getattr(Detector, item.tag)
			elif item.tag in dir(Analyser):
				class_ = getattr(Analyser,item.tag)
			else:
				raise ValueError("Item '{}' not recognized as MJOLNIR detector or analyser.".format(item.tag))
			
			itemSettings = {}
			for attrib in item.keys():
				attribVal = item.get(attrib).strip().split(',')
				if len(attribVal)==1:
					itemSettings[attrib]=float(attribVal[0])
				else:
					if(attrib=='split'):
						#print(type(attribVal))
						itemSettings[attrib]=attribVal
					else:
						itemSettings[attrib]=np.array(attribVal,dtype=float)	
			try:
				temp_item = class_(**itemSettings)
			except TypeError as e:
				print(e.args[0])
				raise ValueError('Item {} misses argument(s):{}'.format(class_,e.args[0].split(':')[1]))
			except ValueError:
				raise ValueError('Item {} not initialized due to error.'.format(class_))
			#print(temp_item)
			temp_wedge.append(temp_item)
			#print()

		#print(str(temp_wedge))
		Instr.append(temp_wedge)
	return Instr

def createXMLString(instrument):
	XMLString = '<?xml version="1.0"?>\n'
	XMLString+= '<Instrument '
	for attrib in instrument.settings:
		XMLString+="{}='{}' ".format(attrib,instrument.settings[attrib])
	XMLString+='>\n'
		
	for wedge in instrument.wedges:
		XMLString+="\t<Wedge "
		for attrib in wedge.settings:
			XMLString+="{}='{}' ".format(attrib,wedge.settings[attrib])
		XMLString+='>\n'

		for item in wedge.analysers + wedge.detectors:
			itemClass = str(item.__class__).split('.')[-1][:-2]
			XMLString+="\t\t<{}".format(itemClass)
			for key in item.__dict__:
				value = item.__getattribute__(key)
				if isinstance(value,type(np.array([0,0,0]))):
					valueStr = ','.join([str(x) for x in item.__getattribute__(key)])
				else:
					valueStr = str(value)
				XMLString+=" {}='{}'".format(str(key)[1:],valueStr)
			XMLString+="></{}>\n".format(itemClass)
			
		
		XMLString+="\t</Wedge>\n"
	XMLString+="</Instrument>\n"
	return XMLString
		

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



