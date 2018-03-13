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

	for wedge in instr_root.getchildren():
		
		if wedge.tag in dir(Wedge):
			Wedgeclass_ = getattr(Wedge, wedge.tag)
		else:
			raise ValueError("Element is supposed to be a Wedge, but got '{}'.".format(wedge.tag))
		wedgeSettings = {}
		
		for attrib in wedge.keys():
			wedgeSettings[attrib]=np.array(wedge.attrib[attrib].strip().split(','),dtype=float)
			
		temp_wedge = Wedgeclass_(**wedgeSettings)
		#print(temp_wedge)
		
		
		
		for item in wedge.getchildren():
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
					itemSettings[attrib]=np.array(attribVal)	
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

def test_parseXML(): # Improve this test!
	Instr = parseXML('MJOLNIR/Geometry/dataTest.xml')
	
	assert(len(Instr.wedges)==2)
	assert(len(Instr.wedges[0].analysers)==5)
	assert(len(Instr.wedges[0].detectors)==5)
	assert(Instr.wedges[0].detectors[0].length==0.5)
	assert(Instr.settings['Author']=='Jakob Lass')
	assert(Instr.settings['Date']=='13/03/18')
	assert(Instr.settings['Name']=='PSI-CAMEA')
	
	testAnalyser = Instr.wedges[0].analysers[3]
	assert(testAnalyser.d_spacing==3.35)
	assert(testAnalyser.mosaicity==60)
	assert(testAnalyser.width==0.05)
	assert(testAnalyser.height==0.1)
	assert(np.all(testAnalyser.position==np.array([0,0.1,0],dtype=float)))
	assert(np.all(testAnalyser.direction==np.array([1,0,0],dtype=float)))
	
	testDetector = Instr.wedges[1].detectors[2]
	assert(testDetector.length==0.5)
	assert(testDetector.pixels==456)
	assert(testDetector.diameter==0.02)
	assert(np.all(testDetector.position==np.array([0,0.0,1],dtype=float)))
	assert(np.all(testDetector.direction==np.array([1,0,0],dtype=float)))

	testWedge = Instr.wedges[0]
	assert(np.all(testWedge.position==np.array([0,0.0,0],dtype=float)))

if __name__ == '__main__':
	import sys
		
	sys.path.append('.')
	sys.path.append('../..')
	import MJOLNIR
		
	Instr = parseXML('/home/lass/Dropbox/PhD/Software/MJOLNIR/MJOLNIR/Geometry/dataTest.xml')
	print(str(Instr))