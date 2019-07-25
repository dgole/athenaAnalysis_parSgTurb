#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors
import sys
sys.path.append('../python')
import athenaTools as tools

thresh = 0.5
class Clump:
	def __init__(self, inArr, doPlan, splitterIds, splitterData):
		#print("###########################################################")
		#print("initializing clump data structure...")
		self.nArr  = inArr[:,0]
		self.idArr = inArr[:,1]
		self.n0    = self.nArr[0]
		self.nLast = self.nArr[-1]
		self.massDict   = {}
		self.posDict    = {}
		self.splitter   = False
		self.splitFrac  = 0.0
		self.persistence = self.nArr.shape / (self.nLast - self.n0 + 1.0)
		nId = 0
		for n in self.nArr:
			for peak in doPlan.peakArrayList[n]:
				if peak[0] == self.idArr[nId]:
					self.massDict[n] = peak[2]
					self.posDict[n]  = peak[4:7]
			nId += 1
		if len(self.massDict.values()) != self.nArr.shape[0]:
			a=1
			#print("couldn't match this clump to a peak at all frames")
		#print("looking for splitters for inArr[0]=" + str(inArr[0]))
		maxSplitFrac = 0.0
		for i in range(splitterIds.shape[0]):
			if splitterIds[i,0] == inArr[0,0] and splitterIds[i,1] == inArr[0,1]:
				loc = i
				splitFrac  = splitterData[loc,0]
				splitFrac1 = splitterData[loc,1]
				#print("at least partial splitter found")
				#print("split fraction is: " + str(splitFrac))
				#print("other split frac is: " + str(splitFrac1))
				#print(splitterIds[i])
				if splitFrac > maxSplitFrac: maxSplitFrac = splitFrac
		self.splitFrac = maxSplitFrac
		#print("split fraction is: " + str(self.splitFrac))
		if self.splitFrac > thresh:
			#print("marking this clump as a splitter")
			self.splitter = True
		else:
			a=1
			#print("NOT marking this clump as a splitter")


def get_initialMasses(clumpObjList):
	massList = []
	for clump in clumpObjList:
		massList.append(clump.massArr[0])
	return massList

def get_currentMasses(clumpObjList, n):
	massList = []
	for clump in clumpObjList:
		if n in clump.massDict.keys():
			massList.append(clump.massDict[n])
		else:
			massList.append(1.e-6)
	return massList

def get_xs(clumpObjList, n):
	xList = []
	for clump in clumpObjList:
		if n in clump.posDict.keys():
			xList.append(clump.posDict[n][0])
		else:
			xList.append(0.1)
	return xList

def get_ys(clumpObjList, n):
	yList = []
	for clump in clumpObjList:
		if n in clump.posDict.keys():
			yList.append(clump.posDict[n][1])
		else:
			yList.append(-0.1)
	return yList

def get_colors(clumpObjList, n):
	colorList = []
	for clump in clumpObjList:
		if clump.splitter:
			colorList.append('b')
		else:
			colorList.append('r')
	return colorList
