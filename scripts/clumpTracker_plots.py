#!/usr/bin/python
import numpy as np
#import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaTools as tools
import planOutputReader as readerPlan
################################################################################
# paths and CL args
pathBase  = str(sys.argv[1])
nStart    = int(sys.argv[2])
nStop     = int(sys.argv[3])
pathPlan  = pathBase + 'planOutput2/'
pathCT    = pathBase + 'planOutput2/clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
pathSave  = pathBase + 'plots/clumpTracking2/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
class Clump:
	def __init__(self, inArr, doPlan):
		print("initializing clump data structure...")
		self.nArr  = inArr[:,0]
		self.idArr = inArr[:,1]
		self.n0    = self.nArr[0]
		self.nLast = self.nArr[-1]
		self.massDict   = {}
		self.posDict    = {}
		nId = 0
		for n in self.nArr:
			for peak in doPlan.peakArrayList[n]:
				if peak[0] == self.idArr[nId]:
					self.massDict[n] = peak[2]
					self.posDict[n]  = peak[4:7]
			nId += 1
		if len(self.massDict.values()) != self.nArr.shape[0]:
			print("couldn't match this clump to a peak at all frames")

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

################################################################################

# read peak files
doPlan = readerPlan.DataPlan(pathPlan, nStart=nStart, nTot=nStop, nPar=512*512*512)
clumpObjList = []
ctFileNameList = os.listdir(pathCT)
for fileName in ctFileNameList:
	inArr = np.load(pathCT + fileName)
	clumpObjList.append(Clump(inArr, doPlan))

#for clump in clumpObjList:
	#plt.semilogy(clump.nArr, clump.massArr)
#plt.show(); plt.clf()

for n in range(nStart, nStop):
	plt.figure(0)
	plt.figure(num=0, figsize=(3,2))
	xs     = get_xs(clumpObjList, n)
	ys     = get_ys(clumpObjList, n)
	masses = get_currentMasses(clumpObjList, n)
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	plt.scatter(xs, ys, s=sizes)
	for i in range(len(clumpObjList)):
		plt.text(xs[i]+0.001, ys[i]+0.001, str(i), fontsize=5)
	plt.xlabel(r'$r/h$')
	plt.ylabel(r'y/h')
	plt.ylim(-0.1, 0.1)
	plt.xlim(-0.1, 0.1)
	plt.tight_layout()
	tools.saveAndClear(pathSave + "scatter_" + str(n) + ".png", figNum=0)
	plt.close('all')



















#
