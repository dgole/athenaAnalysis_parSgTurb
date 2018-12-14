#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import matplotlib.colors as colors

################################################################################
# Data class ###################################################################
################################################################################
class DataPlan:
	def __init__(self, path, dt=0.1):
		print("initializing PLAN data structure from " + path)
		self.path = path
		contentsList = os.listdir(path)
		self.peakFileNameList = []
		for item in contentsList:
			if item[:5]=='peaks': self.peakFileNameList.append(item)
		print('finding ' + str(len(self.peakFileNameList)) + ' peak files')
		self.peakArrayList = [0 for n in range(len(self.peakFileNameList))]
		self.timeList = []
		n = 0
		for peakFileName in self.peakFileNameList:
			n1   = int(peakFileName[9:12])
			n2   = int(peakFileName[13])
			time = float(n1) + dt * float(n2)
			n    = int(np.round(time/dt))
			self.peakArrayList[n] = np.loadtxt(self.path+peakFileName)
		self.timeList = [n*dt for n in range(0,len(self.peakFileNameList))]
		self.nClumpsList = []
		for arr in self.peakArrayList: self.nClumpsList.append(arr.shape[0])
		print('number of clumps found at each data dump: ')
		print(self.nClumpsList)
		self.nClumps = np.asarray(self.nClumpsList)
		self.time    = np.asarray(self.timeList)

def getCumMassHist(do, n, spacing=0.25):
	if do.nClumps[n] != 0:
		masses = []
		for mass in do.peakArrayList[n][:,2]:
			masses.append(mass)
		masses = np.asarray(masses)
		indexStart = np.floor(np.log10(np.amin(masses))) - spacing
		indexEnd   = np.ceil(np.log10(np.amax(masses)))
		lLimList = [np.power(10,index) for index in np.arange(indexStart, indexEnd, spacing)]
		nBins    = len(lLimList)
		massHist = np.zeros(nBins)
		for mass in masses:
			for n in range(nBins):
				if lLimList[n] < mass:
					massHist[n]+=1
		return lLimList, massHist
	else:
		print('no partilce found in this output')
		return 0.0, 0.0

####################################################################
# plotting functions ###############################################
####################################################################









































#
