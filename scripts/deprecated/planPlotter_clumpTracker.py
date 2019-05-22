#!/usr/bin/python
import numpy as np
import matplotlib as m
#m.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
import planOutputReader as readerPlan
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy.polynomial.polynomial as poly
################################################################################
# paths and CL args
pathBase  = str(sys.argv[1])
nStart    = int(sys.argv[2])
nTot      = int(sys.argv[3])
makePlots = True
pathPlan  = pathBase + 'planOutput/'
pathSave  = pathBase + 'plots/clumpTracking/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
doPlan = readerPlan.DataPlan(pathPlan, nStart=nStart, nTot=nTot, nPar=512*512*512)
################################################################################
class Clump:
	massThresh = 1.e20
	posThresh  = 1.e-2
	failCountThresh = 1
	def __init__(self, do, n, nClump):
		print("initializing a clump data structure...")
		self.massList    = [do.peakArrayList[n][nClump, 2]]
		self.posList     = [np.asarray(do.peakArrayList[n][nClump, 4:7])]
		self.velList     = [np.asarray(do.peakArrayList[n][nClump, 7:10])]
		self.dt          = do.dt
		self.theoPosList = []
		self.existsList  = [True]
		self.nClumpList  = [nClump]
		self.deprecatedList = [False]
		self.failCount   = 0
		self.nStart      = n
	def getPeriodicTimeRight(self):
		nSim = self.nStart + len(self.massList) - 0.5
		tSim = nSim * self.dt
		pt   = np.mod(tSim, 2./3.)
		return pt
	def getPeriodicTimeLeft(self):
		nSim = self.nStart + len(self.massList) + 0.5
		tSim = nSim * self.dt
		pt   = np.mod(tSim, 2./3.)
		return pt
	def getShearVelocity(self, pos):
		return pos[0] * (-1.5)
	def xBC(self, theoPos):
		if theoPos[0] >  0.1:
			theoPos[0] -= 0.2
			theoPos[1] -= 1.5 * 1.0 * 0.2 * self.getPeriodicTimeRight()
		elif theoPos[0] < -0.1:
			theoPos[0] += 0.2
			theoPos[1] -= 1.5 * 1.0 * 0.2 * self.getPeriodicTimeLeft()
		return theoPos
	def yBC(self, theoPos):
		if theoPos[1] >  0.1:
			theoPos[1] -= 0.2
		elif theoPos[1] < -0.1:
			theoPos[1] += 0.2
		return theoPos
	def updateTheoPos(self):
		# simple update if not new
		shearVel   = self.getShearVelocity(self.posList[-1])
		totVel     = self.velList[-1]
		totVel[1] += shearVel
		theoPos    = self.posList[-1] + self.dt * totVel
		theoPos    = self.xBC(theoPos)
		theoPos    = self.yBC(theoPos)
		self.theoPosList.append(theoPos)
	def findClump(self, do, n):
		if not self.deprecatedList[-1]:
			self.updateTheoPos()
			nClumps = len(do.peakArrayList[n])
			posDiffArr  = np.zeros(nClumps)
			massDiffArr = np.zeros(nClumps)
			for nClump in range(nClumps):
				peak     = do.peakArrayList[n][nClump]
				pos      = np.asarray(peak[4:7])
				mass     = peak[2]
				vel      = np.asarray(peak[7:10])
				posDiff  = np.sum(np.absolute(self.theoPosList[-1] - pos))
				massDiff = np.absolute(mass - self.massList[-1]) / self.massList[-1]
				posDiffArr[nClump]  = posDiff
				massDiffArr[nClump] = massDiff
			nClumpBestMatchPos  = np.argmin(posDiffArr)
			nClumpBestMatchMass = np.argmin(massDiffArr)
			if posDiffArr[nClumpBestMatchPos] < self.posThresh and massDiffArr[nClumpBestMatchPos] < self.massThresh:
				nBest = nClumpBestMatchPos
				matchFound = True
			else:
				matchFound = False
			#print("minimum diff found at clump number " + str(nClumpBestMatch))
			if matchFound:
				peak     = do.peakArrayList[n][nBest]
				pos      = np.asarray(peak[4:7])
				mass     = peak[2]
				vel      = np.asarray(peak[7:10])
				self.posList.append(pos)
				self.massList.append(mass)
				self.velList.append(vel)
				self.nClumpList.append(nBest)
				self.existsList.append(True)
				self.failCount = 0
				self.deprecatedList.append(False)
				print("match found at n=" + str(n))
				#print("pos was:       " + str(self.posList[-2]))
				#print("velocity was:  " + str(self.velList[-2]))
				#print("expected pos:  " + str(self.theoPosList[-1]))
				#print("actual pos:    " + str(self.posList[-1]))
				#print("previous mass: " + str(self.massList[-2]))
				#print("current mass:  " + str(self.massList[-1]))
				print("pos diff:      " + str(posDiffArr[nBest]))
				print("mass diff:     " + str(massDiffArr[nBest]))
			else:
				print("failed to find clump at n=" + str(n))
				#print("diff array was: ")
				self.posList.append(self.theoPosList[-1])
				self.massList.append(self.massList[-1])
				self.velList.append(self.velList[-1])
				self.nClumpList.append(-1)
				self.existsList.append(False)
				self.failCount += 1
				print("expected pos was " + str(self.theoPosList[-1]))
				print("min pos diff was " + str(posDiffArr[nClumpBestMatchPos]))
				print("with mass diff   " + str(massDiffArr[nClumpBestMatchPos]))
				print("setting pos to   " + str(self.posList[-1]))
				print("setting vel to   " + str(self.velList[-1]))
				print("setting mass to  " + str(self.massList[-1]))
				if self.failCount == self.failCountThresh:
					print("failed to find clump too many times in a row")
					print("deprecating this clump")
					self.deprecatedList.append(True)
				else:
					self.deprecatedList.append(False)
		else:
			self.deprecatedList.append(True)

nEnd = doPlan.nTot
print("##################################################################")
# make a clump object for each clump that exists at the starting time
clumpObjList = []
for nClump in range(len(doPlan.peakArrayList[nStart])):
	clumpObjList.append(Clump(doPlan, nStart, nClump))
print("found " + str(len(clumpObjList)) + " clumps to start at n=" + str(nStart))

# try to track clumps over time
for n in range(nStart+1, nEnd):
	claimedClumps = []
	for clump in clumpObjList:
		clump.findClump(doPlan, n)
		claimedClumps.append(clump.nClumpList[-1])
	if n<nEnd-1:
		for nClump in range(len(doPlan.peakArrayList[n])):
			if nClump not in claimedClumps:
				clumpObjList.append(Clump(doPlan, n, nClump))

# make plots
if makePlots:
	for n in range(nStart, nEnd):
		nMatched=0; nNew =0; nDep=0; nFailed=0; nEver=0;
		# plot the starting clumps over time
		for nClump in range(len(clumpObjList)):
			clump = clumpObjList[nClump]
			nEff  = n-clump.nStart
			if nEff>=0:
				nEver+=1
				if not clump.deprecatedList[nEff]:
					if clump.existsList[nEff]:
						if n>clump.nStart:
							nMatched+=1
							#plt.scatter(clump.posList[nEff-1][0], clump.posList[nEff-1][1], s=80, facecolors='none', edgecolors='b', linestyle=':')
							plt.scatter(clump.posList[nEff][0], clump.posList[nEff][1], s=80, facecolors='none', edgecolors='b')
							plt.scatter(clump.theoPosList[nEff][0], clump.theoPosList[nEff][1], s=80, facecolors='none', edgecolors='b', linestyle='--')
						else:
							plt.scatter(clump.posList[nEff][0], clump.posList[nEff][1], s=80, facecolors='none', edgecolors='r')
							plt.scatter(clump.theoPosList[nEff][0], clump.theoPosList[nEff][1], s=80, facecolors='none', edgecolors='r', linestyle='--')
							nNew+=1
					else:
						nFailed+=1
						if n>clump.nStart:
							a=1
							#plt.scatter(clump.posList[nEff-1][0], clump.posList[nEff-1][1], s=80, facecolors='none', edgecolors='gray', linestyle=':')
						plt.scatter(clump.posList[nEff][0], clump.posList[nEff][1], s=80, facecolors='none', edgecolors='gray')
						plt.scatter(clump.theoPosList[nEff][0], clump.theoPosList[nEff][1], s=80, facecolors='none', edgecolors='gray', linestyle='--')
					plt.text(clump.posList[nEff][0]+0.002, clump.posList[nEff][1]+0.002, str(nClump)+","+str(clump.nClumpList[nEff]), fontsize=5)
				else:
					nFailed+=1
					nDep+=1
		# scatter plot for all clumps
		masses = doPlan.peakArrayList[n][:, 2]
		xs     = doPlan.peakArrayList[n][:, 4]
		ys     = doPlan.peakArrayList[n][:, 5]
		zs     = doPlan.peakArrayList[n][:, 6]
		sizes  = [np.power(1.e4*mass, 1./2.) for mass in masses]
		plt.scatter(xs, ys, s=sizes)

		# other plot stuff
		plt.text(0.102, 0.05,  "looked for: " + str(nMatched+nFailed))
		plt.text(0.102, 0.04,  "matched: " + str(nMatched))
		plt.text(0.102, 0.03,  "failed: " + str(nFailed))
		plt.text(0.102, 0.02,   "new: " + str(nNew))
		plt.text(0.102, 0.01, "deprecated: " + str(nDep))
		plt.text(0.102, 0.00, "total ever: " + str(nEver))
		plt.xlabel(r'$r/h$')
		plt.ylabel(r'y/h')
		plt.ylim(-0.1, 0.1)
		plt.xlim(-0.1, 0.1)
		print("saving plot for n=" + str(n))
		tools.saveAndClear(pathSave + "clumpTracking_"+ str(n) +".png", figNum=0)


# assemble final stats
initMassList = []
n = 260
for clump in clumpObjList:
	nEff = n - clump.nStart
	if not clump.deprecatedList[nEff] and clump.existsList[nEff]:
		initMassList.append(clump.massList[0])


tools.saveAndClear(pathSave + "clumpTracking_"+ str(n) +".png", figNum=0)

















#
