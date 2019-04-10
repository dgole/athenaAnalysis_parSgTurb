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

################################################################################
# Data class ###################################################################
################################################################################
class DataPlan:
	def __init__(self, path, dt=0.1, nPar=3.e5, G=0.05):
		print("initializing PLAN data structure from " + path)
		self.path  = path
		self.nPar  = nPar
		self.G     = G
		self.m0_ceres = G * (720.0/0.05)
		contentsList = os.listdir(path)
		self.peakFileNameList = []
		for item in contentsList:
			if item[:5]=='peaks': self.peakFileNameList.append(item)
		print('finding ' + str(len(self.peakFileNameList)) + ' peak files')
		self.peakArrayList = [0 for n in range(len(self.peakFileNameList)-1)]
		self.timeList = []
		n = 0
		for peakFileName in self.peakFileNameList:
			n1   = int(peakFileName[9:12])
			n2   = int(peakFileName[13])
			time = float(n1) + 1.0 * float(n2)* 0.1
			n    = int(np.round(time/dt))
			n    = min(n, len(self.peakArrayList)-1)
			print(n)
			temp = np.loadtxt(self.path+peakFileName)
			if len(temp.shape) == 2:
				self.peakArrayList[n] = temp
			elif len(temp.shape) == 1 and temp.shape[0] == 13:
				self.peakArrayList[n] = temp.reshape((1, 13))
			else:
				self.peakArrayList[n] = np.zeros((1,13))
		self.peakArrayList = self.peakArrayList[:-1]
		for n in range(len(self.peakArrayList)):
			self.peakArrayList[n][:,2] *= self.m0_ceres
		self.timeList = [n*dt for n in range(0,len(self.peakArrayList))]
		self.nClumpsList = []
		for arr in self.peakArrayList:
			if arr[0,0] == 0 and arr[0,1] == 0:
				self.nClumpsList.append(0)
			else:
				self.nClumpsList.append(arr.shape[0])
			n+=1
		print('number of clumps found at each data dump: ')
		print(self.nClumpsList)
		self.nClumps = np.asarray(self.nClumpsList)
		self.time    = np.asarray(self.timeList)
		if np.sum(self.nClumps)==0: print('NO CLUMPS AT ANY TIME!')
		mTestPar = self.peakArrayList[-1][0,2]
		nTestPar = self.peakArrayList[-1][0,1]
		self.mParTot = self.nPar * (mTestPar/nTestPar)
		for n in range(0,10000):
			if self.nClumps[n]>0:
				self.nFirstClump = n
				break
		self.nt = len(self.timeList)

def getCumMassHist(do, n, spacing=0.05):
	if do.nClumps[n] != 0:
		masses = []
		for mass in do.peakArrayList[n][:,2]:
			masses.append(mass)
		masses     = np.asarray(masses)
		indexStart = -4
		indexEnd   = 1
		lLimList   = [np.power(10,index) for index in np.arange(indexStart, indexEnd, spacing)]
		nBins      = len(lLimList)
		massHist   = np.zeros(nBins)
		for mass in masses:
			for n in range(nBins):
				if lLimList[n] < mass:
					massHist[n]+=1
		return lLimList, massHist
	else:
		print('no partilce found in this output')
		return 0, 0

def getDiffMassHist(do, n, spacing=0.05):
	if do.nClumps[n] != 0:
		masses = []
		for mass in do.peakArrayList[n][:,2]:
			masses.append(mass)
		masses     = np.asarray(masses)
		indexStart = -4
		indexEnd   = 1
		lLimList   = [np.power(10.0,index) for index in np.arange(indexStart, indexEnd, spacing)]
		uLimList   = [np.power(10.0,index+spacing) for index in np.arange(indexStart, indexEnd, spacing)]
		avgMassList= [np.power(10.0,index+spacing/2.0) for index in np.arange(indexStart, indexEnd, spacing)]
		nBins      = len(lLimList)
		dm         = np.zeros(nBins)
		dn         = np.zeros(nBins)
		for mass in masses:
			for n in range(nBins):
				if lLimList[n] < mass < uLimList[n]:
					dn[n]+=1
		for n in range(nBins): dm[n] = uLimList[n] - lLimList[n]
		return avgMassList, dn/dm
	else:
		print('no partilce found in this output')
		return 0, 0

####################################################################
# plotting functions ###############################################
####################################################################
def nClumpsTimeEvo(do, figNum=0, colorOption='k'):
	print('plotting number of clumps vs time for ' + do.path)
	plt.figure(figNum)
	plt.plot(do.time, do.nClumps, colorOption)
	plt.xlabel(r'$t \Omega$')
	plt.ylabel(r'$N_{clumps}$')

def plotCumMassHist(do, nStart=None, nEnd=None, spacing=0.05, figNum=0, legendLabel=None, colorOption='ko'):
	print('plotting time-averaged cumulative mass hist for ' + do.path)
	plt.figure(figNum)
	if nStart is None: nStart = do.nt // 2
	if nEnd   is None: nEnd   = do.nt
	count=0
	for n in range(nStart, nEnd):
		bins, hist = getCumMassHist(do, n, spacing=spacing)
		if not isinstance(bins, int):
			if count==0:
				masterHist  = hist
				plotBins    = bins
			else:
				masterHist += hist
			count+=1
	masterHist /= count
	plt.loglog(plotBins, masterHist, colorOption, label=legendLabel)
	plt.xlabel(r'$M_p$')
	plt.ylabel(r'$N(>M_p)$')

def plotDiffMassHist(do, nStart=None, nEnd=None, spacing=0.05, figNum=0, legendLabel=None, colorOption='ko'):
	print('plotting time-averaged cumulative mass hist for ' + do.path)
	plt.figure(figNum)
	if nStart is None: nStart = do.nt - 100
	if nEnd   is None: nEnd   = do.nt
	count=0
	for n in range(nStart, nEnd):
		bins, hist = getDiffMassHist(do, n, spacing=spacing)
		if not isinstance(bins, int):
			if count==0:
				masterHist  = hist
				plotBins    = bins
			else:
				masterHist += hist
			count+=1
	masterHist /= count
	plt.loglog(plotBins, masterHist, colorOption, label=legendLabel)
	plt.xlabel(r'$M_p$')
	plt.ylabel(r'$dN/dM_p$')
	minModelErr = 1.e40
	minModelC   = 10000
	for c in [np.power(10.0,i) for i in np.arange(-3,0)]:
		model = c * np.power(bins, -1.6)
		plt.loglog(bins, model, linestyle='--', color='k')

def planMassFracTimeEvo(do, figNum=0, legendLabel=None, colorOption='k'):
	print('plotting mass frac in planetesimals for ' + do.path)
	plt.figure(figNum)
	mFracList = []
	for n in range(0, len(do.peakArrayList)):
		if n < do.nFirstClump: mNow = 0
		else: mNow = np.sum(do.peakArrayList[n][:,2])
		mFracNow = mNow / do.mParTot
		mFracList.append(mFracNow)
	plt.plot(do.timeList, mFracList, colorOption, label=legendLabel)
	plt.xlabel(r'$t \Omega$')
	plt.ylabel(r'$M_{plan} / M_{par}$')

def scatterPlotXZ(do, n, figNum=0):
	print('plotting XZ scatter plot for ' + do.path + ', n=' +  str(n))
	plt.figure(figNum)
	masses = do.peakArrayList[n][:, 2]
	xs     = do.peakArrayList[n][:, 4]
	ys     = do.peakArrayList[n][:, 5]
	zs     = do.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	############################################################################
	plt.scatter(xs, zs, s=sizes)
	plt.xlabel(r'x')
	plt.ylabel(r'z')
	plt.ylim(-0.1, 0.1)
	plt.xlim(-0.1, 0.1)
	plt.title(r'$t=$' + str(np.round(do.timeList[n],1)))

def scatterPlotXY(do, n, figNum=0):
	print('plotting XY scatter plot for ' + do.path + ', n=' +  str(n))
	plt.figure(figNum)
	masses = do.peakArrayList[n][:, 2]
	xs     = do.peakArrayList[n][:, 4]
	ys     = do.peakArrayList[n][:, 5]
	zs     = do.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	############################################################################
	plt.scatter(xs, ys, s=sizes)
	plt.xlabel(r'x')
	plt.ylabel(r'y')
	plt.ylim(-0.1, 0.1)
	plt.xlim(-0.1, 0.1)
	plt.title(r'$t=$' + str(np.round(do.timeList[n],1)))

def scatterPlotXYZ(do, n, figNum=0):
	print('plotting XYZ scatter plot for ' + do.path + ', n=' +  str(n))
	plt.figure(figNum)
	masses = do.peakArrayList[n][:, 2]
	xs     = do.peakArrayList[n][:, 4]
	ys     = do.peakArrayList[n][:, 5]
	zs     = do.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	fig = plt.figure(figNum)
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(ys, xs, zs, s=sizes)
	ax.set_xlabel(r'$\phi$')
	ax.set_ylabel(r'$r$')
	ax.set_zlabel(r'$z$')
	ax.set_xlim(-0.1, 0.1)
	ax.set_ylim(-0.1, 0.1)
	ax.set_zlim(-0.1, 0.1)
	plt.title(r'$t=$' + str(np.round(do.timeList[n],1)))

def get_p(do, n):
	nplan   = do.nClumpsList[n]
	if nplan > 2:
		minMass = np.amin(do.peakArrayList[n][:,2])
		sum     = 0
		for mass in do.peakArrayList[n][:,2]:
			sum += np.log(mass / minMass)
		p   = 1 + nplan * np.power(sum, -1)
		err = (p-1)/np.sqrt(nplan)
		return p, err
	else:
		return 0.0, 0.0

def get_p_avg(do, nStart=None, nEnd=None):
	if nStart is None: nStart = do.nt-100
	if nEnd   is None: nEnd   = do.nt
	pList   = []
	errList = []
	for n in range(nStart, nEnd):
		p, err = get_p(do, n)
		pList.append(p)
		errList.append(err)
	pArr   = np.asarray(pList)
	errArr = np.asarray(errList)
	pAvg   = np.round(np.mean(pArr),   decimals=4)
	errAvg = np.round(np.mean(errArr), decimals=4)
	return pAvg, errAvg

def pValuePlot(do, nStart=None, nEnd=None, figNum=0):
	plt.figure(figNum)
	pList   = []
	errList = []
	for n in range(do.nFirstClump, do.nt):
		p, err = get_p(do, n)
		pList.append(p)
		errList.append(err)
	pArr   = np.asarray(pList)
	errArr = np.asarray(errList)
	plt.plot(do.timeList[do.nFirstClump:], pArr)
	plt.plot(do.timeList[do.nFirstClump:], pArr-errArr)
	plt.plot(do.timeList[do.nFirstClump:], pArr+errArr)
	plt.ylim(1.0, 3.0)
	pAvg, errAvg = get_p_avg(do, nStart=nStart, nEnd=nEnd)
	plt.ylabel(r'$p$')
	plt.xlabel(r'$t \Omega$')
	plt.title(r'$p_{end}=$' + str(pAvg) + r'$\pm$' + str(errAvg))



































#
