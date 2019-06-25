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
nb        = int(sys.argv[4])
pathPlan  = pathBase + 'planOutput2/'
pathCT    = pathBase + 'planOutput2/clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
pathSave  = pathBase + 'plots/clumpTracking3/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################

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

################################################################################

# read peak files
doPlan = readerPlan.DataPlan(pathPlan, nStart=nStart, nTot=nStop, nPar=512*512*512)

# read splitter ids and data
splitterIds  = np.load(pathCT + "splitterIds.npy", allow_pickle=True)
splitterData  = np.load(pathCT + "splitterData.npy", allow_pickle=True)

clumpObjList = []
ctFileNameList = os.listdir(pathCT)
for fileName in ctFileNameList:
	if fileName[0:8] != "splitter":
		inArr = np.load(pathCT + fileName, allow_pickle=True)
		#print(fileName)
		#print(inArr)
		clumpObjList.append(Clump(inArr, doPlan, splitterIds, splitterData))



################################################################################
# get initial mass spectrum
nSplitters=0; nStats=0;
mp = []
for clump in clumpObjList:
	m0 = clump.massDict[clump.n0]
	addToList = True
	# apply conditions to clumps to make the list that get their stats taken
	#if len(clump.massDict.keys())<10: addToList = False
	#if m0 < 1.e-4: addToList = False
	#if m0 > 1.0: addToList = False
	if clump.splitter: addToList = False; nSplitters+=1;
	if addToList: mp.append(m0); nStats+=1;
print("total clumps ever: " + str(len(clumpObjList)))
print("splitters: " + str(nSplitters))
print("doing stats on: " + str(nStats))
print("average mass: " + str(np.mean(np.asarray(mp))))

# get split fraction distribution
splitFracList = []
massOfSplittersList = []
persistenceList = []
for clump in clumpObjList:
	if clump.splitFrac > -1.0:
		splitFracList.append(clump.splitFrac)
		massOfSplittersList.append(clump.massDict[clump.n0])
		persistenceList.append(clump.persistence)
plt.figure(0)
plt.hist(np.asarray(splitFracList), bins=20)
plt.xlim(0,1.0)
plt.ylabel("Count")
plt.xlabel("Split Fraction")
tools.saveAndClear(pathSave + "hist_splitFracs_ct.png", figNum=0)

plt.figure(0)
plt.semilogx(massOfSplittersList, splitFracList, 'ko', markersize=2)
#plt.xlim(0,1.0)
plt.xlabel("Mp (initial)")
plt.ylabel("Split Frac")
#plt.axhline(0.9, color='gray', linestyle="--")
plt.axhline(0.5, color='gray', linestyle="--")
#plt.axhline(0.2, color='gray', linestyle="--")
plt.axvline(1.e-2, color='gray', linestyle="--")
countBL = np.sum(
		  np.where(np.asarray(splitFracList)<0.5, 1, 0) *
		  np.where(np.asarray(massOfSplittersList)<1.e-2,
		  1, 0))
countTL = np.sum(
		  np.where(np.asarray(splitFracList)>0.5, 1, 0) *
		  np.where(np.asarray(massOfSplittersList)<1.e-2,
		  1, 0))
countBR = np.sum(
		  np.where(np.asarray(splitFracList)<0.5, 1, 0) *
		  np.where(np.asarray(massOfSplittersList)>1.e-2,
		  1, 0))
countTR = np.sum(
		  np.where(np.asarray(splitFracList)>0.5, 1, 0) *
		  np.where(np.asarray(massOfSplittersList)>1.e-2,
		  1, 0))
plt.text(5.e-3, 0.25, str(countBL))
plt.text(5.e-3, 0.75, str(countTL))
plt.text(1.e-1, 0.25, str(countBR))
plt.text(1.e-1, 0.75, str(countTR))

tools.saveAndClear(pathSave + "splitFrac_vs_mass.png", figNum=0)

plt.figure(0)
plt.semilogx(massOfSplittersList, persistenceList, 'ko', markersize=2)
#plt.xlim(0,1.0)
plt.xlabel("Mp (initial)")
plt.ylabel("Persistence")
tools.saveAndClear(pathSave + "persistence_vs_mass.png", figNum=0)



################################################################################
# do advanced stats on IMS
# make hist
mp1, ngtm = readerPlan.getCumMassHist2(mp)
minMass = np.amin(mp1); maxMass = np.amax(mp1);
nm = mp1.shape[0]
plt.figure(0)






################################################################################
# diff hist
mp, dndmp = readerPlan.getDiffMassHist2(mp)
mp = np.asarray(mp); dndmp = np.asarray(dndmp);
min = 1.e20
preFacArr = np.arange(1.0,200.0,0.1)
# Single PL MLE
p_mle, err_mle = readerPlan.get_p_mle2(mp1)
for i in range(preFacArr.shape[0]):
	preFac   = preFacArr[i]
	model    = preFac * np.power(mp, -p_mle)
	logDndmp = np.log10(dndmp)
	logModel = np.log10(model)
	diff     = np.sum(np.square(logDndmp - logModel))
	if diff < min:
		min  = diff
		bestPreFac = preFac
		#print("new min found")
		#print(min)
		#print(bestPreFac)

plt.figure(0)
plt.loglog(mp, dndmp, 'ko', ms=2)
mp2 = np.logspace(-10.0, 10, num=500)
model = bestPreFac * np.power(mp2, -p_mle)
plt.loglog(mp2, model, linestyle='--', linewidth=1, color='k')
# plot labels etc.
plt.xlabel(r'$M_p$')
plt.ylabel(r'$dN/dM_p$')
plt.ylim(3.e0,  1.e5)
plt.xlim(1e-3, 1.e0)
tools.saveAndClear(pathSave + "hist_differential_oneFrameExample_" + ".png", figNum=0)















################################################################################

pathSave1 = pathSave
# Single PL MLE
p_mle, err_mle = readerPlan.get_p_mle2(mp1)
# Fit with spl
means_spl, errsPlus_spl, errsMinus_spl, maxLike_spl = readerPlan.bootstrap(mp1, readerPlan.fit_spl, 1, nb=nb, pathSave=pathSave1)
# Fit with stpl
means_stpl, errsPlus_stpl, errsMinus_stpl, maxLike_stpl = readerPlan.bootstrap(mp1, readerPlan.fit_stpl, 2, nb=nb, pathSave=pathSave1)
# Fit with vtpl
means_vtpl, errsPlus_vtpl, errsMinus_vtpl, maxLike_vtpl =	readerPlan.bootstrap(mp1, readerPlan.fit_vtpl, 3, nb=nb, pathSave=pathSave1)
# Fit with bcpl
means_bcpl, errsPlus_bcpl, errsMinus_bcpl, maxLike_bcpl =	readerPlan.bootstrap(mp1, readerPlan.fit_bcpl, 3, nb=nb, pathSave=pathSave1)
# Fit with tpl
means_tpl, errsPlus_tpl, errsMinus_tpl, maxLike_tpl     = readerPlan.bootstrap(mp1, readerPlan.fit_tpl, 2, nb=nb, pathSave=pathSave1)
# Fit with bpl
means_bpl, errsPlus_bpl, errsMinus_bpl, maxLike_bpl     = readerPlan.bootstrap(mp1, readerPlan.fit_bpl, 3, nb=nb, pathSave=pathSave1)

################################################################################
# information criteria analysis stuff

N = mp1.shape[0]

K = 1
bic_spl = readerPlan.BIC(K, N, maxLike_spl)
aic_spl = readerPlan.AIC(K, N, maxLike_spl)
K = 2
bic_stpl = readerPlan.BIC(K, N, maxLike_stpl)
aic_stpl = readerPlan.AIC(K, N, maxLike_stpl)
K = 3
bic_vtpl = readerPlan.BIC(K, N, maxLike_vtpl)
aic_vtpl = readerPlan.AIC(K, N, maxLike_vtpl)
K = 3
bic_bcpl = readerPlan.BIC(K, N, maxLike_bcpl)
aic_bcpl = readerPlan.AIC(K, N, maxLike_bcpl)
K = 2
bic_tpl = readerPlan.BIC(K, N, maxLike_tpl)
aic_tpl = readerPlan.AIC(K, N, maxLike_tpl)
K = 3
bic_bpl = readerPlan.BIC(K, N, maxLike_bpl)
aic_bpl = readerPlan.AIC(K, N, maxLike_bpl)

bic_min = min(bic_spl, bic_stpl, bic_vtpl, bic_bcpl, bic_tpl, bic_bpl)
aic_min = min(aic_spl, aic_stpl, aic_vtpl, aic_bcpl, aic_tpl, aic_bpl)

dbic_spl  = bic_spl  - bic_min
daic_spl  = aic_spl  - aic_min
dbic_stpl = bic_stpl - bic_min
daic_stpl = aic_stpl - aic_min
dbic_vtpl = bic_vtpl - bic_min
daic_vtpl = aic_vtpl - aic_min
dbic_bcpl = bic_bcpl - bic_min
daic_bcpl = aic_bcpl - aic_min
dbic_tpl  = bic_tpl  - bic_min
daic_tpl  = aic_tpl  - aic_min
dbic_bpl  = bic_bpl  - bic_min
daic_bpl  = aic_bpl  - aic_min

################################################################################
# print out all fits / params / errors in a nice format
#print("########################################################")
#print("Nclumps: " + str(doPlan.nClumpsList[n]))
mp1_min = np.amin(mp1)
print("########################################################")
print("Single PL MLE Slope")
print(str(np.round(p_mle,  3)) + " pm "
	+ str(np.round(err_mle,3)))
readerPlan.reportParams("SPL", ["alpha"],
		   means_spl, errsPlus_spl, errsMinus_spl,
		   maxLike_spl, dbic_spl, daic_spl,
		   mp1_min)
readerPlan.reportParams("STPL", ["alpha", "x_exp"],
		   means_stpl, errsPlus_stpl, errsMinus_stpl,
		   maxLike_stpl, dbic_stpl, daic_stpl,
		   mp1_min)
readerPlan.reportParams("VTPL", ["alpha", "beta", "x_exp"],
		   means_vtpl, errsPlus_vtpl, errsMinus_vtpl,
		   maxLike_vtpl, dbic_vtpl, daic_vtpl,
		   mp1_min)
readerPlan.reportParams("BCPL", ["alpha1", "alpha2", "x_br"],
		   means_bcpl, errsPlus_bcpl, errsMinus_bcpl,
		   maxLike_bcpl, dbic_bcpl, daic_bcpl,
		   mp1_min)
readerPlan.reportParams("TPL", ["alpha", "x_tr"],
		   means_tpl, errsPlus_tpl, errsMinus_tpl,
		   maxLike_tpl, dbic_tpl, daic_tpl,
		   mp1_min)
readerPlan.reportParams("BPL", ["alpha1", "alpha2", "x_br"],
		   means_bpl, errsPlus_bpl, errsMinus_bpl,
		   maxLike_bpl, dbic_bpl, daic_bpl,
		   mp1_min)
print("########################################################")
################################################################################
# cumulative hist plot
logMinMass = np.log10(minMass); logMaxMass = np.log10(maxMass);
#fakeMp1 = np.logspace(logMinMass-0.01 , logMaxMass+0.01, num=mp1.shape[0])
fakeMp1 = np.logspace(logMinMass , logMaxMass, num=1000)
#print(fakeMp1)
#print(fakeMp1.shape)
# SPL
ngtm_spl = readerPlan.P_spl(fakeMp1, means_spl)
fac = 1.e5
plt.loglog(mp1, fac*ngtm, color=(0,0,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_spl, color=(0,0,0,1), label='SPL', linestyle='--')
# STPL
ngtm_stpl = readerPlan.P_stpl(fakeMp1, means_stpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(1,0,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(1,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_stpl, color=(1,0,0,1), label='STPL', linestyle='--')
arg = np.argmin(np.absolute(minMass*np.exp(means_stpl[1])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_stpl[arg], color=(1,0,0,1), marker="s", linewidth=0, markersize=7)
# VTPL
ngtm_vtpl = readerPlan.P_vtpl(fakeMp1, means_vtpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(0,1,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,1,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_vtpl, color=(0,1,0,1), label='VTPL', linestyle='--')
#arg = np.argmin(np.absolute(minMass*np.exp(means_vtpl[2])-fakeMp1))
#plt.loglog(fakeMp1[arg], nm*fac*ngtm_vtpl[arg], color=(0,1,0,1), marker="s", linewidth=0, markersize=7)
# BCPL
ngtm_bcpl = readerPlan.P_bcpl(fakeMp1, means_bcpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(0,0,1,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,0,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_bcpl, color=(0,0,1,1.0), label='BCPL', linestyle='--')
arg = np.argmin(np.absolute(minMass*np.exp(means_bcpl[2])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_stpl[arg], color=(0,0,1,1), marker="s", linewidth=0, markersize=7)
# TPL
ngtm_tpl = readerPlan.P_tpl(fakeMp1, means_tpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(1,0,1,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(1,0,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_tpl, color=(1,0,1,1), label='TPL', linestyle='--')
# BPL
ngtm_bpl = readerPlan.P_bpl(fakeMp1, means_bpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(0,1,1,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,1,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_bpl, color=(0,1,1,1), label='BPL', linestyle='--')
arg = np.argmin(np.absolute(minMass*np.exp(means_bpl[2])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_bpl[arg], color=(0,1,1,1), marker="s", linewidth=0, markersize=7)
# plot labels etc.
plt.xlabel(r'$M_p [M_G]$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(6.e-1, 1.e8)
plt.xlim(1e-3, 1.e1)
plt.legend(prop={'size':10}, loc='upper right')
tools.saveAndClear(pathSave + "hist_cumulative_ct.png", figNum=0)
















#plt.figure(0)
#for clump in clumpObjList:
	#plt.semilogy(clump.nArr, clump.massArr)
#plt.show(); plt.clf()
'''
# scatter plots
for n in range(nStart, nStop):
	plt.figure(0)
	plt.figure(num=0, figsize=(3,2))
	xs     = get_xs(clumpObjList, n)
	ys     = get_ys(clumpObjList, n)
	colors = get_colors(clumpObjList, n)
	masses = get_currentMasses(clumpObjList, n)
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	plt.scatter(xs, ys, s=sizes, color=colors)
	for i in range(len(clumpObjList)):
		plt.text(xs[i]+0.001, ys[i]+0.001, str(i), fontsize=5)
	plt.xlabel(r'$r/h$')
	plt.ylabel(r'y/h')
	plt.ylim(-0.1, 0.1)
	plt.xlim(-0.1, 0.1)
	plt.tight_layout()
	tools.saveAndClear(pathSave + "scatter_" + str(n) + ".png", figNum=0)
	plt.close('all')
'''










#
