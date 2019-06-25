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
	def __init__(self, inArr, doPlan, splitterIds):
		print("initializing clump data structure...")
		self.nArr  = inArr[:,0]
		self.idArr = inArr[:,1]
		self.n0    = self.nArr[0]
		self.nLast = self.nArr[-1]
		self.massDict   = {}
		self.posDict    = {}
		self.splitter   = False
		nId = 0
		for n in self.nArr:
			for peak in doPlan.peakArrayList[n]:
				if peak[0] == self.idArr[nId]:
					self.massDict[n] = peak[2]
					self.posDict[n]  = peak[4:7]
			nId += 1
		if len(self.massDict.values()) != self.nArr.shape[0]:
			print("couldn't match this clump to a peak at all frames")
		if self.idArr[0] in splitterIds:
			print("marking this clump as a splitter")
			self.splitter = True

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
splitterIds = np.load(pathCT + "splitterIds.npy")
clumpObjList = []
ctFileNameList = os.listdir(pathCT)
for fileName in ctFileNameList:
	inArr = np.load(pathCT + fileName)
	clumpObjList.append(Clump(inArr, doPlan, splitterIds))

#for clump in clumpObjList:
	#plt.semilogy(clump.nArr, clump.massArr)
#plt.show(); plt.clf()

# scatter plots
'''
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
'''






# get initial mass spectrum
nSplitters=0; nStats=0;
mp = []
for clump in clumpObjList:
	m0 = clump.massDict[clump.n0]
	addToList = True
	# apply conditions to clumps to make the list that get their stats taken
	#if len(clump.massDict.keys())<10: addToList = False
	#if m0 < 1.e-4: addToList = False
	if clump.splitter: addToList = False; nSplitters+=1;
	if m0 > 1:         addToList = False; nSplitters+=1;
	if addToList: mp.append(m0); nStats+=1;

print("total clumps ever: " + str(len(clumpObjList)))
print("splitters: " + str(nSplitters))
print("doing stats on: " + str(nStats))







# do advanced stats on IMS

################################################################################
# make hist
#mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
#mp = np.asarray(mp); dndmp = np.asarray(dndmp);
mp1, ngtm = readerPlan.getCumMassHist2(mp)
minMass = np.amin(mp1); maxMass = np.amax(mp1);
nm = mp1.shape[0]
plt.figure(0)
################################################################################
nb = 1
# Single PL MLE
p_mle, err_mle = readerPlan.get_p_mle2(mp1)
# Fit with spl
means_spl, errsPlus_spl, errsMinus_spl, maxLike_spl = readerPlan.bootstrap(mp1, readerPlan.fit_spl, 1, nb=nb)
# Fit with stpl
means_stpl, errsPlus_stpl, errsMinus_stpl, maxLike_stpl = readerPlan.bootstrap(mp1, readerPlan.fit_stpl, 2, nb=nb)
# Fit with vtpl
means_vtpl, errsPlus_vtpl, errsMinus_vtpl, maxLike_vtpl =	readerPlan.bootstrap(mp1, readerPlan.fit_vtpl, 3, nb=nb)
# Fit with bcpl
means_bcpl, errsPlus_bcpl, errsMinus_bcpl, maxLike_bcpl =	readerPlan.bootstrap(mp1, readerPlan.fit_bcpl, 3, nb=nb)
# Fit with tpl
means_tpl, errsPlus_tpl, errsMinus_tpl, maxLike_tpl     = readerPlan.bootstrap(mp1, readerPlan.fit_tpl, 2, nb=nb)
# Fit with bpl
means_bpl, errsPlus_bpl, errsMinus_bpl, maxLike_bpl     = readerPlan.bootstrap(mp1, readerPlan.fit_bpl, 3, nb=nb)

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
print(mp1_min)
################################################################################
# cumulative hist plot
logMinMass = np.log10(minMass); logMaxMass = np.log10(maxMass);
#fakeMp1 = np.logspace(logMinMass-0.01 , logMaxMass+0.01, num=mp1.shape[0])
#fakeMp1 = np.logspace(logMinMass , logMaxMass+0.01, num=1000)
fakeMp1 = mp1
#fakeMp1[-1]+=0.01
#print(fakeMp1)
#print(fakeMp1.shape)
# BPL

print(means_bpl)

means_bpl = (-0.814, 1.3, 3.016)
ngtm_bpl = readerPlan.P_bpl(fakeMp1, means_bpl)

fac = 1.0
plt.loglog(mp1, fac*ngtm, color=(0,0,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, 1.01*nm*fac*ngtm_bpl, color=(1,0,0,1), label='BPL', linestyle='--', linewidth=1)
arg = np.argmin(np.absolute(minMass*np.exp(means_bpl[2])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_bpl[arg], color=(1,0,0,1), marker="s", linewidth=0, markersize=5)


# plot labels etc.
plt.xlabel(r'$M_p [M_G]$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(3.e-1, 1.e3)
plt.xlim(1e-3, 1.e0)
#plt.legend(prop={'size':10})
tools.saveAndClear(pathSave + "hist_cumulative_ct_jake.png", figNum=0)

np.savetxt("./mp_weak.csv", mp1)
np.savetxt("./ngtm_weak.csv", ngtm)














#
