f#!/usr/bin/python
import numpy as np
#import matplotlib as m
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaTools as tools
import planOutputReader as readerPlan
import ctReader as readerCt
import plan_stats as pstats
plt.style.use('seaborn-dark-palette')
################################################################################
# paths and CL args
pathBase  = str(sys.argv[1])
nStart    = int(sys.argv[2])
nStop     = int(sys.argv[3])
pathPlan  = pathBase + 'planOutput2/'
pathCT    = pathBase + 'planOutput2/clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
pathSave  = pathBase + 'plots/clumpTracking_mcmc_snaps/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
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
		clumpObjList.append(readerCt.Clump(inArr, doPlan, splitterIds, splitterData))

################################################################################
fitInfo  = pstats.fitInfo_bpl
fits     = []
n1List   = [220, 230, 240, 250, 260]
n2List   = [230, 240, 250, 260, 270]
minMass  = 1000.0
maxMass  = 0.0
mp1List  = []
ngtmList = []
for i in range(len(n1List)):
	nSplitters=0; nStats=0;
	mp = []
	for clump in clumpObjList:
		m0 = clump.massDict[clump.n0]
		addToList = True
		# apply conditions to clumps to make the list that get their stats taken
		#if len(clump.massDict.keys())<10: addToList = False
		#if m0 < 1.e-4: addToList = False
		#if m0 > 1.0: addToList = False
		if clump.n0 < n1List[i] or clump.n0 > n2List[i]: addToList = False;
		if clump.splitter: addToList = False; nSplitters+=1;
		if addToList: mp.append(m0); nStats+=1;
	print("total clumps ever: " + str(len(clumpObjList)))
	print("splitters: " + str(nSplitters))
	print("doing stats on: " + str(nStats))
	print("average mass: " + str(np.mean(np.asarray(mp))))
	print("min mass: " + str(np.amin(mp)))
	# make hist
	mp1, ngtm = readerPlan.getCumMassHist2(mp)
	fits.append(pstats.Fit_Pipeline(mp1, fitInfo))
	mp1List.append(mp1)
	ngtmList.append(ngtm)
	if np.amin(mp1) < minMass: minMass = np.amin(mp1);
	if np.amax(mp1) > maxMass: maxMass = np.amax(mp1);
	print(fits[-1].params_mcmc)
################################################################################
# cumulative hist plot
plt.figure(num=0, figsize=(6,10))
plotFrac   = np.power(10,1.5)
fac        = np.power(plotFrac,(len(fits)-1))
for i in range(len(fits)):
	fit     = fits[i]
	label   = labelList[i]
	nm      = len(mp1List[i])
	mp1     = mp1List[i]
	ngtm    = ngtmList[i]
	nm      = len(mp1)

	logMinMass = np.log10(np.amin(mp1))
	logMaxMass = np.log10(np.amax(mp1))
	fakeMp     = np.logspace(logMinMass, logMaxMass+3, num=1000)
	arg1       = np.argmin(np.absolute(fakeMp-minMass))
	arg2       = np.argmin(np.absolute(fakeMp-maxMass))+1

	p       = fitInfo.pFunc(fit.params_mcmc[0], fakeMp)
	P       = pstats.p_to_P(fakeMp, p)

	plt.loglog(fakeMp[arg1:arg2], fac*nm*P[arg1:arg2],
			   color=(0,0,0,1), linewidth=1, linestyle='--',
			   label=label)
	plt.loglog(mp1, fac*ngtm,
			   color=(0,0,0,1), marker=".", linewidth=0, markersize=2)
	plt.loglog(mp1, fac*ngtm,
			   color=(0,0,0,0.2), linewidth=5)
	plt.text(1.0, fac*ngtm[-1], str(np.round(np.transpose(fit.params_mcmc),3)))
	fac /= plotFrac

# plot labels etc.
plt.xlabel(r'$M_p [M_G]$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(6.e-1, 500*np.power(plotFrac,(len(fits)-1)))
plt.xlim(1e-3, 3.0)
#plt.text()
tools.saveAndClear(pathSave + "hist_cumulative_ct_snaps.png", figNum=0)
















#
