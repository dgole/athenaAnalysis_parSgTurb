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
import ctReader as readerCt
import plan_stats as pstats
plt.style.use('seaborn-dark-palette')
################################################################################
# paths and CL args
runId        = 100
pathBase = "../../data/prodRuns/run100/"
nStart     = 220
nStop      = 260
pathPlanList = [pathBase + 'planOutput_lowDensThresh/', pathBase + 'planOutput_medDensThresh/', pathBase + 'planOutput_highDensThresh/']
pathCTList   = [pathPlanList[i] + 'clumpTracking_' + str(nStart) + '_' + str(nStop) + '/' for i in range(len(pathPlanList)) ]
pathSave   = pathBase + 'plots/clumpTracking_mcmc_compareThresholds/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
plotFrac   = np.power(10,1.5)
fac        = np.power(plotFrac,(len(pathPlanList)-1))
color1List = [(1,0,0,1),(0,0,0,1),(0,0,1,1)]
color2List = [(1,0,0,0.2),(0,0,0,0.2),(0,0,1,0.2)]
strList    = []
str1       = (r"$\alpha_1 = -2.36^{+0.32}_{-0.35}$" + "\n" +
			  r"$\alpha_2 =  1.04^{+0.06}_{-0.05}$" + "\n" +
			  r"$M_{br}   =  0.0081^{+0.0005}_{-0.0005}$" )
strList.append(str1)
str1       = (r"$\alpha_1 = -2.06^{+0.24}_{-0.26}$" + "\n" +
			  r"$\alpha_2 =  1.71^{+0.17}_{-0.16}$" + "\n" +
			  r"$M_{br}   = 0.0397^{+0.0028}_{-0.0021}$" )
strList.append(str1)
str1       = (r"$\alpha_1 = -1.96^{+0.26}_{-0.30}$" + "\n" +
			  r"$\alpha_2 =  1.68^{+0.22}_{-0.18}$" + "\n" +
			  r"$M_{br}   = 0.0636^{+0.0049}_{-0.0054}$" )
strList.append(str1)

labelList  = ["low threshold", "normal threshold", "high threshold"]
# read peak files
for i in range(len(pathPlanList)):
	pathPlan = pathPlanList[i]
	pathCT   = pathCTList[i]
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
	print("min mass: " + str(np.amin(mp)))
	################################################################################
	# do advanced stats on IMS
	# make hist
	mp1, ngtm = readerPlan.getCumMassHist2(mp)
	minMass = np.amin(mp1); maxMass = np.amax(mp1);
	nm = mp1.shape[0]
	plt.figure(num=0, figsize=(10,10))
	################################################################################
	fitInfo = pstats.fitInfo_bpl
	fit     = pstats.Fit_Pipeline(mp1, fitInfo)
	################################################################################
	# cumulative hist plot
	plt.figure(0)
	logMinMass = np.log10(minMass)
	logMaxMass = np.log10(maxMass)
	fakeMp     = np.logspace(logMinMass, logMaxMass+3, num=1000)
	arg1       = np.argmin(np.absolute(fakeMp-minMass))
	arg2       = np.argmin(np.absolute(fakeMp-maxMass))+1
	print('making histogram for ' + fitInfo.name)
	p       = fitInfo.pFunc(fit.params_mcmc[0], fakeMp)
	P       = pstats.p_to_P(fakeMp, p)
	plt.loglog(fakeMp[arg1:arg2], fac*nm*P[arg1:arg2],
			   color=color1List[i], linewidth=1, linestyle='--',
			   label=labelList[i])
	plt.loglog(mp1, fac*ngtm,
			   color=color1List[i], marker=".", linewidth=0, markersize=3)
	plt.loglog(mp1, fac*ngtm,
			   color=color2List[i], linewidth=8)
	#plt.loglog(mp1, fac*minMass*np.exp(fit.params_mcmc[0,2]),
			   #color=color1List[i], marker=".", linewidth=0, markersize=10)
	plt.text(1.0, fac*ngtm[-1], strList[i])
	fac /= plotFrac

# plot labels etc.
plt.xlabel(r'$M_p [M_G]$', fontsize=14)
plt.ylabel(r'$N(>M_p)$', fontsize=14)
plt.ylim(1.e-2, 1000*np.power(plotFrac,(len(pathPlanList)-1)))
plt.xlim(2e-3, 8.e0)
plt.legend(prop={'size':13}, loc='lower left')
tools.saveAndClear(pathSave + "hist_cumulative_ct_compare.png", figNum=0)
















#
