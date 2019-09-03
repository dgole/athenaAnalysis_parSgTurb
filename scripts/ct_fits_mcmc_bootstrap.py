#!/usr/bin/python
import numpy as np
import time
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
################################################################################
# paths and CL args
pathBase  = str(sys.argv[1])
nStart    = int(sys.argv[2])
nStop     = int(sys.argv[3])
nb        = int(sys.argv[4])
#pathPlan  = pathBase + 'planOutput2/'
pathPlan  = pathBase
pathCT    = pathPlan + 'clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
pathSave  = pathBase + 'plots/clumpTracking_mcmc_bootstrap/'
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
		#print(fileName)
		#print(inArr)
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
print("total clumps ever: " + str(len(clumpObjList)))
print("splitters: " + str(nSplitters))
print("doing stats on: " + str(nStats))
print("average mass: " + str(np.mean(np.asarray(mp))))
print("\n")

################################################################################
# do advanced stats on IMS
# make hist
mp1, ngtm = readerPlan.getCumMassHist2(mp)
minMass = np.amin(mp1); maxMass = np.amax(mp1);
nm = mp1.shape[0]
plt.figure(0)

################################################################################

fitInfoList = [
			   pstats.fitInfo_spl,
			   pstats.fitInfo_stpl,
			   #pstats.fitInfo_tpl,
			   pstats.fitInfo_bcpl,
			   pstats.fitInfo_bpl,
			   #pstats.fitInfo_vtpl,
			   #pstats.fitInfo_tspl
			   ]

finalStats = []
for fitInfo in fitInfoList:
	fits = []
	paramsAll = np.zeros([nb, len(fitInfo.params0)])
	t0 = time.time()
	n  = 0
	while n < nb:
		if n%1==0:
			print("##########################################################################")
			str1 = "n = " + str(n) + " of " + str(nb)
			str2 = "average time per sample = " + str(np.round((time.time()-t0)/float(n+0.01),3))
			print(str1 + ", " + str2)
		sample          = np.random.choice(mp1, size=nm, replace=True)
		fits.append(pstats.Fit_Pipeline(sample, fitInfo, verbose=False, makePlots=False))
		try:
			print(fits[-1].fitInfo.name, fits[-1].params_opt)
			paramsAll[n] = fits[-1].params_opt
			n+=1
		except:
			print("fit failed, not advancing n")
	np.save(pathSave+fitInfo.name+"_"+str(nb)+".npy", paramsAll)
	finalStats    = np.zeros([len(fitInfo.params0), 3])
	finalStats[:,0] = np.percentile(paramsAll, 50, axis=0)
	finalStats[:,1] = np.percentile(paramsAll, 84, axis=0)
	finalStats[:,2] = np.percentile(paramsAll, 16, axis=0)
	np.savetxt(pathSave+fitInfo.name+"_"+str(nb)+".txt", finalStats)



################################################################################
