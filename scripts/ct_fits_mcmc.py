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
runId        = int(sys.argv[1])
densThresh   = None
if len(sys.argv) > 2: densThresh = str(sys.argv[2])
pathBaseDict = {100: "../../data/prodRuns/run100/",
				101: "../../data/prodRuns/run101/",
				102: "../../data/prodRuns/run102/",
				103: "../../data/prodRuns/run103/"
				}
nDict        = {100: [220, 275],
				101: [300, 390],
				102: [000, 000],
				103: [220, 290],
				}
pathBase   = pathBaseDict[runId]
nStart     = nDict[runId][0]
nStop      = nDict[runId][1]
pathPlan   = pathBase + 'planOutput2/'
pathCT     = pathBase + 'planOutput2/clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
pathSave   = pathBase + 'plots/clumpTracking_mcmc/'
if densThresh is not None:
	pathPlan   = pathBase + 'planOutput_' + densThresh + "/"
	pathCT     = pathPlan + 'clumpTracking_' + str(nStart) + '_' + str(nStop) + '/'
	pathSave   = pathBase + 'plots/clumpTracking_mcmc_' + densThresh + "/"
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
file = open(pathSave+"report.txt", 'w')
file.write("total clumps ever: " + str(len(clumpObjList)) + "\n")
file.write("splitters: " + str(nSplitters) + "\n")
file.write("doing stats on: " + str(nStats) + "\n")
file.write("average mass: " + str(np.mean(np.asarray(mp))) + "\n")
file.write("min mass: " + str(np.amin(mp)) + "\n")
file.write("total mass: " + str(np.sum(np.asarray(mp))))

################################################################################
# do advanced stats on IMS
# make hist
mp1, ngtm = readerPlan.getCumMassHist2(mp)
minMass = np.amin(mp1); maxMass = np.amax(mp1);
nm = mp1.shape[0]
plt.figure(num=0, figsize=(6,10))

################################################################################

fitInfoList = [pstats.fitInfo_spl,
			   pstats.fitInfo_stpl,
			   #pstats.fitInfo_tpl, # not working
			   pstats.fitInfo_bcpl,
			   pstats.fitInfo_bpl,
			   #pstats.fitInfo_vtpl, # not working
			   pstats.fitInfo_tspl
			   ]

fits = []
for fitInfo in fitInfoList:
	fits.append(pstats.Fit_Pipeline(mp1, fitInfo))
	fits[-1].cp1.savefig(pathSave + "cp1_"+ fitInfo.name + ".png")
	#fits[-1].cp2.savefig(pathSave + "cp2_"+ fitInfo.name + ".png")


################################################################################
# information criteria analysis

bicList = [fit.bic for fit in fits]
aicList = [fit.aic for fit in fits]
bic_min = min(bicList)
aic_min = min(aicList)
for fit in fits:
	fit.dbic = fit.bic - bic_min
	fit.daic = fit.aic - aic_min

file.write("\n\n#########################\n")
file.write("FINAL REPORT\n")
for fit in fits:
	file.write("#########################\n")
	file.write(fit.fitInfo.name + "\n")
	for i in range(len(fit.fitInfo.paramNames)):
		file.write(fit.fitInfo.paramNames[i].ljust(10) + " = ")
		file.write(str(fit.params_opt[i]) + "\n")
	file.write("ICS: ")
	file.write(str(fit.dbic) + ", ")
	file.write(str(fit.daic) + "\n")
file.close()

################################################################################
# cumulative hist plot
plt.figure(0)
logMinMass = np.log10(minMass)
logMaxMass = np.log10(maxMass)
fakeMp     = np.logspace(logMinMass, logMaxMass+3, num=1000)
arg1       = np.argmin(np.absolute(fakeMp-minMass))
arg2       = np.argmin(np.absolute(fakeMp-maxMass))+1
plotFrac   = np.power(10,1.5)
fac        = np.power(plotFrac,(len(fits)-1))
for i in range(len(fits)):
	fit     = fits[i]
	fitInfo = fitInfoList[i]
	print('making histogram for ' + fitInfo.name)
	p       = fitInfo.pFunc(fit.params_opt, fakeMp)
	P       = pstats.p_to_P(fakeMp, p)
	plt.loglog(fakeMp[arg1:arg2], fac*nm*P[arg1:arg2],
			   color=fitInfo.color1, linewidth=1, linestyle='--',
			   label=fitInfo.name.upper())
	plt.loglog(mp1, fac*ngtm,
			   color=fitInfo.color1, marker=".", linewidth=0, markersize=3)
	plt.loglog(mp1, fac*ngtm,
			   color=fitInfo.color2, linewidth=8)
	fac /= plotFrac


# plot labels etc.
plt.xlabel(r'$M_p [M_G]$', fontsize=14)
plt.ylabel(r'$N(>M_p)$', fontsize=14)
plt.ylim(1.e-2, 1000*np.power(plotFrac,(len(fits)-1)))
plt.xlim(2e-3, 8.e-1)
plt.legend(prop={'size':13}, loc='lower left')
tools.saveAndClear(pathSave + "hist_cumulative_ct.png", figNum=0)
















#
