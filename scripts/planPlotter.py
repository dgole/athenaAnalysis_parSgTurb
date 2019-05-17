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
pathBase = str(sys.argv[1])
np1      = int(sys.argv[2])
np2      = int(sys.argv[3])
nTot     = int(sys.argv[4])
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/plan/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
# set up data object and other params
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan = readerPlan.DataPlan(pathPlan, nStart=200, nTot=nTot, nPar=128*128*128)
tp1    = np1*doPlan.dt
tp2    = np2*doPlan.dt


################################################################################
# MLE
################################################################################

# calculate the PL slope with the MLE at any times with enough planetesimals
pList=[]; errList=[]; tList=[]
for n in range(0, nTot):
	try:
		mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
		if len(mp)>4:
			p, err = readerPlan.get_p(doPlan, n)
			pList.append(p)
			errList.append(err)
			tList.append(doPlan.timeList[n])
	except:
		pass
pArr = np.asarray(pList); errArr = np.asarray(errList)

# calculate average p value from MLE over desired time range
pList2=[]; errList2=[];
for n in range(len(tList)):
	if np1 < int(tList[n]/doPlan.dt) < np2:
		pList2.append(pArr[n])
		errList2.append(errArr[n])
pArr2 = np.asarray(pList2); errArr2 = np.asarray(errList2)
p = np.median(pArr2); err1 = np.std(pArr2); err2 = np.median(errList2)

# plot p_mle value over time
plt.plot(tList, pArr, 'k')
plt.fill_between(tList, pArr-errArr, pArr+errArr, color='gray', alpha=0.5)
plt.ylim(1.0, 2.5)
plt.ylabel(r'$p$')
plt.xlabel(r'$t \Omega$')
plt.plot((tp1, tp2), (p,     p    ), color='r', linestyle='--')
plt.plot((tp1, tp2), (p+err2, p+err2), color='b', linestyle='--')
plt.plot((tp1, tp2), (p-err2, p-err2), color='b', linestyle='--')
plt.axhline(1.7, color=(0,0,0,0.2), linestyle='--')
plt.axhline(1.4, color=(0,0,0,0.2), linestyle='--')
plt.title(r'$p=$'  + str(np.round(p,   2)) +
          r'$\pm$' + str(np.round(err2,2)))
tools.saveAndClear(pathSave + "p_vs_t.png", figNum=0)


################################################################################
# Fit the histogram
################################################################################

# calculate p by fitting the diff hist
pFitList=[]; cFitList=[];
for n in range(0, nTot):
	try:
		mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
		if len(mp)>4:
			mp = np.asarray(mp); dndmp = np.asarray(dndmp);
			x = np.log10(mp) ; y = np.log10(dndmp)
			coefs = poly.polyfit(x, y, 1)
			this_p = -coefs[1]; this_c = np.power(10,coefs[0])
			pFitList.append(this_p); cFitList.append(this_c);
			#if n in np.arange(320,450,10):
				#plt.loglog(mp, dndmp, 'ko')
				#plt.loglog(mp, this_c*np.power(mp, -this_p), 'b-')
				#plt.show(); plt.clf();
	except:
		pass

# calculate average p value from fit over desired time range
pFitList2=[]; cFitList2=[];
for n in range(len(tList)):
	if np1 < int(tList[n]/doPlan.dt) < np2:
		pFitList2.append(pFitList[n])
		cFitList2.append(cFitList[n])
pFitArr2=np.asarray(pFitList2); cFitArr2=np.asarray(cFitList2);
pFit = np.median(pFitArr2); pFitErr1 = np.std(pArr2);
cFit = np.median(cFitArr2);

# plot pFit value over time
plt.figure(0)
plt.plot(tList, pFitList, 'k')
plt.ylim(0.0, 2.5)
plt.ylabel(r'$p$')
plt.xlabel(r'$t \Omega$')
plt.plot((tp1, tp2), (pFit,     pFit    ), color='r', linestyle='--')
plt.plot((tp1, tp2), (pFit+pFitErr1, pFit+pFitErr1), color='b', linestyle='--')
plt.plot((tp1, tp2), (pFit-pFitErr1, pFit-pFitErr1), color='b', linestyle='--')
plt.axhline(1.7, color=(0,0,0,0.2), linestyle='--')
plt.axhline(1.4, color=(0,0,0,0.2), linestyle='--')
plt.title(r'$p=$'  + str(np.round(pFit,     2)) +
          r'$\pm$' + str(np.round(pFitErr1, 2)))
tools.saveAndClear(pathSave + "pFit_vs_t.png", figNum=0)


################################################################################
# plot diff hist with both models
################################################################################

# gather all mp, dndmp fromm relevant time period and average
mp_interp         = np.logspace(-2.5, -0.5, num=1000)
dndmp_interp_list = []
for n in range(np1, np2):
	mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
	mp = np.asarray(mp); dndmp = np.asarray(dndmp);
	this_interp = np.interp(mp_interp, mp, dndmp)
	dndmp_interp_list.append(this_interp)
	plt.loglog(mp, dndmp, 'ko', ms=2)
dndmp_interp = np.mean(np.asarray(dndmp_interp_list), axis=0)
plt.loglog(mp_interp, dndmp_interp, 'b-')

# plot MLE model slope
mp = np.arange(1.e-3, 1, 1.e-3)
for i in np.arange(-5,5,0.5):
	preFactor = np.power(10,i)
	model = preFactor*np.power(mp, -p)
	plt.loglog(mp, model, color=(0,0,0,0.2), linestyle='--', linewidth=1)

# plot average fit to histogram
def modelFunc(mp1):
	return cFit * np.power(mp1, -pFit)
model = modelFunc(mp)
plt.loglog(mp, model, color=(0,0,0,0.2), linestyle='-', linewidth=1)

# plot labels etc.
plt.xlabel(r'$M_p$')
plt.ylabel(r'$dN/dM_p$')
plt.ylim(0.8*np.amin(dndmp_interp), 1.2*np.amax(dndmp_interp))
plt.xlim(0.8*np.amin(mp_interp),    1.2*np.amax(mp_interp))
tools.saveAndClear(pathSave + "hist_differential.png", figNum=0)


################################################################################
# Other useful plots
################################################################################

# plot nClumps over time
plt.plot(doPlan.time, doPlan.nClumpsList, 'k', linewidth=2)
plt.ylabel(r'$N_{clumps}$')
plt.xlabel(r'$t \Omega$')
tools.saveAndClear(pathSave + "nClumps.png", figNum=0)

# plot mass frac in planetesimals over time
mFracList = []
for item in doPlan.peakArrayList:
	try:
		mNow = np.sum(item[:,2])
	except:
		mNow = 0.0
	mFracNow = mNow / doPlan.mParTot
	mFracList.append(mFracNow)
plt.plot(doPlan.timeList, mFracList, 'k', linewidth=2)
plt.xlabel(r'$t \Omega$')
plt.ylabel(r'$M_{plan} / M_{par}$')
tools.saveAndClear(pathSave + "massFrac.png", figNum=0)














#
