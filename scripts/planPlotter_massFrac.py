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
pathSave = '../../plots/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
# set up data object and other params
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlanList = [readerPlan.DataPlan("../../data/prodRuns/run100/planOutput2/", nStart=220, nTot=275),
			  readerPlan.DataPlan("../../data/prodRuns/run103/planOutput2/", nStart=220, nTot=299),
			  readerPlan.DataPlan("../../data/prodRuns/run101/planOutput2/", nStart=300, nTot=390)]
nStartList = [220, 220, 300]
################################################################################
# plot mass frac in planetesimals over time
colorList = ['k', 'b', 'r', 'g']
labelList = ['control', r'$\alpha=10^{-4}$', r'$\alpha=10^{-3.5}$', r'$\alpha=10^{-3}$']
for n in range(4):
	try:
		doPlan = doPlanList[n]
		mFracList = []
		for item in doPlan.peakArrayList:
			try:    mNow = np.sum(item[:,2])
			except: mNow = 0.0
			mFracNow = mNow / doPlan.mParTot
			mFracList.append(mFracNow)
		plt.plot(doPlan.time, mFracList, color=colorList[n], linewidth=2, label=labelList[n])
		print(mFracList[-1])
	except:
		plt.plot(np.arange(100), np.zeros(100), color=colorList[n], linewidth=2, label=labelList[n])
	n+=1

plt.xlim(0, 19)
plt.legend(fontsize=15)
plt.ylabel(r'$M_{plan} / M_{par}$', fontsize=15)
plt.xlabel(r'$t \Omega$', fontsize=15)
plt.tight_layout()
plt.savefig(pathSave + "massFrac.png", bbox_inches='tight')
plt.close('all')




# plot nClumps over time
#axNum = 2
#ax[axNum].plot(doPlan1.time[nStart:], doPlan1.nClumpsList[nStart:], 'k', linewidth=2)
#ax[axNum].set_ylabel(r'$N_{clumps}$')
#ax[axNum].set_xlabel(r'$t \Omega$')
#ax[axNum].set_xlim(xLim1, doPlan1.tMax)
#ax[axNum].axvline(t_vline1, color=(0,0,1,1), linestyle='--')
#ax[axNum].axvline(t_vline2, color=(1,0,0,1), linestyle='--')










#
