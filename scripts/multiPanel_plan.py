#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
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
################################################################################
# paths
pathBase = str(sys.argv[1])
G        = float(sys.argv[2])
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/plan/anim1/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan = readerPlan.DataPlan(pathPlan, G=G)
if len(sys.argv)>3: nStart = int(sys.argv[3])
else:               nStart = 0
################################################################################
pList   = []
errList = []
for n in range(0, doPlan.nt):
	p, err = readerPlan.get_p(doPlan, n)
	pList.append(p)
	errList.append(err)
pArr   = np.asarray(pList)
errArr = np.asarray(errList)
################################################################################
def makeAnimFrame(self, n):
	print('saving anim frame for n = ' + str(n))
	fig = plt.figure(figsize=(8.03,10.5), dpi=80)
	ax = []
	ax.append(plt.subplot2grid((7, 2), (0, 0), rowspan=2))
	ax.append(plt.subplot2grid((7, 2), (0, 1), rowspan=2))
	ax.append(plt.subplot2grid((7, 2), (2, 0), rowspan=2))
	ax.append(plt.subplot2grid((7, 2), (2, 1), rowspan=2))
	ax.append(plt.subplot2grid((7, 2), (4, 0), colspan=2))
	ax.append(plt.subplot2grid((7, 2), (5, 0), colspan=2))
	ax.append(plt.subplot2grid((7, 2), (6, 0), colspan=2))

	tPlot = doPlan.time

	# Nclumps
	axNum = 6
	ax[axNum].plot(tPlot[n:], doPlan.nClumpsList[n:], 'gray', linewidth=1)
	ax[axNum].plot(tPlot[:n], doPlan.nClumpsList[:n], 'k', linewidth=2)
	ax[axNum].plot(tPlot[n], doPlan.nClumpsList[n], 'ro', markersize=5)
	ax[axNum].set_ylabel(r'$N_{clumps}$')
	ax[axNum].set_xlabel(r'$t \Omega$')
	ax[axNum].set_xlim(0, doPlan.timeList[-1])

	# mass frac
	axNum = 5
	mFracList = []
	for n1 in range(0, len(doPlan.peakArrayList)):
		mNow = np.sum(doPlan.peakArrayList[n1][:,2])
		mFracNow = mNow / doPlan.mParTot
		mFracList.append(mFracNow)
	ax[axNum].plot(doPlan.timeList[n:], mFracList[n:], 'gray', linewidth=1)
	ax[axNum].plot(doPlan.timeList[:n], mFracList[:n], 'k', linewidth=2)
	ax[axNum].plot(doPlan.timeList[n], mFracList[n], 'ro', markersize=5)
	#ax[axNum].get_xaxis().set_visible(False)
	#ax[axNum].set_ylim(-0.05, 1.02)
	ax[axNum].set_ylabel(r'$M_{plan} / M_{par}$')
	ax[axNum].set_xlim(0, doPlan.timeList[-1])

	# p
	axNum = 4
	#ax[axNum].plot(tPlot[n:], pArr[n:], 'gray', linewidth=1)
	ax[axNum].plot(tPlot[:n], pArr[:n], 'k', linewidth=2)
	ax[axNum].plot(tPlot[n], pArr[n], 'ro', markersize=5)
	ax[axNum].fill_between(tPlot, pArr-errArr, pArr+errArr, color='gray', alpha=0.5)
	ax[axNum].set_ylabel(r'$p$')
	lim1    = np.amin(pArr - errArr)
	lim2    = np.amax(pArr + errArr)
	if lim2 >= 5.0: lim2 = 5.0
	if lim1 == 0.0: lim1 = -lim2*0.06
	ax[axNum].set_ylim(lim1, lim2)
	pMean   = np.round(np.mean(pArr[max(n-50,0):n]),2)
	errMean = np.round(np.mean(errArr[max(n-50,0):n]),2)
	if str(pMean) != 'nan':
		fig.text(0.45, 0.425, r'$p_5=$' + str(pMean) + r'$\pm$' + str(errMean))
	#ax[axNum].axhline(y=1.6, linestyle='--')
	ax[axNum].set_xlim(0, doPlan.timeList[-1])

	# cumulative hist
	axNum = 2
	bins, hist = readerPlan.getCumMassHist(doPlan, n)
	if not isinstance(bins, int): ax[axNum].loglog(bins, hist, 'ko', markersize=3)
	ax[axNum].set_xlabel(r'$M_p$')
	ax[axNum].set_ylabel(r'$N(>M_p)$')
	ax[axNum].set_yscale('log')
	ax[axNum].set_xscale('log')
	ax[axNum].set_ylim(8.e-1, 1.1*max(doPlan.nClumpsList))
	ax[axNum].set_xlim(1.e-4, 1.e1)

	# differential hist
	axNum = 3
	bins, hist = readerPlan.getDiffMassHist(doPlan, n)
	if not isinstance(bins, int) and np.sum(hist>0)>2:
		ax[axNum].loglog(bins, hist, 'ko', markersize=3)
		p, err   = readerPlan.get_p(doPlan, n)
		err2best = 1.e100
		for i in np.arange(-10, 0, 0.01):
			scale = np.power(10.0, i)
			model = scale * np.power(bins, -p)
			logModel = np.log10(model)
			logHist  = np.log10(hist)
			for i in range(len(logHist)):
				if logHist[i]<-1.e-10:
					logHist[i] = 10.0
			err2  = np.square(logModel - logHist)
			err2  = err2 * (hist > 0)
			err2  = np.sum(err2)
			if err2 < err2best:
				err2best  = err2
				bestModel = model
		ax[axNum].loglog(bins, bestModel, 'gray', linestyle='--')
		fig.text(0.865, 0.665, r'$p=$'+str(np.round(p,2)))
	ax[axNum].set_xlabel(r'$M_p$')
	ax[axNum].set_ylabel(r'$dN/dM_p$')
	ax[axNum].set_yscale('log')
	ax[axNum].set_xscale('log')
	ax[axNum].set_ylim(1.e-2, 1.e4)
	ax[axNum].set_xlim(1.e-4, 1.e1)

	# xy and xz scatters
	masses = doPlan.peakArrayList[n][:, 2]
	xs     = doPlan.peakArrayList[n][:, 4]
	ys     = doPlan.peakArrayList[n][:, 5]
	zs     = doPlan.peakArrayList[n][:, 6]
	sizes  = [np.power(3.e3*mass, 1./2.) for mass in masses]
	axNum = 0
	ax[axNum].scatter(xs, ys, s=sizes)
	ax[axNum].set_xlabel(r'$r/h$')
	ax[axNum].set_ylabel(r'y/h')
	ax[axNum].set_ylim(-0.1, 0.1)
	ax[axNum].set_xlim(-0.1, 0.1)
	axNum = 1
	ax[axNum].scatter(xs, zs, s=sizes)
	ax[axNum].set_xlabel(r'$r/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_ylim(-0.025, 0.025)
	ax[axNum].set_xlim(-0.1, 0.1)

	plt.tight_layout()
	plt.savefig(pathSave + "anim_" + str(n) + ".png", bbox_inches='tight')
	plt.close('all')
################################################################################
#if nStart == 0: nStart = doPlan.nFirstClump
for n in range(nStart, doPlan.nt):
	makeAnimFrame(doPlan, n)














#
