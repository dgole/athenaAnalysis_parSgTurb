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
import athenaReader3d as reader3d
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
################################################################################
# paths
pathBase = str(sys.argv[1])
G        = float(sys.argv[2])
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/planAndGasAnim/'
path3d   = pathBase + '3d/'
do3d     = reader3d.Data3d(path3d)
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan = readerPlan.DataPlan(pathPlan, G=G)
nStart = 50
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
	ax[axNum].set_ylim(-0.05, 1.02)
	ax[axNum].set_ylabel(r'$M_{plan} / M_{par}$')

	# p
	axNum = 4
	#ax[axNum].plot(tPlot[n:], pArr[n:], 'gray', linewidth=1)
	ax[axNum].plot(tPlot[:n], pArr[:n], 'k', linewidth=2)
	ax[axNum].plot(tPlot[n], pArr[n], 'ro', markersize=5)
	ax[axNum].fill_between(tPlot, pArr-errArr, pArr+errArr, color='gray', alpha=0.5)
	ax[axNum].set_ylabel(r'$p$')
	ax[axNum].set_ylim(0.5, 3.0)
	pMean   = np.round(np.mean(pArr[max(n-50,0):n]),2)
	errMean = np.round(np.mean(errArr[max(n-50,0):n]),2)
	if str(pMean) != 'nan':
		fig.text(0.45, 0.39, r'$p_5=$' + str(pMean) + r'$\pm$' + str(errMean))
	#ax[axNum].axhline(y=1.6, linestyle='--')

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
	if not isinstance(bins, int) and np.sum(hist>0)>4:
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
	axNum = 0
	data3d = do3d.get3d('dpar', n)
	data2d = np.mean(data3d, axis=2)
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	aspect  = 0.77
	plotData = np.transpose(np.fliplr(data2d))
	cmapType = 'viridis'
	plotData = np.clip(plotData, 0.01, 2)
	#norm     = colors.LogNorm()
	#ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType), norm=norm)
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	axNum = 1
	data3d = do3d.get3d('dpar', n)
	data2d = np.mean(data3d, axis=1)
	extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
	aspect  = 0.77
	plotData = np.transpose(np.fliplr(data2d))
	plotData = np.clip(plotData, 0.01, 10)
	cmapType = 'viridis'
	#norm     = colors.LogNorm()
	#ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType), norm=norm)
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	masses = doPlan.peakArrayList[n][:, 2]
	xs     = doPlan.peakArrayList[n][:, 4]
	ys     = doPlan.peakArrayList[n][:, 5]
	zs     = doPlan.peakArrayList[n][:, 6]
	sizes  = [np.power(1.e5*mass, 1./2.) for mass in masses]
	axNum = 0
	ax[axNum].scatter(xs, ys, s=sizes, facecolors='none', edgecolors='r')
	ax[axNum].set_xlabel(r'$r/h$')
	ax[axNum].set_ylabel(r'y/h')
	ax[axNum].set_ylim(-0.1, 0.1)
	ax[axNum].set_xlim(-0.1, 0.1)
	axNum = 1
	ax[axNum].scatter(xs, zs, s=sizes, facecolors='none', edgecolors='r')
	ax[axNum].set_xlabel(r'$r/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_ylim(-0.1, 0.1)
	ax[axNum].set_xlim(-0.1, 0.1)



	plt.tight_layout()
	plt.savefig(pathSave + "anim_" + str(n) + ".png", bbox_inches='tight')
	plt.close('all')
################################################################################
for n in range(nStart, doPlan.nt):
	makeAnimFrame(doPlan, n)














#
