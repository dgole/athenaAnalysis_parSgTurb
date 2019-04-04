#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import math
import sys
sys.path.append('../python')
import athenaReader3d as reader3d
import athenaTools as tools
import athenaReaderPhst as readerPhst
import planOutputReader as readerPlan
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp
################################################################################
myNpc        = int(sys.argv[1])
G            = float(sys.argv[2])
pathBase     = str(sys.argv[3])
dt3d         = 1.0
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/masterAnim/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d, dt=dt3d)
pathPlan = pathBase + 'planOutput/'
if do3d.nt > 15: nAvgEnd = 10*int(1.0/dt3d)
else:            nAvgEnd = 2
################################################################################
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
planPlots = True
try:
	doPlan = readerPlan.DataPlan(pathPlan, G=G)
except:
	planPlots = False
if np.sum(do3d.get3d('dpar', 0)) > 1.e-20:
	dparPlots = True
else:
	dparPlots = False
################################################################################
if planPlots:
	pList   = []
	errList = []
	for n in range(0, doPlan.nt):
		p, err = readerPlan.get_p(doPlan, n)
		pList.append(p)
		errList.append(err)
	pArr   = np.asarray(pList)
	errArr = np.asarray(errList)
################################################################################
def prepPlotData(plotData, lim1, lim2):
	plotData = np.transpose(np.fliplr(plotData))
	plotData = np.clip(plotData, lim1, lim2)
	plotData[0,0] = lim1
	plotData[0,1] = lim2
	return plotData
################################################################################
plotDataDict = {}
avgDataDict  = {}
for key in ['dvx', 'dvy', 'dvz', 'dv', 'drho']:
	plotDataDict[key] = np.zeros(do3d.nt)
	avgDataDict[key]  = np.zeros(do3d.nz)
	nAvg = 0
	for n in range(do3d.nt):
		data = do3d.get3d(key, n)
		plotDataDict[key][n] = np.mean(np.absolute(data))
		if n >= do3d.nt-nAvgEnd:
			avgDataDict[key]+= np.mean(np.absolute(data), axis=(0,1))
			nAvg+=1
	avgDataDict[key]/=nAvg
for key in ['vx', 'vy', 'vz']:
	avgDataDict[key]  = np.zeros(do3d.nz)
	nAvg = 0
	for n in range(do3d.nt-nAvgEnd, do3d.nt):
		data = do3d.get3d(key, n)
		avgDataDict[key]+= np.mean(data, axis=(0,1))
		nAvg+=1
	avgDataDict[key]/=nAvg

key = 'pspec'
psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx', nStart=do3d.nt-nAvgEnd)
psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy', nStart=do3d.nt-nAvgEnd)
psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz', nStart=do3d.nt-nAvgEnd)
psk  = psk_vx  + psk_vy  + psk_vz
avgDataDict[key] = psk

################################################################################

def makeAnimFrame(self, n):
	print('saving anim frame for n = ' + str(n))
	sizeFactor = 2.5
	fig = plt.figure(figsize=(8.545*sizeFactor, 5*sizeFactor), dpi=120)
	ax = []

	ax.append(plt.subplot2grid((5, 7), (0, 0), rowspan=1)) # row 1 first 3
	ax.append(plt.subplot2grid((5, 7), (0, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (0, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (1, 0), rowspan=1)) # row 2 first 3
	ax.append(plt.subplot2grid((5, 7), (1, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (1, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (2, 0), rowspan=1)) # row 3 first 3
	ax.append(plt.subplot2grid((5, 7), (2, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (2, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (3, 0), rowspan=1)) # row 4 first 3
	ax.append(plt.subplot2grid((5, 7), (3, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (3, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (4, 0), colspan=3)) # bottom time evo

	ax.append(plt.subplot2grid((5, 7), (0, 3), rowspan=1)) # col 4 going down
	ax.append(plt.subplot2grid((5, 7), (1, 3), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (2, 3), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (3, 3), rowspan=1))
	ax.append(plt.subplot2grid((5, 7), (4, 3), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (2, 4), colspan=3)) # 3 time evos on right
	ax.append(plt.subplot2grid((5, 7), (3, 4), colspan=3))
	ax.append(plt.subplot2grid((5, 7), (4, 4), colspan=3))

	ax.append(plt.subplot2grid((5, 7), (0, 4), rowspan=1)) # above right time evos
	ax.append(plt.subplot2grid((5, 7), (1, 4), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (0, 5), rowspan=1)) # above right time evos
	ax.append(plt.subplot2grid((5, 7), (1, 5), rowspan=1))

	ax.append(plt.subplot2grid((5, 7), (0, 6), rowspan=1)) # above right time evos
	ax.append(plt.subplot2grid((5, 7), (1, 6), rowspan=1))

	axNumDict = {'vx':0, 'vy':1, 'vz':2}
	for key in ['vx', 'vy', 'vz']:
		axNum = axNumDict[key]
		ax[axNum].set_xlabel(r'$x/h$')
		ax[axNum].set_ylabel(r'$y/h$')
		ax[axNum].set_title(do3d.header[key])
		data3d = do3d.get3d(key, n)
		data2d = data3d[:,:,do3d.nz//2]
		extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
		aspect  = 0.84
		plotData = prepPlotData(data2d, np.amin(data2d), np.amax(data2d))
		cmapType = 'coolwarm'
		ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	axNumDict = {'vx':3, 'vy':4, 'vz':5}
	for key in ['vx', 'vy', 'vz']:
		axNum = axNumDict[key]
		ax[axNum].set_xlabel(r'$x/h$')
		ax[axNum].set_ylabel(r'$z/h$')
		ax[axNum].set_title(do3d.header[key])
		data3d = do3d.get3d(key, n)
		data2d = data3d[:,do3d.ny//2,:]
		extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
		plotData = prepPlotData(data2d, np.amin(data2d), np.amax(data2d))
		cmapType = 'coolwarm'
		ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	axNumDict = {'dvx':6, 'dvy':7, 'dvz':8}
	for key in ['dvx', 'dvy', 'dvz']:
		axNum = axNumDict[key]
		ax[axNum].set_xlabel(r'$z/h$')
		ax[axNum].set_ylabel(do3d.header[key])
		plotData = np.mean(np.absolute(do3d.get3d(key, n)), axis=(0,1))
		ax[axNum].semilogy(do3d.z, plotData, 'k', linewidth=2)
		ax[axNum].semilogy(do3d.z, avgDataDict[key], 'gray', linewidth=1)
		ax[axNum].set_ylim(1.e-4, 3.e-1)
		#ax[axNum].axhline(y=np.sqrt(alphaDesired/3.0), linestyle='--', color='gray')
		if key == 'dvz': parhTheo = np.mean(plotData)/np.sqrt(0.3)

	axNumDict = {'vx':9, 'vy':10, 'vz':11}
	for key in ['vx', 'vy', 'vz']:
		axNum = axNumDict[key]
		ax[axNum].set_xlabel(r'$z/h$')
		ax[axNum].set_ylabel(do3d.header[key])
		plotData = np.mean(do3d.get3d(key, n), axis=(0,1))
		ax[axNum].plot(do3d.z, plotData, 'k', linewidth=2)
		ax[axNum].plot(do3d.z, avgDataDict[key], 'gray', linewidth=1)
		ax[axNum].set_ylim(-0.2,0.2)

	axNum = 12
	ax[axNum].set_xlabel(r'$t \Omega$')
	ax[axNum].set_ylabel(r'$\delta v$')
	colors = {'dvx':'r', 'dvy':'g', 'dvz': 'b', 'dv': 'k'}
	for key in ['dvx', 'dvy', 'dvz', 'dv']:
		plotData = plotDataDict[key]
		ax[axNum].semilogy(do3d.t[n:], plotData[n:], 'gray', linewidth=1)
		ax[axNum].semilogy(do3d.t[:n+1], plotData[:n+1], colors[key], linewidth=2)
		ax[axNum].semilogy(do3d.t[n], plotData[n], colors[key]+'o', markersize=5, label=do3d.header[key])
	#ax[axNum].axhline(y=np.sqrt(alphaDesired), linestyle='--', color='gray')
	ax[axNum].set_ylim(8.e-5, 3.e-1)
	ax[axNum].set_xlim(0.0, np.amax(do3d.t))
	ax[axNum].legend(loc=(0.5,0.05), ncol=4)

	axNum = 13
	key   = 'drho'
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$y/h$')
	ax[axNum].set_title(do3d.header[key])
	data3d = do3d.get3d(key, n)
	data2d = data3d[:,:,do3d.nz//2]
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	plotData = prepPlotData(data2d, np.amin(data2d), np.amax(data2d))
	cmapType = 'PuOr'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	axNum = 14
	key   = 'drho'
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_title(do3d.header[key])
	data3d = do3d.get3d(key, n)
	data2d = data3d[:,do3d.ny//2,:]
	extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
	plotData = prepPlotData(data2d, np.amin(data2d), np.amax(data2d))
	cmapType = 'PuOr'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	axNum  = 15
	key    = 'drho'
	ax[axNum].set_xlabel(r'$z/h$')
	ax[axNum].set_ylabel(do3d.header[key])
	plotData = np.mean(np.absolute(do3d.get3d(key, n)), axis=(0,1))
	ax[axNum].semilogy(do3d.z, plotData, 'k', linewidth=2)
	ax[axNum].semilogy(do3d.z, avgDataDict[key], 'gray', linewidth=1)
	ax[axNum].set_ylim(1.e-5, 1.e-2)

	nPspec = int((n//2)*2.0)
	axNum  = 16
	key    = 'drho'
	acf = reader3d.acfMean(do3d, key, nStart=nPspec, nEnd=nPspec+1)
	plotData = prepPlotData(acf, -1.0, 1.0)
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	ax[axNum].set_xlabel(r'$\delta x/h$')
	ax[axNum].set_ylabel(r'$\delta y/h$')
	ax[axNum].set_title(do3d.header[key] + ' ACF')
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap('PuOr'))

	axNum  = 17
	key    = 'pspec'
	psk_vx, freqs = reader3d.psProfile(do3d, 'rootRhoDvx', n=nPspec)
	psk_vy, freqs = reader3d.psProfile(do3d, 'rootRhoDvy', n=nPspec)
	psk_vz, freqs = reader3d.psProfile(do3d, 'rootRhoDvz', n=nPspec)
	psk  = psk_vx  + psk_vy  + psk_vz
	ax[axNum].set_xlabel(r'$k$')
	ax[axNum].set_ylabel('Power')
	ax[axNum].loglog(freqs[1:], psk[1:], 'k', linewidth=2)
	ax[axNum].loglog(freqs[1:], avgDataDict[key][1:], 'gray', linewidth=1)
	for i in range(-100,100):
		factor = np.power(10.0,float(i)*0.5)
		ax[axNum].loglog(freqs[1:], factor*np.power(freqs[1:], -2./3.),
						 color=(1,0,0,0.3), linewidth=1, linestyle='-')
		ax[axNum].loglog(freqs[1:], factor*np.power(freqs[1:], -5./3.),
						 color=(0,1,0,0.3), linewidth=1, linestyle='-')
		ax[axNum].loglog(freqs[1:], factor*np.power(freqs[1:], -8./3.),
						 color=(0,0,1,0.3), linewidth=1, linestyle='-')
	ax[axNum].set_ylim(0.3*np.amin(avgDataDict[key][1:]), 3.0*np.amax(avgDataDict[key][1:]))
	ax[axNum].set_xlim(freqs[1], freqs[-1])

	#nPlan = int(n//10)

	# p
	if planPlots:
		axNum = 18
		#ax[axNum].plot(tPlot[n:], pArr[n:], 'gray', linewidth=1)
		ax[axNum].plot(doPlan.time[:n+1], pArr[:n+1], 'k', linewidth=2)
		ax[axNum].plot(doPlan.time[n], pArr[n], 'ro', markersize=5)
		ax[axNum].fill_between(doPlan.time, pArr-errArr, pArr+errArr, color='gray', alpha=0.5)
		ax[axNum].set_ylabel(r'$p$')
		lim1    = np.amin(pArr - errArr)
		lim2    = np.amax(pArr + errArr)
		if lim2 >= 4.0: lim2 = 4.0
		if lim1 == 0.0: lim1 = -lim2*0.06
		ax[axNum].set_ylim(lim1, lim2)
		pMean   = np.round(np.mean(pArr[max(n-50,0):n]),2)
		errMean = np.round(np.mean(errArr[max(n-50,0):n]),2)
		if str(pMean) != 'nan':
			fig.text(0.78, 0.56, r'$p_5=$' + str(pMean) + r'$\pm$' + str(errMean))
		ax[axNum].axhline(y=1.6, linestyle='--', color='k')
		ax[axNum].set_xlim(0, doPlan.timeList[-1])

		# mass frac
		axNum = 19
		mFracList = []
		for n1 in range(0, len(doPlan.peakArrayList)):
			mNow = np.sum(doPlan.peakArrayList[n1][:,2])
			mFracNow = mNow / doPlan.mParTot
			mFracList.append(mFracNow)
		ax[axNum].plot(doPlan.timeList[n:], mFracList[n:], 'gray', linewidth=1)
		ax[axNum].plot(doPlan.timeList[:n+1], mFracList[:n+1], 'k', linewidth=2)
		ax[axNum].plot(doPlan.timeList[n], mFracList[n], 'ro', markersize=5)
		#ax[axNum].get_xaxis().set_visible(False)
		ax[axNum].set_ylim(-0.05, 1.02)
		ax[axNum].set_ylabel(r'$M_{plan} / M_{par}$')
		ax[axNum].set_xlim(0, doPlan.timeList[-1])

		# Nclumps
		axNum = 20
		ax[axNum].plot(doPlan.time[n:], doPlan.nClumpsList[n:], 'gray', linewidth=1)
		ax[axNum].plot(doPlan.time[:n+1], doPlan.nClumpsList[:n+1], 'k', linewidth=2)
		ax[axNum].plot(doPlan.time[n], doPlan.nClumpsList[n], 'ro', markersize=5)
		ax[axNum].set_ylabel(r'$N_{clumps}$')
		ax[axNum].set_xlabel(r'$t \Omega$')
		ax[axNum].set_xlim(0, doPlan.timeList[-1])

		# xy and xz scatters
		masses = doPlan.peakArrayList[n][:, 2]
		xs     = doPlan.peakArrayList[n][:, 4]
		ys     = doPlan.peakArrayList[n][:, 5]
		zs     = doPlan.peakArrayList[n][:, 6]
		sizes  = [np.power(1.e4*mass, 1./2.) for mass in masses]
		axNum = 23
		ax[axNum].scatter(xs, ys, s=sizes)
		ax[axNum].set_xlabel(r'$r/h$')
		ax[axNum].set_ylabel(r'y/h')
		ax[axNum].set_ylim(-do3d.ymax, do3d.ymax)
		ax[axNum].set_xlim(-do3d.xmax, do3d.xmax)
		ax[axNum].set_title('Clumps from PLAN')
		axNum = 24
		ax[axNum].scatter(xs, zs, s=sizes)
		ax[axNum].set_xlabel(r'$r/h$')
		ax[axNum].set_ylabel(r'$z/h$')
		ax[axNum].set_ylim(-do3d.zmax, do3d.zmax)
		ax[axNum].set_xlim(-do3d.xmax, do3d.xmax)
		ax[axNum].set_title('Clumps from PLAN')

		# cumulative hist
		axNum = 25
		bins, hist = readerPlan.getCumMassHist(doPlan, n)
		if not isinstance(bins, int): ax[axNum].loglog(bins, hist, 'ko', markersize=3)
		ax[axNum].set_xlabel(r'$M_p$')
		ax[axNum].set_ylabel(r'$N(>M_p)$')
		ax[axNum].set_yscale('log')
		ax[axNum].set_xscale('log')
		ax[axNum].set_ylim(8.e-1, 1.1*max(doPlan.nClumpsList))
		ax[axNum].set_xlim(1.e-4, 1.e1)

		# differential hist
		axNum = 26
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
			#ax[axNum].loglog(bins, bestModel, 'gray', linestyle='--')
			fig.text(0.9, 0.66, r'$p=$'+str(np.round(p,2)))
		ax[axNum].set_xlabel(r'$M_p$')
		ax[axNum].set_ylabel(r'$dN/dM_p$')
		ax[axNum].set_yscale('log')
		ax[axNum].set_xscale('log')
		ax[axNum].set_ylim(1.e-2, 1.e4)
		ax[axNum].set_xlim(1.e-4, 1.e1)

	if dparPlots == True:
		axNum = 21
		key   = 'dpar'
		ax[axNum].set_xlabel(r'$x/h$')
		ax[axNum].set_ylabel(r'$y/h$')
		ax[axNum].set_title(do3d.header[key])
		data3d = do3d.get3d(key, n)
		data2d = np.mean(data3d, axis=2)
		extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
		plotData = prepPlotData(data2d, np.percentile(data2d,95)/10.0, np.percentile(data2d,95)*2.0)
		cmapType = 'viridis'
		ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

		axNum = 22
		key   = 'dpar'
		ax[axNum].set_xlabel(r'$x/h$')
		ax[axNum].set_ylabel(r'$z/h$')
		ax[axNum].set_title(do3d.header[key])
		data3d = do3d.get3d(key, n)
		data2d = np.mean(data3d, axis=1)
		extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
		plotData = prepPlotData(data2d, np.percentile(data2d,95)/10.0, np.percentile(data2d,95)*2.0)
		cmapType = 'viridis'
		ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
		data1d  = np.mean(data3d, axis=(0,1))
		sum=0; i=1;
		while sum < 0.68*np.sum(data1d):
			sum += data1d[do3d.nz//2 - i] + data1d[do3d.nz//2 + (i-1)]
			parhAct = i*do3d.dx
			i+=1
		ax[axNum].axhline(y=-parhTheo, color=(1,0,0,0.5), linestyle='--');
		ax[axNum].axhline(y=parhTheo, color=(1,0,0,0.5), linestyle='--');
		ax[axNum].axhline(y=-parhAct, color=(1,0,0,0.5), linestyle='-');
		ax[axNum].axhline(y=parhAct, color=(1,0,0,0.5), linestyle='-');

	plt.tight_layout()
	plt.savefig(pathSave + "anim_" + str(n) + ".png", bbox_inches='tight')
	plt.close('all')

################################################################################
jobList = []
for n in range(0, do3d.nt, 1):
	job = mp.Process(target=makeAnimFrame, args=(do3d, n))
	jobList.append(job)

while len(jobList)>0:
	theseJobs = []
	for n in range(myNpc):
		try: theseJobs.append(jobList.pop(0))
		except: a=1
	for job in theseJobs:
		job.start()
	for job in theseJobs:
		job.join()



#
