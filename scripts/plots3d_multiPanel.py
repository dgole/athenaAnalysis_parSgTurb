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
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing as mp
################################################################################
myNpc       = int(sys.argv[1])
#parInitProc = str(sys.argv[2])
parInitProc = 'False'
alphaDesired = float(sys.argv[2])
pathBase  	= str(sys.argv[3])
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/plots3dAnim/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
################################################################################
#plotDataDict = {}

if parInitProc=='True':
	manager     = mp.Manager()
	dataDict    = manager.dict()
	for key in ['dvx', 'dvy', 'dvz', 'dv']:
		dataDict[key] = manager.list(range(do3d.nt))
	def processOne(do3d, key, n, result):
		result[n] = np.mean(np.absolute(do3d.get3d(key, n)))
	jobList = []
	for key in ['dvx', 'dvy', 'dvz', 'dv']:
		for n in range(do3d.nt):
			job = mp.Process(target=processOne, args=(do3d, key, n, dataDict[key]))
			jobList.append(job)
	while len(jobList)>0:
		print(len(jobList))
		theseJobs = []
		for n in range(myNpc):
			try: theseJobs.append(jobList.pop(0))
			except: a=1
		for job in theseJobs:
			job.start()
		for job in theseJobs:
			job.join()
	plotDataDict = {}
	for key in ['dvx', 'dvy', 'dvz', 'dv']:
		plotDataDict[key] = np.asarray(dataDict[key])
else:
	plotDataDict = {}
	avgDataDict  = {}
	for key in ['dvx', 'dvy', 'dvz', 'dv']:
		plotDataDict[key] = np.zeros(do3d.nt)
		avgDataDict[key]  = np.zeros(do3d.nz)
		nAvg = 0
		for n in range(do3d.nt):
			data = do3d.get3d(key, n)
			plotDataDict[key][n] = np.mean(np.absolute(data))
			if n >= do3d.nt-10:
				avgDataDict[key]+= np.mean(np.absolute(data), axis=(0,1))
				nAvg+=1
		avgDataDict[key]/=nAvg
	for key in ['vx', 'vy', 'vz']:
		avgDataDict[key]  = np.zeros(do3d.nz)
		nAvg = 0
		for n in range(do3d.nt-10, do3d.nt):
			data = do3d.get3d(key, n)
			avgDataDict[key]+= np.mean(data, axis=(0,1))
			nAvg+=1
		avgDataDict[key]/=nAvg



################################################################################

def makeAnimFrame(self, n):
	print('saving anim frame for n = ' + str(n))
	sizeFactor = 2.5
	fig = plt.figure(figsize=(5*sizeFactor, 5*sizeFactor), dpi=80)
	ax = []

	ax.append(plt.subplot2grid((5, 4), (0, 0), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (0, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (0, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 4), (1, 0), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (1, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (1, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 4), (2, 0), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (2, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (2, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 4), (3, 0), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (3, 1), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (3, 2), rowspan=1))

	ax.append(plt.subplot2grid((5, 4), (4, 0), colspan=3))

	ax.append(plt.subplot2grid((5, 4), (0, 3), rowspan=1))
	ax.append(plt.subplot2grid((5, 4), (1, 3), rowspan=1))
	#ax.append(plt.subplot2grid((4, 4), (2, 3), rowspan=1))
	#ax.append(plt.subplot2grid((4, 4), (3, 3), rowspan=1))


	#
	axNum = 12
	ax[axNum].set_xlabel(r'$t \Omega$')
	ax[axNum].set_ylabel(r'$\delta v$')
	colors = {'dvx':'r', 'dvy':'g', 'dvz': 'b', 'dv': 'k'}
	for key in ['dvx', 'dvy', 'dvz', 'dv']:
		plotData = plotDataDict[key]
		ax[axNum].semilogy(do3d.t[n:], plotData[n:], 'gray', linewidth=1)
		ax[axNum].semilogy(do3d.t[:n+1], plotData[:n+1], colors[key], linewidth=2)
		ax[axNum].semilogy(do3d.t[n], plotData[n], colors[key]+'o', markersize=5, label=do3d.header[key])
	#limBase = np.mean(plotData[do3d.nt//2:])
	limBase  = 1.e-2
	ax[axNum].axhline(y=np.sqrt(alphaDesired), linestyle='--', color='gray')
	ax[axNum].set_ylim(limBase/30.0, limBase*5.0)
	ax[axNum].legend(loc=(0.90,0.05))

	axNumDict = {'vx':0, 'vy':1, 'vz':2}
	for key in ['vx', 'vy', 'vz']:
		axNum = axNumDict[key]
		ax[axNum].set_xlabel(r'$x/h$')
		ax[axNum].set_ylabel(r'$y/h$')
		ax[axNum].set_title(do3d.header[key])
		data3d = do3d.get3d(key, n)
		data2d = data3d[:,:,do3d.nz//2]
		extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
		aspect  = 0.85
		plotData = np.transpose(np.fliplr(data2d))
		plotData = np.clip(plotData, -limBase*5.0, limBase*5.0)
		plotData[0,0] = -limBase*5.0
		plotData[0,1] =  limBase*5.0
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
		aspect  = 0.85
		plotData = np.transpose(np.fliplr(data2d))
		plotData = np.clip(plotData, -limBase*5.0, limBase*5.0)
		plotData[0,0] = -limBase*5.0
		plotData[0,1] =  limBase*5.0
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
		ax[axNum].set_ylim(limBase/30.0, 3.0*limBase)
		ax[axNum].axhline(y=np.sqrt(alphaDesired/3.0), linestyle='--', color='gray')

	axNumDict = {'vx':9, 'vy':10, 'vz':11}
	for key in ['vx', 'vy', 'vz']:
		axNum = axNumDict[key]
		ax[axNum].set_xlabel(r'$z/h$')
		ax[axNum].set_ylabel(do3d.header[key])
		plotData = np.mean(do3d.get3d(key, n), axis=(0,1))
		ax[axNum].plot(do3d.z, plotData, 'k', linewidth=2)
		ax[axNum].plot(do3d.z, avgDataDict[key], 'gray', linewidth=1)
		ax[axNum].set_ylim(-limBase*2.0, limBase*5.0)

	axNum = 13
	key   = 'dpar'
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$y/h$')
	ax[axNum].set_title(do3d.header[key])
	data3d = do3d.get3d(key, n)
	data2d = np.mean(data3d, axis=2)
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	aspect  = 0.85
	plotData = np.transpose(np.fliplr(data2d))
	plotData = np.clip(plotData, 1.e-1, 3.0)
	plotData[0,0] = 1.e-1
	plotData[0,1] = 3.0
	cmapType = 'viridis'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))

	axNum = 14
	key   = 'dpar'
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_title(do3d.header[key])
	data3d = do3d.get3d(key, n)
	data2d = np.mean(data3d, axis=1)
	extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
	aspect  = 0.85
	plotData = np.transpose(np.fliplr(data2d))
	plotData = np.clip(plotData, 1.e-1, 6.0)
	plotData[0,0] = 1.e-1
	plotData[0,1] = 6.0
	cmapType = 'viridis'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))




	plt.tight_layout()
	plt.savefig(pathSave + "anim_" + str(n) + ".png", bbox_inches='tight')
	plt.close('all')

################################################################################
jobList = []
for n in range(0, do3d.nt):
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
