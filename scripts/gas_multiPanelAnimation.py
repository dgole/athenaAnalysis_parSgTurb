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
pathSave = pathBase + 'plots/gasAnim/'
path3d   = pathBase + '3d/'
do3d     = reader3d.Data3d(path3d)
pathParHst = pathBase + '/'
doPhst = readerPhst.DataPhst(pathParHst)
path1d   = pathBase + '1d/'
do1d     = reader1d.Data1d(path1d)
do1d.addCol(reader1d.dv, 'dv', r'$\delta v$')
do1d.addCol(reader1d.dvz, 'dvz', r'$\delta v_z$')
do1d.addCol(reader1d.alpha,  'alpha',  r'$\alpha$'  )
do1d.addCol(reader1d.alphaz, 'alphaz', r'$\alpha_z$')
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
nStart   = 0
alpha_in = float(sys.argv[2])
ts       = 0.3
aspect   = 0.8
################################################################################
def applyLims(data, lim1, lim2):
	np.clip(data, lim1, lim2)
	data[0,0]=lim1
	data[0,1]=lim2
	return data
################################################################################
def makeAnimFrame(n):
	print('saving anim frame for n = ' + str(n))
	fig = plt.figure(figsize=(8.03,13.0), dpi=100)
	ax = []
	ax.append(plt.subplot2grid((8, 2), (0, 0), rowspan=2))
	ax.append(plt.subplot2grid((8, 2), (0, 1), rowspan=2))
	ax.append(plt.subplot2grid((8, 2), (2, 0), rowspan=2))
	ax.append(plt.subplot2grid((8, 2), (2, 1), rowspan=2))
	ax.append(plt.subplot2grid((8, 2), (4, 0), rowspan=2))
	ax.append(plt.subplot2grid((8, 2), (4, 1), rowspan=2))
	ax.append(plt.subplot2grid((8, 2), (6, 0), colspan=2))
	ax.append(plt.subplot2grid((8, 2), (7, 0), colspan=2))

	tPlot = do1d.t

	# Alphas
	axNum = 7
	alphaAvg = np.mean(do1d.data['alpha'], axis=1)
	alphazAvg = np.mean(do1d.data['alphaz'], axis=1)
	ax[axNum].semilogy(tPlot[n:], alphaAvg[n:], 'gray', linewidth=1)
	ax[axNum].semilogy(tPlot[:n], alphaAvg[:n], 'k', linewidth=2)
	ax[axNum].semilogy(tPlot[n],  alphaAvg[n], 'ro', markersize=5, label=r'$\alpha$')
	ax[axNum].semilogy(tPlot[n:], alphazAvg[n:], 'gray', linewidth=1)
	ax[axNum].semilogy(tPlot[:n], alphazAvg[:n], 'k', linewidth=2)
	ax[axNum].semilogy(tPlot[n],  alphazAvg[n], 'bo', markersize=5, label=r'$\alpha_z$')
	ax[axNum].axhline(y=alpha_in, linestyle=':', color='k', label=r'$\alpha_{in}$')
	ax[axNum].set_ylim(3.e-2*alpha_in, 2.0*alpha_in)
	ax[axNum].set_ylabel(r'$\alpha$')
	ax[axNum].set_xlabel(r'$t\Omega$')
	#ax[axNum].legend()

	# particle scale height
	axNum = 6
	parh    = doPhst.data['zvar']
	parhEst = np.sqrt(np.mean(alphazAvg[do1d.nt//2:])/ts)
	ax[axNum].semilogy(tPlot[n:], parh[n:], 'gray', linewidth=1)
	ax[axNum].semilogy(tPlot[:n], parh[:n], 'k', linewidth=2)
	ax[axNum].semilogy(tPlot[n],  parh[n], 'ro', markersize=5)
	ax[axNum].axhline(y=parhEst, linestyle=':', color='k')
	ax[axNum].set_ylim(0.2/100, 0.2)
	ax[axNum].set_ylabel(r'$h_p$/h')
	ax[axNum].set_xlabel(r'$t\Omega$')

	# dpar
	axNum = 0
	data3d = do3d.get3d('dpar', do3d.gettindex(tPlot[n]))
	data2d = np.mean(data3d, axis=2)
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	plotData = np.transpose(np.fliplr(data2d))
	plotData = applyLims(plotData, 0.0, 1.5)
	cmapType = 'viridis'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$y/h$')
	ax[axNum].set_title('dpar (column)')

	axNum = 1
	data3d = do3d.get3d('dpar', do3d.gettindex(tPlot[n]))
	data2d = np.mean(data3d, axis=1)
	extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
	plotData = np.transpose(np.fliplr(data2d))
	plotData = applyLims(plotData, 0.0, 1.5)
	cmapType = 'viridis'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_title('dpar (column)')

	# rho (gas)
	axNum = 2
	data3d = do3d.get3d('rho', do3d.gettindex(tPlot[n]))
	data2d = data3d[:,:,do3d.nz//2]
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	plotData = np.transpose(np.fliplr(data2d))
	plotData = applyLims(plotData, 0.98, 1.0)
	cmapType = 'viridis'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$y/h$')
	ax[axNum].set_title(r'$\rho$' + ' (slice)')

	axNum = 3
	data3d = do3d.get3d('rho', do3d.gettindex(tPlot[n]))
	data2d = data3d[:,do3d.ny//2,:]
	extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
	plotData = np.transpose(np.fliplr(data2d))
	plotData = applyLims(plotData, 0.98, 1.0)
	cmapType = 'viridis'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_title(r'$\rho$' + ' (slice)')

	# vz (gas)
	axNum = 4
	data3d = do3d.get3d('vz', do3d.gettindex(tPlot[n]))
	data2d = data3d[:,:,do3d.nz//2]
	extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
	plotData = np.transpose(np.fliplr(data2d))
	plotData = applyLims(plotData, -0.03, 0.03)
	cmapType = 'coolwarm'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$y/h$')
	ax[axNum].set_title(r'$v_z$' + ' (slice)')

	axNum = 5
	data3d = do3d.get3d('vz', do3d.gettindex(tPlot[n]))
	data2d = data3d[:,do3d.ny//2,:]
	extent = [-do3d.xmax, do3d.xmax, -do3d.zmax, do3d.zmax]
	plotData = np.transpose(np.fliplr(data2d))
	plotData = applyLims(plotData, -0.03, 0.03)
	cmapType = 'coolwarm'
	ax[axNum].imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
	ax[axNum].set_xlabel(r'$x/h$')
	ax[axNum].set_ylabel(r'$z/h$')
	ax[axNum].set_title(r'$v_z$' + ' (slice)')


	plt.tight_layout()
	plt.savefig(pathSave + "anim_" + str(n) + ".png", bbox_inches='tight')
	plt.close('all')
################################################################################
for n in range(nStart, do1d.nt, 10):
	makeAnimFrame(n)














#
