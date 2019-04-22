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
################################################################################
pathBase = str(sys.argv[1])
dt        = 0.1
nStart    = 0
nEnd      = 49
nSamples  = 20
nEachCorr = (nEnd-nStart)-nSamples
nZones    = 16
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/diffusion/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d, dt=dt)
################################################################################
# read in vz
vzList = []
for n in range(nStart, nEnd):
	vz =  do3d.get3d('vz', n)
	vzList.append(vz)
vzMaster = np.asarray(vzList)
################################################################################
# calculate mean dvz^2 for estimation
vzMean2   = np.mean(np.square(vzMaster[:,:,:,32-nZones//2:32+nZones//2]))
################################################################################
# calculate diffusion coeffecient properly with autocorr of velocity
def autoCorr(x):
	result = dt*np.correlate(x, x, mode='full')
	return result[int(result.size/2):]
count = 0
for n in range(0,nSamples):
	for i in range(0,64):
		for j in range(0,64):
			for k in range(32-nZones//2, 32+nZones//2):
				ac = autoCorr(vzMaster[n:n+nEachCorr,i,j,k])
				if n==0 and i==0 and j==0 and k==32-nZones//2:
					acMaster  = ac
					count+=1
				else:
					acMaster += ac
					count+=1
acMaster /= float(count)
time   = dt*np.arange(nStart, nEachCorr)

dtUp   = dt/100.0
timeUp = dtUp*np.arange(nStart, nEachCorr*(dt/dtUp))
acUp   = np.interp(timeUp, time, acMaster)
integral = np.sum(acUp) * dtUp
################################################################################
# plot
plt.figure(0)
plt.plot(timeUp, acUp, 'k');
plt.title(r'$D=$' + str(np.round(integral,6)) + '\n' + r'$vz^2=$' + str(np.round(vzMean2,6)))
tools.saveAndClear(pathSave+'acf.png', figNum=0)
################################################################################








#
