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
pathBase    = '../../data/kspaceTest/'
runNameList = ['run30', 'run31', 'run32']
kHighList   = [32, 16, 8]
colorList   = ['ko', 'bo', 'go']
colorList2  = ['k', 'b', 'g']
pathSave    = pathBase + 'plots/'
ms = 2
################################################################################
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
def addFiveThirdsToFig():
	for i in range(-2, 3):
		plt.loglog(freqs, 10**i*np.power(freqs, -5.0/3.0), color='tab:gray',
 				   linestyle='--', linewidth=0.5)
################################################################################
do3dList = []
for n in range(len(runNameList)):
	path3d     = pathBase + runNameList[n] + '/3d/'
	do3dList.append(reader3d.Data3d(path3d))

for n in range(len(do3dList)):
	do3d    = do3dList[n]
	color   = colorList[n]
	kHigh   = kHighList[n] * (1./do.xmax)
	psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoVx')
	psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoVy')
	psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVz')
	psk  = psk_vx  + psk_vy  + psk_vz
	plt.loglog(freqs, psk/psk[1], color, markersize=ms, label='kHighCode='+str(kHigh))
	plt.axvline(x=kHigh, color=colorList2[n])
	#plt.axvline(x=kIn[0]*(6.28/0.2), color=color)
	#plt.axvline(x=kIn[1]*(6.28/0.2), color=color)
	################################################################################
	#psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx')
	#psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy')
	#psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz')
	#psk  = psk_vx  + psk_vy  + psk_vz
	#plt.loglog(freqs, psk/psk[1], 'bo', markersize=ms, label='pert')
	################################################################################
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + 'keSpec.png', figNum=0)


















#
