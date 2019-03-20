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
pathBase    = '../../data/prodRuns/'
runNameList = ['run320', 'run321', 'run322', 'run323', 'run324']
colorList   = ['r', 'g', 'b', 'k', 'm']
pathSave = pathBase + 'plots/pspec_320/'
################################################################################
#pathBase    = '../../data/prodRuns/'
#runNameList = ['run310', 'run311', 'run312']
#colorList   = ['r', 'g', 'b', 'k', 'm']
#pathSave = pathBase + 'plots/pspec_310/'
################################################################################
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
	psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx')
	psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy')
	psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz')
	psk  = psk_vx  + psk_vy  + psk_vz
	plt.loglog(freqs, psk/psk[1], color, label=runNameList[n])
	################################################################################
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + 'multi_keSpec.png', figNum=0)


















#
