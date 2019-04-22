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
#pathBase    = '../../data/kspaceTest/'
#runNameList = ['run150', 'run151', 'run152', 'run153', 'run154']
#colorList   = ['r', 'g', 'm', 'b', 'k']
#pathSave = pathBase + 'plots/pspec_150/'
################################################################################
#pathBase    = '../../data/kspaceTest/'
#runNameList = ['run160', 'run161', 'run162', 'run163', 'run164']
#colorList   = ['r', 'g', 'm', 'b', 'k']
#pathSave = pathBase + 'plots/pspec_160/'
################################################################################
pathBase    = '../../data/kspaceTest/'
runNameList = ['run300', 'run301', 'run302', 'run303']
colorList   = ['r', 'g', 'm', 'b', 'k']
pathSave = pathBase + 'plots/pspec_300/'
################################################################################
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
################################################################################
do3dList = []
for n in range(len(runNameList)):
    path3d     = pathBase + runNameList[n] + '/3d/'
    do3dList.append(reader3d.Data3d(path3d))
################################################################################
vExpo = -1.833
eExpo = 2.0*(vExpo)+2.0
nStart = 2

for n in range(len(do3dList)):
    do3d    = do3dList[n]
    color   = colorList[n]
    psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx', nStart=nStart)
    psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy', nStart=nStart)
    psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz', nStart=nStart)
    psk  = psk_vx  + psk_vy  + psk_vz
    psk*=np.power(freqs, -eExpo)
    psk/=np.mean(psk)
    plt.loglog(freqs, psk, color=color, label=runNameList[n])
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e-3,1.e3)
plt.xlim(freqs[1],freqs[-1])
plt.legend()

index = 2
for prefactor in np.arange(-10,10,0.5):
    plt.plot(freqs, np.power(10,prefactor)*np.power(freqs, index), color=(1,0,0,0.2))
index = 1
for prefactor in np.arange(-10,10,0.5):
    plt.plot(freqs, np.power(10,prefactor)*np.power(freqs, index), color=(0,1,0,0.2))
index = -1
for prefactor in np.arange(-10,10,0.5):
    plt.plot(freqs, np.power(10,prefactor)*np.power(freqs, index), color=(0,0,1,0.2))
index = -2
for prefactor in np.arange(-10,10,0.5):
    plt.plot(freqs, np.power(10,prefactor)*np.power(freqs, index), color=(0,0,0,0.2))


tools.saveAndClear(pathSave + 'adjustedPspecSpheresSumPerts.png', figNum=0)


















#
