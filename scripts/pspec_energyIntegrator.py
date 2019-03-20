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
ms   = 2
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
plt.figure(0)
nStart = 0
nEnd   = do3d.nt
################################################################################
def addFiveThirdsToFig():
    for i in range(-2, 3):
        plt.loglog(freqs, 10**i*np.power(freqs, -5.0/3.0), color='tab:gray',
                   linestyle='--', linewidth=0.5)
################################################################################
powerSum = []
keSum = []
ratio = []
ps3d_vz, freqs = reader3d.calcPs(do3d, 'rootRhoVz', 10)
dkVol = np.power(freqs[1]-freqs[0],3)
for n in range(30,80,4):
    ps3d_vz, freqs = reader3d.calcPs(do3d, 'rootRhoDVz', n)
    kez            = np.square(do3d.get3d('rootRhoDvz', n))
    powerSum.append(np.sum(ps3d_vz*dkVol))
    keSum.append   (np.sum(kez))
    ratio.append   (np.sum(ps3d_vz*dkVol)/np.sum(kez))
print(powerSum)
print(keSum)
print(ratio)
################################################################################
#plt.semilogy(range(30,80,4), powerSum)
#plt.semilogy(range(30,80,4), keSum)
#tools.saveAndClear(pathSave + 'keSpec_energyIntegrator.png', figNum=0)
#plt.semilogy(range(30,80,4), ratio)
#tools.saveAndClear(pathSave + 'keSpec_ratio.png', figNum=0)


















#plt.xlabel(r'$|\mathbf{k}|$')
#plt.ylabel('Power')
#tools.saveAndClear(pathSave + 'keSpec_energyIntegrator.png', figNum=0)
################################################################################
#print(np.sum(kez))
#print(np.sum(psk_vz*dk))
#print(np.sum(np.sqrt(psk_vz)*dk))
#print(np.sum(np.sqrt(psk_vz*dk)))
#print(np.sum(psk_vz))
















#
