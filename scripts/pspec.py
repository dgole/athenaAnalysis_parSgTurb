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
ms = 2
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
def addFiveThirdsToFig():
    for i in range(-2, 3):
        plt.loglog(freqs, 10**i*np.power(freqs, -5.0/3.0), color='tab:gray',
                   linestyle='--', linewidth=0.5)
################################################################################
psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoVx')
psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoVy')
psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVz')
psk  = psk_vx  + psk_vy  + psk_vz
plt.loglog(freqs, psk/psk[1], 'ko', markersize=ms, label='normal')
################################################################################
psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx')
psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy')
psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz')
psk  = psk_vx  + psk_vy  + psk_vz
plt.loglog(freqs, psk/psk[1], 'bo', markersize=ms, label='pert')
################################################################################
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + 'keSpec.png', figNum=0)


















#
