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
for n in range(nStart, nEnd, 2):
    psk_vx, freqs = reader3d.psProfiles(do3d, 'rootRhoVx', n)
    psk_vy, freqs = reader3d.psProfiles(do3d, 'rootRhoVy', n)
    psk_vz, freqs = reader3d.psProfiles(do3d, 'rootRhoVz', n)
    psk  = psk_vx  + psk_vy  + psk_vz
    plt.loglog(freqs, psk/psk[1], color=tools.getColor(n, 0, nEnd), label='n='+str(n))
################################################################################
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
addFiveThirdsToFig()
#plt.legend(loc=(1.01,0.0))
plt.ylim(1.e-6, 1.e1)
tools.saveAndClear(pathSave + 'keSpec_timeEvo.png', figNum=0)


















#
