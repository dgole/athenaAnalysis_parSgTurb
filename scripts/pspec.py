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
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d     = reader3d.Data3d(path3d)
nStart   = 4
nEnd     = 10
plt.figure(0)
################################################################################
vExpo = -1.833
eExpo = 2.0*(vExpo)+2.0
print("velocity spectrum PL exponent is " + str(vExpo))
print("KE power spectrum PL exponent is " + str(eExpo))
################################################################################
psk_vx, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx', nStart=nStart, nEnd=nEnd)
psk_vy, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy', nStart=nStart, nEnd=nEnd)
psk_vz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz', nStart=nStart, nEnd=nEnd)
psk  = psk_vx  + psk_vy  + psk_vz
psk *= np.power(freqs, -eExpo)
psk /= np.mean(psk)
plt.loglog(freqs, psk)
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e-2,1.e2)
plt.xlim(freqs[1],freqs[-1])
ipeak = np.argmax(psk)
k2 = ipeak + np.argmin(np.absolute(psk[ipeak:]-psk[2]))
plt.axhline(y=psk[2],    color=(0,0,0,0.2), linestyle='--')
plt.axvline(x=freqs[2],  color=(0,0,0,0.2), linestyle='--')
plt.axvline(x=freqs[k2], color=(0,0,0,0.2), linestyle='--')
tools.saveAndClear(pathSave + 'adjustedPspecSpheresSumPerts.png', figNum=0)








#
