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
################################################################################
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
vExpo = -1.33
eExpo = 2.0*(vExpo)+1.0
def addFiveThirdsToFig():
    for i in np.arange(-40, 40, 0.25):
        plt.loglog(freqs, 10**i*np.power(freqs, eExpo), color='tab:gray',
                   linestyle='--', linewidth=0.5)
################################################################################
nStart = 2
nEnd   = do3d.nt

for key in ['rootRhoDvx','rootRhoDvy','rootRhoDvz','rootRhoVx','rootRhoVy','rootRhoVz']:
    count=0
    for n in range(nStart, nEnd, 2):
        if n == nStart:
            ps, freqs = reader3d.psProfile(do3d, key, n=n)
            count+=1
        else:
            ps1, freqs1 = reader3d.psProfile(do3d, key, n=n)
            ps+=ps1
            count+=1
    ps/=count
    plt.loglog(freqs, ps)
    plt.xlabel(r'$|\mathbf{k}|$')
    plt.ylabel('Power')
    plt.ylim(1.e-12,1.e-6)
    addFiveThirdsToFig()
    plt.legend()
    tools.saveAndClear(pathSave + 'pspecSpheresSum_' + key + '.png', figNum=0)

n1=0
colorList = ['r','g','b']
for key in ['rootRhoDvx','rootRhoDvy','rootRhoDvz']:
    color = colorList[n1]
    count=0
    for n in range(nStart, nEnd, 2):
        if n == nStart:
            ps, freqs = reader3d.psProfile(do3d, key, n=n)
            count+=1
        else:
            ps1, freqs1 = reader3d.psProfile(do3d, key, n=n)
            ps+=ps1
            count+=1
    ps/=count
    ps*=np.power(freqs, -eExpo)
    ps/=np.mean(ps)
    plt.loglog(freqs, ps, label=key, color=color)
    n1+=1
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e-2,1.e2)
plt.legend()
tools.saveAndClear(pathSave + 'adjustedPspecSpheresSumPerts.png', figNum=0)

n1=0
colorList = ['r','g','b']
for key in ['rootRhoVx','rootRhoVy','rootRhoVz']:
    color = colorList[n1]
    count=0
    for n in range(nStart, nEnd, 2):
        if n == nStart:
            ps, freqs = reader3d.psProfile(do3d, key, n=n)
            count+=1
        else:
            ps1, freqs1 = reader3d.psProfile(do3d, key, n=n)
            ps+=ps1
            count+=1
    ps/=count
    ps*=np.power(freqs, -eExpo)
    ps/=np.mean(ps)
    plt.loglog(freqs, ps, label=key, color=color)
    n1+=1
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e-2,1.e2)
plt.legend()
tools.saveAndClear(pathSave + 'adjustedPspecSpheresSumTots.png', figNum=0)








#
