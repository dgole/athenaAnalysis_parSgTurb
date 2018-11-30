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
################################################################################
def addFiveThirdsToFig():
    for i in range(-1, 3):
        plt.loglog(freqs, 10**i*np.power(freqs, -5.0/3.0), color='tab:gray',
                   linestyle='--', linewidth=0.5)
################################################################################

title = 'keSpec'
vxpsk, vxpskx, vxpsky, vxpskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVx')
vypsk, vypskx, vypsky, vypskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVy')
vzpsk, vzpskx, vzpsky, vzpskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVz')
psk  = vxpsk  + vypsk  + vzpsk
pskx = vxpskx + vypskx + vzpskx
psky = vxpsky + vypsky + vzpsky
pskz = vxpskz + vypskz + vzpskz
norm = psk[1]
plt.figure(0)
plt.xlabel(r'$k$')
plt.ylabel('Power' + title)
plt.loglog(freqs, psk/norm,  'ko', label='all', markersize=ms)
plt.loglog(freqs, pskx/norm, 'ro', label='kx',  markersize=ms)
plt.loglog(freqs, psky/norm, 'go', label='ky',  markersize=ms)
plt.loglog(freqs, pskz/norm, 'bo', label='kz',  markersize=ms)
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + title + '.png', figNum=0)
plt.figure(1)
plt.loglog(freqs, psk/norm,  'bo', label=title, markersize=ms)

################################################################################

title = 'kePertSpec'
vxpsk, vxpskx, vxpsky, vxpskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvx')
print(vxpsk.shape, freqs.shape)
vypsk, vypskx, vypsky, vypskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvy')
vzpsk, vzpskx, vzpsky, vzpskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoDvz')
psk  = vxpsk  + vypsk  + vzpsk
pskx = vxpskx + vypskx + vzpskx
psky = vxpsky + vypsky + vzpsky
pskz = vxpskz + vypskz + vzpskz
norm = psk[1]
plt.figure(0)
plt.xlabel(r'$k$')
plt.ylabel('Power' + title)
plt.loglog(freqs, psk/norm,  'ko', label='all', markersize=ms)
plt.loglog(freqs, pskx/norm, 'ro', label='kx',  markersize=ms)
plt.loglog(freqs, psky/norm, 'go', label='ky',  markersize=ms)
plt.loglog(freqs, pskz/norm, 'bo', label='kz',  markersize=ms)
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + title + '.png', figNum=0)
plt.figure(1)
plt.loglog(freqs, psk/norm,  'ro', label=title, markersize=ms)

################################################################################

plt.figure(1)
plt.xlabel(r'$k$')
plt.ylabel('Power')
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + 'compare' + '.png', figNum=1)

















#
