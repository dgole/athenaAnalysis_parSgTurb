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

title = 'KE'
psk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'KE')
plt.figure(0)
plt.xlabel(r'$k$')
plt.ylabel('Power ' + title)
norm = psk[1]
plt.loglog(freqs, psk/norm,  'ko', label='all', markersize=ms)
plt.loglog(freqs, pskx/norm, 'ro', label='kx', markersize=ms)
plt.loglog(freqs, psky/norm, 'go', label='ky', markersize=ms)
plt.loglog(freqs, pskz/norm, 'bo', label='kz', markersize=ms)
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + title + '.png', figNum=0)
plt.figure(1)
plt.loglog(freqs, psk/norm,  'ko', label=title, markersize=ms)

################################################################################

title = 'dKE'
psk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'dKE')
plt.figure(0)
plt.xlabel(r'$k$')
plt.ylabel('Power ' + title)
norm = psk[1]
plt.loglog(freqs, psk/norm,  'ko', label='all', markersize=ms)
plt.loglog(freqs, pskx/norm, 'ro', label='kx', markersize=ms)
plt.loglog(freqs, psky/norm, 'go', label='ky', markersize=ms)
plt.loglog(freqs, pskz/norm, 'bo', label='kz', markersize=ms)
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + title + '.png', figNum=0)
plt.figure(1)
plt.loglog(freqs, psk/norm,  'go', label=title, markersize=ms)

################################################################################

title = 'root of KE comps FFTs summed'
vxpsk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'rootKEx')
vypsk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'rootKEy')
vzpsk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'rootKEz')
psk = vxpsk + vypsk + vzpsk
norm = psk[1]
plt.figure(0)
plt.xlabel(r'$k$')
plt.ylabel('Power' + title)
plt.loglog(freqs, psk/norm, 'ko', markersize=ms)
addFiveThirdsToFig()
tools.saveAndClear(pathSave + title + '.png', figNum=0)
plt.figure(1)
plt.loglog(freqs, psk/norm,  'bo', label=title, markersize=ms)

################################################################################

title = 'root(0.5*rho)*v FFTs summed'
vxpsk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVx')
vypsk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVy')
vzpsk, pskx, psky, pskz, freqs = reader3d.psProfileMean(do3d, 'rootRhoVz')
psk = vxpsk + vypsk + vzpsk
norm = psk[1]
plt.figure(0)
plt.xlabel(r'$k$')
plt.ylabel('Power' + title)
plt.loglog(freqs, psk/norm, 'ko', markersize=ms)
addFiveThirdsToFig()
tools.saveAndClear(pathSave + title + '.png', figNum=0)
plt.figure(1)
plt.loglog(freqs, psk/norm,  'ro', label=title, markersize=ms)

################################################################################

plt.figure(1)
plt.xlabel(r'$k$')
plt.ylabel('Power')
addFiveThirdsToFig()
plt.legend()
tools.saveAndClear(pathSave + 'comparison' + '.png', figNum=1)

















#
