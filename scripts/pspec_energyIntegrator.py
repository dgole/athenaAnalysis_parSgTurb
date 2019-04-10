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
powerSum  = []
powerSum2 = []
keSum     = []
dVol = do3d.dx*do3d.dx*do3d.dx
key = 'rootRhoVx'
for n in range(2,11,2):
    data       = do3d.get3d(key, n)
    totKez     = np.sum(np.square(data)*dVol)
    keSum.append(np.sum(totKez))
    # linear FFT
    freqs      = np.fft.fftfreq(data.shape[0], d=do3d.dx)
    dk3        = np.power(freqs[1]-freqs[0],3)
    fft        = np.fft.fftn(data*dVol)
    ps3d       = np.square(np.absolute(fft))*dk3
    totPower   = np.sum(ps3d)
    powerSum.append(totPower)
    # spherical FFT
    psk       = np.zeros(int(ps3d.shape[0]*np.sqrt(3.0))+1)
    count     = np.zeros_like(psk)
    for i in range(ps3d.shape[0]):
        for j in range(ps3d.shape[1]):
            for k in range(ps3d.shape[2]):
                dist        = np.sqrt(i*i+j*j+k*k)
                index       = int(np.floor(dist))
                psk[index] += ps3d[i,j,k]
                count[index] += 1
    totPower2 = np.sum(psk)
    powerSum2.append(totPower2)
    psk /= count

print(keSum)
print(powerSum)
print(powerSum2)
################################################################################






















#
