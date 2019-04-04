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
name = str(sys.argv[1])
pathBase = '../../data/prodRuns/'
runName = 'run310'
pathSave = pathBase + runName + "/plots/pspec_test/"
if not os.path.exists(pathSave): os.makedirs(pathSave)
path3d = pathBase + runName + '/3d/'
do3d = reader3d.Data3d(path3d)
plt.figure(0)
def addFiveThirdsToFig():
    for i in np.arange(-10, 20, 0.5):
        plt.loglog(freqs, 10**i*np.power(freqs, -5.0/3.0), color='tab:gray',
                   linestyle='--', linewidth=0.5)
n=6
newVx = np.zeros_like(do3d.get3d('vx', n))
################################################################################
# whie noise
if name == 'whiteNoise':
    newVx = np.random.random_sample(newVx.shape)-0.5
################################################################################
# gaussian noise
if name == 'gasussianNoise':
    newVx = np.random.normal(0.0, 0.5, size=newVx.shape)
################################################################################
# sine wave in all directions
if name == '2modes':
    lBox        = 0.2
    kBox        = (1./lBox)*2.*3.14159
    def vxFunc(amp, modesPerBox, phase, x,y,z):
        vx = amp*np.sin(kBox*modesPerBox*x+phase)*np.sin(kBox*modesPerBox*y+phase)*np.sin(kBox*modesPerBox*z+phase)
        return vx
    for modesPerBox in [2.0, 10.0]:
        print(modesPerBox)
        #phase = 3.1415*2.*np.random.random_sample()
        phase = 0.
        #amp = np.sqrt(np.power(modesPerBox,-5./3.))
        amp = 1.0
        for i in range(do3d.nx):
            for j in range(do3d.ny):
                for k in range(do3d.nz):
                    newVx[i,j,k] += vxFunc(amp, modesPerBox, phase, do3d.x[i],do3d.y[j],do3d.z[k])
    newVx+=np.random.random_sample(newVx.shape)-0.5

################################################################################
# many sine waves
if name == 'fiveThirds':
    lBox        = 0.2
    kBox        = (1./lBox)*2.*3.14159
    def vxFunc(amp, modesPerBox, phase, x,y,z):
        vx = amp*np.sin(kBox*modesPerBox*x+phase)*np.sin(kBox*modesPerBox*y+phase)*np.sin(kBox*modesPerBox*z+phase)
        return vx
    for modesPerBox in np.arange(1.,32.,1.0):
        print(modesPerBox)
        #phase = 3.1415*2.*np.random.random_sample()
        phase = 0.
        amp = 0.1*np.sqrt(np.power(modesPerBox,-5./3.))
        for i in range(do3d.nx):
            for j in range(do3d.ny):
                for k in range(do3d.nz):
                    newVx[i,j,k] += vxFunc(amp, modesPerBox, phase, do3d.x[i],do3d.y[j],do3d.z[k])

################################################################################
data = newVx
freqs      = np.fft.fftfreq(data.shape[0], d=do3d.dx)
nFreqs     = freqs.shape[0]//2
freqs      = freqs[:nFreqs]
fft        = np.fft.fftn(data)[:nFreqs,:nFreqs,:nFreqs]
ps         = np.square(np.absolute(fft))
#return ps, freqs
psk       = np.zeros(ps.shape[0])
count     = np.zeros_like(psk)
for i in range(ps.shape[0]):
    for j in range(ps.shape[1]):
        for k in range(ps.shape[2]):
            dist          = np.sqrt(i*i+j*j+k*k)
            index         = int(np.floor(dist))
            if index < freqs.shape[0]:
                psk[index]   += ps[i,j,k]
                count[index] += 1
#psk /= count
#return psk, freqs
################################################################################
# ps along 1 axis
plt.loglog(freqs, np.sum(ps, axis=(1,2)), color='r')
plt.loglog(freqs, np.sum(ps, axis=(0,2)), color='g')
plt.loglog(freqs, np.sum(ps, axis=(0,1)), color='b')
plt.ylim(1.e4,1.e10)
addFiveThirdsToFig()
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.title(np.sum(ps))
tools.saveAndClear(pathSave + name + '_pspec.png', figNum=0)
################################################################################
# spherical shells not normalized
plt.loglog(freqs, psk)
plt.ylim(1.e4,1.e10)
addFiveThirdsToFig()
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.title(np.sum(psk))
tools.saveAndClear(pathSave + name + '_pspecShells.png', figNum=0)
################################################################################
# spherical shells
plt.loglog(freqs, psk/count)
plt.ylim(1.e4,1.e10)
addFiveThirdsToFig()
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.title(np.sum(psk))
tools.saveAndClear(pathSave + name + '_pspecShellsNorm.png', figNum=0)
################################################################################
data2d = newVx[:, :, do3d.nz//2]
extent = [-do3d.xmax, do3d.xmax, -do3d.ymax, do3d.ymax]
aspect  = 1.0
plotData = np.transpose(np.fliplr(data2d))
cmapType = 'coolwarm'
plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
plt.colorbar()
plt.tight_layout()
tools.saveAndClear(pathSave + name + '_slice.png', figNum=0)
################################################################################














#
