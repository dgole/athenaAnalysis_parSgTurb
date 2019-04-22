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
pathBase = '../../data/kspaceTest/run141/'
pathSave = pathBase + 'plots/turbForcingDump/'
vPertSlope = 5./6.
fileName = 'out_dump.txt'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
nx=32;ny=32;nz=32;
nGz = 5
setupPlots = False
amplPlots  = False
################################################################################
def getData(pathBase, fileName, lookupStr, gzOffset, ind):
    inFile   = open(pathBase+fileName, 'r')
    data1 = np.zeros([nx,ny,nz])
    count=0
    for line in inFile:
        if line[0:len(lookupStr)] == lookupStr:
            split=line.split()
            #if count%1.e0==0:
                #print(split)
                #print(count/(nx*ny*nz))
            i=int(split[1])-gzOffset
            j=int(split[2])-gzOffset
            k=int(split[3])-gzOffset
            data1[i,j,k]=float(split[4+ind])
            count+=1
            if count == nx*ny*nz: break
    return data1
################################################################################
def addFiveThirdsToFig():
    thisSlope = 2*(vPertSlope-1.)
    print(-thisSlope)
    for i in np.arange(-20, 20, 0.5):
        plt.loglog(freqs, 10**i*np.power(freqs, -thisSlope), color='tab:gray',
                   linestyle='--', linewidth=0.5)
################################################################################
if setupPlots == True:
    data = getData(pathBase, fileName, 'q1q2q3_rand', nGz, 0)
    plt.hist(data.flatten(), bins=100)
    tools.saveAndClear(pathSave + 'q1_rand.png', figNum=0)
    data = getData(pathBase, fileName, 'q1q2q3_rand', nGz, 1)
    plt.hist(data.flatten(), bins=100)
    tools.saveAndClear(pathSave + 'q2_rand.png', figNum=0)
    data = getData(pathBase, fileName, 'q1q2q3_rand', nGz, 2)
    plt.hist(data.flatten(), bins=100)
    tools.saveAndClear(pathSave + 'q3_rand.png', figNum=0)
    data = getData(pathBase, fileName, 'q3', nGz, 0)
    plt.hist(data.flatten(), bins=100)
    plt.ylim(0,200)
    tools.saveAndClear(pathSave + 'q3.png', figNum=0)
    data = getData(pathBase, fileName, 'q3', nGz, 1)
    data = getData(pathBase, fileName, 'ampl_before', 0, 0)
    plt.hist(data.flatten(), bins=100)
    tools.saveAndClear(pathSave + 'ampl0_before.png', figNum=0)
    data = getData(pathBase, fileName, 'ampl_before', 0, 1)
    plt.hist(data.flatten(), bins=100)
    tools.saveAndClear(pathSave + 'ampl1_before.png', figNum=0)
################################################################################
if amplPlots == True:
    extent = [0, nx//2, 0, ny//2]
    aspect  = 1.0
    cmapType = 'coolwarm'
    data = getData(pathBase, fileName, 'ampl_after', 0, 1)
    data2d = data[:nx//2, :ny//2, 0]
    plotData = np.transpose(np.fliplr(data2d))
    plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
    plt.colorbar()
    plt.tight_layout()
    tools.saveAndClear(pathSave + 'ampl0_colorMap.png', figNum=0)
    data = getData(pathBase, fileName, 'ampl_after', 0, 2)
    data2d = data[:, :, 0]
    plotData = np.transpose(np.fliplr(data2d))
    plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
    plt.colorbar()
    plt.tight_layout()
    tools.saveAndClear(pathSave + 'ampl1_colorMap.png', figNum=0)
################################################################################
data = getData(pathBase, fileName, 'perturb', 5, 0)
plt.plot(data[:,ny//2,nz//2])
tools.saveAndClear(pathSave + 'perturb_vx.png', figNum=0)
data = getData(pathBase, fileName, 'perturb', 5, 1)
plt.plot(data[:,ny//2,nz//2])
tools.saveAndClear(pathSave + 'perturb_vy.png', figNum=0)
data = getData(pathBase, fileName, 'perturb', 5, 2)
plt.plot(data[:,ny//2,nz//2])
tools.saveAndClear(pathSave + 'perturb_vz.png', figNum=0)
################################################################################
freqs  = np.fft.fftfreq(data.shape[0], d=0.2/float(nx))
nFreqs = freqs.shape[0]//2
freqs  = freqs[:nFreqs]
ps     = np.zeros(nFreqs)
count  = 0
for j in range(0,ny):
    fft    = np.fft.fftn(data[:,j,0])[:nFreqs]
    ps    += np.square(np.absolute(fft))
    count += 1
plt.loglog(freqs[3:], ps[3:]/float(count), 'ko')
addFiveThirdsToFig()
plt.ylim(1.e-7,1.e-2)
plt.xlim(freqs[3],freqs[-1])
tools.saveAndClear(pathSave + '1testPlot.png', figNum=0)



################################################################################
freqs      = np.fft.fftfreq(data.shape[0], d=0.2/float(nx))
nFreqs     = freqs.shape[0]//2
freqs      = freqs[:nFreqs]
fft        = np.fft.fftn(data)[:nFreqs,:nFreqs,:nFreqs]
ps         = np.square(np.absolute(fft))
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
# ps colorMap at k=0
extent = [0, nx//2, 0, ny//2]
aspect  = 1.0
cmapType = 'coolwarm'
data2d = ps[:, :, 0]
plotData = np.transpose(np.fliplr(data2d))
plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
plt.colorbar()
plt.tight_layout()
tools.saveAndClear(pathSave + 'pspec_colorMap.png', figNum=0)
# spherical shells
#plt.loglog(freqs, psk/count)
plt.loglog(freqs[3:], psk[3:], 'ko')
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e0,1.e4)
plt.xlim(freqs[3],freqs[-1])
addFiveThirdsToFig()
tools.saveAndClear(pathSave + 'pspecShells.png', figNum=0)
# slice
data2d = data[:, :, nz//2]
extent = [-0.2, 0.2, -0.2, 0.2]
aspect  = 1.0
plotData = np.transpose(np.fliplr(data2d))
cmapType = 'coolwarm'
plt.imshow(plotData, extent=extent, aspect=aspect, cmap=plt.get_cmap(cmapType))
plt.colorbar()
plt.tight_layout()
tools.saveAndClear(pathSave + 'slice.png', figNum=0)
################################################################################








#
