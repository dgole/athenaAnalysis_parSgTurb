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
nStart   = int(sys.argv[2])
nEnd     = int(sys.argv[3])
path3d   = pathBase + '3d/'
pathSave = pathBase + 'plots/pspec/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
do3d     = reader3d.Data3d(path3d)
plt.figure(0)
################################################################################
vExpo = -1.833
eExpo = 2.0*(vExpo)+2.0
print("velocity spectrum PL exponent is " + str(vExpo))
print("KE power spectrum PL exponent is " + str(eExpo))
################################################################################

def calcPs(do, key, n):
	data       = do.get3d(key, n)
	data       = data[:,:,(do.nz//2-10):(do.nz//2)+10]
	freqs      = [np.fft.fftfreq(shape, d=do.dx) for shape in data.shape]
	nFreqs     = [freqs[0].shape[0]//2,freqs[1].shape[0]//2,freqs[2].shape[0]//2]
	fft        = np.fft.fftn(data*np.power(do.dx,3))[:nFreqs[0],:nFreqs[1],:nFreqs[2]]
	ps         = np.square(np.absolute(fft))
	return ps, freqs

def psProfile(do, key, n, norm=False):
	ps, freqs = calcPs(do, key, n)
	dk        = freqs[0][1]-freqs[0][0]
	#psk       = np.zeros(int(ps.shape[0]*np.sqrt(3))+1)
	#freqs     = np.arange(0,dk*len(psk),dk)
	psk       = np.zeros(ps.shape[0])
	freqs     = np.arange(0,dk*len(psk),dk)
	count     = np.zeros_like(psk)
	for i in range(ps.shape[0]):
		for j in range(ps.shape[1]):
			for k in range(ps.shape[2]):
				dist          = np.sqrt(i*i+j*j+k*k)
				index         = int(np.floor(dist))
				if index < freqs.shape[0]:
					psk[index]   += ps[i,j,k]
					count[index] += 1
	if norm==True: psk /= count
	return psk, freqs

def psProfileMean(do, key, nStart=None, nEnd=None, norm=False):
	#if nStart is None: nStart = do.nt // 2
	if nStart is None: nStart = do.nt - 10*int(1.0/do.dt)
	if nEnd   is None: nEnd   = do.nt
	pskList = []
	count   = 0
	for n in range(nStart, nEnd):
		if np.mod(do.t[n], do.shearTime)<1.e-8:
			psk, freqs = psProfile(do, key, n, norm=norm)
			pskList.append(psk)
			count += 1
	psk = np.asarray(np.sum(pskList,  axis=0))/float(count)
	return psk, freqs

#calcPs(do3d, 'vx', 6)
#psProfile(do3d, 'vx', 6)
#psProfileMean(do3d, 'vx', 6, 11)

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
