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
kStart   = int(sys.argv[4])
kEnd     = int(sys.argv[5])
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

def calcPs_2D(do, data2d):
	#data       = do.get3d(key, n)
	#data2d     = data[:,:,k]
	freqs      = np.fft.fftfreq(data2d.shape[0], d=do.dx)
	nFreqs     = freqs.shape[0]//2
	freqs      = freqs[:nFreqs]
	fft        = np.fft.fftn(data2d*np.power(do.dx,2))[:nFreqs,:nFreqs]
	ps         = np.square(np.absolute(fft))
	return ps, freqs

def psProfile_2D(do, data2d, norm=False):
	ps, freqs = calcPs_2D(do, data2d)
	dk        = freqs[1]-freqs[0]
	freqsNew  = np.arange(0, dk*ps.shape[0], dk)
	psk       = np.zeros_like(freqsNew)
	count     = np.zeros_like(freqsNew)
	for i in range(ps.shape[0]):
		for j in range(ps.shape[1]):
			totFreq		  = np.sqrt(np.square(freqs[i]) + np.square(freqs[j]))
			index         = np.argmin(np.absolute(freqsNew-totFreq))
			if index < freqsNew.shape[0]:
				#print(totFreq, index)
				psk[index]   += ps[i,j]
				count[index] += 1
	if norm==True: psk /= count
	return psk, freqsNew

def psProfile_ztAvg(do, key, nStart, nEnd, kStart, kEnd, norm=False):
	pskList = []
	count   = 0
	for n in range(nStart, nEnd):
		if np.mod(do.t[n], do.shearTime)<1.e-8:
			data = do.get3d(key, n)
			for k in range(kStart, kEnd):
				data2d = data[:,:,k]
				psk, freqs = psProfile_2D(do, data2d, norm=norm)
				pskList.append(psk)
				count += 1
	psk = np.asarray(np.sum(pskList,  axis=0))/float(count)
	return psk, freqs

#calcPs_2D(do3d, 'rootRhoDvx', 10, 32)
#psProfile_2D(do3d, 'rootRhoDvx', 10, 32)
#psProfileMean_2D(do3d, 'rootRhoDvx', 9, 11, 30, 34)

psk_vx, freqs = psProfile_ztAvg(do3d, 'rootRhoDvx', nStart, nEnd, kStart, kEnd)
psk_vy, freqs = psProfile_ztAvg(do3d, 'rootRhoDvy', nStart, nEnd, kStart, kEnd)
psk_vz, freqs = psProfile_ztAvg(do3d, 'rootRhoDvz', nStart, nEnd, kStart, kEnd)
psk  = psk_vx  + psk_vy  + psk_vz
psk *= np.power(freqs, -eExpo)
psk /= np.mean(psk)
print(freqs)
print(psk)
plt.loglog(freqs, psk)
plt.xlabel(r'$|\mathbf{k}|$')
plt.ylabel('Power')
plt.ylim(1.e-2,1.e2)
plt.xlim(freqs[1],freqs[-1])
tools.saveAndClear(pathSave + 'adjustedPspecSpheresSumPerts.png', figNum=0)









#
