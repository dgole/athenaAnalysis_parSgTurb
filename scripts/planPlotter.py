#!/usr/bin/python
import numpy as np
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import os
import math
import sys
sys.path.append('../python')
import athenaReader1d as reader1d
import athenaReaderPhst as readerPhst
import athenaTools as tools
import planOutputReader as readerPlan
from matplotlib.backends.backend_pdf import PdfPages
################################################################################
# paths
pathBase = str(sys.argv[1])
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/plan/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
################################################################################
# 1    2    3      4    5    6
# npar mass r_hill xcom ycom zcom
doPlan = readerPlan.DataPlan(pathPlan)
################################################################################
plt.figure(0)
plt.plot(doPlan.time, doPlan.nClumps, 'ko')
plt.xlabel(r'$t \Omega$')
plt.ylabel(r'$N_{clumps}$')
tools.saveAndClear(pathSave + "nClumps_timeEvo.png", figNum=0)
################################################################################
plt.figure(0)
n = len(doPlan.time)-1
bins, hist = readerPlan.getCumMassHist(doPlan, n)
plt.loglog(bins, hist, 'ko', label='at end')
n = np.argmax(doPlan.nClumps)
bins, hist = readerPlan.getCumMassHist(doPlan, n)
plt.loglog(bins, hist, 'bo', label='at max Nclumps')
plt.xlabel(r'$M_P$')
plt.ylabel(r'$N(>M_P)$')
plt.legend()
tools.saveAndClear(pathSave + "hist.png", figNum=0)
################################################################################
plt.figure(0)
for n in range(0,10000):
	if doPlan.nClumps[n]>0:
		nStart = n
		break;
nEnd = nStart + 50
for n in range(nStart,nEnd,5):
	bins, hist = readerPlan.getCumMassHist(doPlan, n)
	rgb = tools.getColor(n, nStart, nEnd)
	plt.loglog(bins, hist, color=rgb, label='n='+str(n))
plt.xlabel(r'$M_P$')
plt.ylabel(r'$N(>M_P)$')
plt.legend()
tools.saveAndClear(pathSave + "hist_timeEvo.png", figNum=0)
################################################################################









#
