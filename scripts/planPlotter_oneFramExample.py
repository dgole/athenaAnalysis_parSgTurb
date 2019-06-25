#!/usr/bin/python
import numpy as np
import matplotlib as m
#m.use('Agg')
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
from mpl_toolkits.mplot3d import Axes3D
import numpy.polynomial.polynomial as poly
################################################################################
# paths and CL args
pathBase = str(sys.argv[1])
n        = int(sys.argv[2])
nb       = 1
color    = 'b'
pathPlan = pathBase + 'planOutput2/'
pathSave = pathBase + 'plots/planForPaper/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
doPlan = readerPlan.DataPlan(pathPlan, nStart=n, nTot=n+1, nPar=512*512*512)
################################################################################

# make hist
mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
mp = np.asarray(mp); dndmp = np.asarray(dndmp);
minMass = np.amin(mp); maxMass = np.amax(mp);
mp1, ngtm = readerPlan.getCumMassHist(doPlan, n)
nm = mp1.shape[0]

################################################################################
# Single PL MLE
p_mle, err_mle = readerPlan.get_p_mle(doPlan, n)
# Fit with spl
means_spl, errsPlus_spl, errsMinus_spl, maxLike_spl = readerPlan.bootstrap(mp1, readerPlan.fit_spl, 1, nb=nb)



################################################################################
# diff hist
min = 1.e20
preFacArr = np.arange(1.0,200.0,0.1)
for i in range(preFacArr.shape[0]):
	preFac   = preFacArr[i]
	model    = preFac * np.power(mp, -p_mle)
	logDndmp = np.log10(dndmp)
	logModel = np.log10(model)
	diff     = np.sum(np.square(logDndmp - logModel))
	if diff < min:
		min  = diff
		bestPreFac = preFac
		#print("new min found")
		#print(min)
		#print(bestPreFac)

plt.figure(0)
plt.loglog(mp, dndmp, color+'o', ms=2)
mp2 = np.logspace(-10.0, 10, num=500)
model = bestPreFac * np.power(mp2, -p_mle)
plt.loglog(mp2, model, linestyle='--', linewidth=1, color=color)
# plot labels etc.
plt.xlabel(r'$M_p$')
plt.ylabel(r'$dN/dM_p$')
plt.ylim(1.e-1,  1.e5)
plt.xlim(1e-3, 1.e1)
tools.saveAndClear(pathSave + "hist_differential_oneFrameExample_" + str(n) + ".png", figNum=0)


################################################################################
'''
# cumulative hist plot
# SPL
ngtm_spl = readerPlan.P_spl(mp1, means_spl)
plt.loglog(mp1, ngtm, color=(0,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, ngtm, color=(0,0,0,1.0), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, nm*ngtm_spl, color=(0,0,0,1.0), label='SPL', linestyle='--', linewidth=1)
# plot labels etc.
plt.xlabel(r'$M_p$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(1.e-1, 1.e3)
plt.xlim(1e-3, 1.e1)
tools.saveAndClear(pathSave + "hist_cumulative_oneFrameExample.png", figNum=0)
'''


















#
