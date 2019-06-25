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
nb       = int(sys.argv[3])
pathPlan = pathBase + 'planOutput2/'
pathSave = pathBase + 'plots/planBplDebugging/'
if not os.path.exists(pathSave): os.makedirs(pathSave)
plt.figure(0)
doPlan = readerPlan.DataPlan(pathPlan, nStart=n, nTot=n+1, nPar=512*512*512)
################################################################################
# make hist
mp, dndmp = readerPlan.getDiffMassHist(doPlan, n)
mp = np.asarray(mp); dndmp = np.asarray(dndmp);
mp1, ngtm = readerPlan.getCumMassHist(doPlan, n)
minMass = np.amin(mp1); maxMass = np.amax(mp1);
nm = mp1.shape[0]

print(mp.shape)
print(mp1.shape)

################################################################################
pathSave1 = None
# Fit with bpl
means_bpl, errsPlus_bpl, errsMinus_bpl, maxLike_bpl     = readerPlan.bootstrap(mp1, readerPlan.fit_bpl, 3, nb=nb, pathSave=pathSave1)

################################################################################
# information criteria analysis stuff
mp1_min = np.amin(mp1)

readerPlan.reportParams("BPL", ["alpha1", "alpha2", "x_br"],
		   means_bpl, errsPlus_bpl, errsMinus_bpl,
		   maxLike_bpl, 0.0, 0.0,
		   mp1_min)
print("########################################################")

################################################################################

# cumulative hist plot
logMinMass = np.log10(minMass); logMaxMass = np.log10(maxMass);
#fakeMp1 = np.logspace(logMinMass-0.01 , logMaxMass+0.01, num=mp1.shape[0])
fakeMp1 = np.logspace(logMinMass , logMaxMass, num=1000)

# BPL
'''
def p_bpl(masses, params):
	a1, a2, xb = params
	c0     = np.power( (1./a1)+((1./a2)-(1./a1))*np.power(np.exp(xb),-a1), -1)
	xArr   = convert_to_x(masses)
	ltxbr  = np.where(xArr <= xb, 1.0, 0.0)
	gtxbr  = np.where(xArr >  xb, 1.0, 0.0)
	ltxbr *= c0 * np.exp(-a1 * xArr)
	gtxbr *= c0 * np.exp((a2-a1)*xb - a2*xArr)
	return gtxbr + ltxbr
'''
#ngtm_bpl = readerPlan.P_bpl(fakeMp1, means_bpl)
p_bpl = readerPlan.p_bpl(fakeMp1, [-1.0, 1.0, 2.0])
P_bpl = readerPlan.p_to_P(fakeMp1, p_bpl)
plt.loglog(fakeMp1, p_bpl, 'ro', label='BPL', markersize=1)
print(fakeMp1)
print(P_bpl)
plt.loglog(fakeMp1, P_bpl, 'ko', label='BPL', markersize=1)


#fac = 1
#plt.loglog(mp1, fac*ngtm, color=(0,1,1,1), marker=".", linewidth=0, markersize=2)
#plt.loglog(mp1, fac*ngtm, color=(0,1,1,0.2), marker=".", linewidth=5, markersize=2)
#plt.loglog(fakeMp1, nm*fac*ngtm_bpl, color=(0,1,1,1), label='BPL', linestyle='--')
#arg = np.argmin(np.absolute(minMass*np.exp(means_bpl[2])-fakeMp1))
#plt.loglog(fakeMp1[arg], nm*fac*ngtm_bpl[arg], color=(0,1,1,1), marker="s", linewidth=0, markersize=7)


# plot labels etc.
plt.xlabel(r'$M_p [M_G]$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(1.e-3, 1.e4)
plt.xlim(1e-3, 1.e1)
plt.legend(prop={'size':10})
tools.saveAndClear(pathSave + "hist_cumulative_" + str(n) +".png", figNum=0)





#plt.figure(0)
#plt.loglog(mp, dndmp, 'ko', ms=2)
#p_spl = readerPlan.p_spl(mp, means_spl)

#dndmp_theo = readerPlan.pdf_to_dndmp(nm, mp, ngtm)
#plt.loglog(mp, dndmp_theo, linestyle='--', linewidth=1, color='k')
#mp2 = np.logspace(-10.0, 10, num=500)
#model = bestPreFac * np.power(mp2, -p_mle)
#plt.loglog(mp2, model, linestyle='--', linewidth=1, color=color)
# plot labels etc.
#plt.xlabel(r'$M_p$')
#plt.ylabel(r'$dN/dM_p$')
#plt.ylim(6.e-1,  1.e4)
#plt.xlim(1e-3, 1.e1)
#tools.saveAndClear(pathSave + "hist_differential_" + str(n) + ".png", figNum=0)






























#np.savetxt("./mp_control_253.csv", mp1)
#np.savetxt("./ngtm_control_253.csv", ngtm)

#np.savetxt("./mp_weak_291.csv", mp1)
#np.savetxt("./ngtm_weak_291.csv", ngtm)

#np.savetxt("./mp_moderate_367.csv", mp1)
#np.savetxt("./ngtm_moderate_367.csv", ngtm)





















#
