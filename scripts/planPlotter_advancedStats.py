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
pathSave = pathBase + 'plots/plan_fitTesting/'
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
pathSave1 = pathSave
# Single PL MLE
p_mle, err_mle = readerPlan.get_p_mle(doPlan, n)
# Fit the histogram
p_fit, err_fit, c_fit = readerPlan.get_p_fit(doPlan, n)
# Fit with spl
means_spl, errsPlus_spl, errsMinus_spl, maxLike_spl = readerPlan.bootstrap(mp1, readerPlan.fit_spl, 1, nb=nb, pathSave=pathSave1)
# Fit with stpl
means_stpl, errsPlus_stpl, errsMinus_stpl, maxLike_stpl = readerPlan.bootstrap(mp1, readerPlan.fit_stpl, 2, nb=nb, pathSave=pathSave1)
# Fit with vtpl
means_vtpl, errsPlus_vtpl, errsMinus_vtpl, maxLike_vtpl =	readerPlan.bootstrap(mp1, readerPlan.fit_vtpl, 3, nb=nb, pathSave=pathSave1)
# Fit with bcpl
means_bcpl, errsPlus_bcpl, errsMinus_bcpl, maxLike_bcpl =	readerPlan.bootstrap(mp1, readerPlan.fit_bcpl, 3, nb=nb, pathSave=pathSave1)
# Fit with tpl
means_tpl, errsPlus_tpl, errsMinus_tpl, maxLike_tpl     = readerPlan.bootstrap(mp1, readerPlan.fit_tpl, 2, nb=nb, pathSave=pathSave1)
# Fit with bpl
means_bpl, errsPlus_bpl, errsMinus_bpl, maxLike_bpl     = readerPlan.bootstrap(mp1, readerPlan.fit_bpl, 3, nb=nb, pathSave=pathSave1)

################################################################################
# information criteria analysis stuff

N = mp.shape[0]

K = 1
bic_spl = readerPlan.BIC(K, N, maxLike_spl)
aic_spl = readerPlan.AIC(K, N, maxLike_spl)
K = 2
bic_stpl = readerPlan.BIC(K, N, maxLike_stpl)
aic_stpl = readerPlan.AIC(K, N, maxLike_stpl)
K = 3
bic_vtpl = readerPlan.BIC(K, N, maxLike_vtpl)
aic_vtpl = readerPlan.AIC(K, N, maxLike_vtpl)
K = 3
bic_bcpl = readerPlan.BIC(K, N, maxLike_bcpl)
aic_bcpl = readerPlan.AIC(K, N, maxLike_bcpl)
K = 2
bic_tpl = readerPlan.BIC(K, N, maxLike_tpl)
aic_tpl = readerPlan.AIC(K, N, maxLike_tpl)
K = 3
bic_bpl = readerPlan.BIC(K, N, maxLike_bpl)
aic_bpl = readerPlan.AIC(K, N, maxLike_bpl)

bic_min = min(bic_spl, bic_stpl, bic_vtpl, bic_bcpl, bic_tpl, bic_bpl)
aic_min = min(aic_spl, aic_stpl, aic_vtpl, aic_bcpl, aic_tpl, aic_bpl)

dbic_spl  = bic_spl  - bic_min
daic_spl  = aic_spl  - aic_min
dbic_stpl = bic_stpl - bic_min
daic_stpl = aic_stpl - aic_min
dbic_vtpl = bic_vtpl - bic_min
daic_vtpl = aic_vtpl - aic_min
dbic_bcpl = bic_bcpl - bic_min
daic_bcpl = aic_bcpl - aic_min
dbic_tpl  = bic_tpl  - bic_min
daic_tpl  = aic_tpl  - aic_min
dbic_bpl  = bic_bpl  - bic_min
daic_bpl  = aic_bpl  - aic_min

################################################################################
# print out all fits / params / errors in a nice format
#print("########################################################")
#print("Nclumps: " + str(doPlan.nClumpsList[n]))
mp1_min = np.amin(mp1)
print("########################################################")
print("Single PL MLE Slope")
print(str(np.round(p_mle,  3)) + " pm "
	+ str(np.round(err_mle,3)))
readerPlan.reportParams("SPL", ["alpha"],
		   means_spl, errsPlus_spl, errsMinus_spl,
		   maxLike_spl, dbic_spl, daic_spl,
		   mp1_min)
readerPlan.reportParams("STPL", ["alpha", "x_exp"],
		   means_stpl, errsPlus_stpl, errsMinus_stpl,
		   maxLike_stpl, dbic_stpl, daic_stpl,
		   mp1_min)
readerPlan.reportParams("VTPL", ["alpha", "beta", "x_exp"],
		   means_vtpl, errsPlus_vtpl, errsMinus_vtpl,
		   maxLike_vtpl, dbic_vtpl, daic_vtpl,
		   mp1_min)
readerPlan.reportParams("BCPL", ["alpha1", "alpha2", "x_br"],
		   means_bcpl, errsPlus_bcpl, errsMinus_bcpl,
		   maxLike_bcpl, dbic_bcpl, daic_bcpl,
		   mp1_min)
readerPlan.reportParams("TPL", ["alpha", "x_tr"],
		   means_tpl, errsPlus_tpl, errsMinus_tpl,
		   maxLike_tpl, dbic_tpl, daic_tpl,
		   mp1_min)
readerPlan.reportParams("BPL", ["alpha1", "alpha2", "x_br"],
		   means_bpl, errsPlus_bpl, errsMinus_bpl,
		   maxLike_bpl, dbic_bpl, daic_bpl,
		   mp1_min)
print("########################################################")
print("mean mass: " + str(np.mean(mp1)))
print("min mass: " + str(np.mean(mp1_min)))


################################################################################
# cumulative hist plot
logMinMass = np.log10(minMass); logMaxMass = np.log10(maxMass);
#fakeMp1 = np.logspace(logMinMass-0.01 , logMaxMass+0.01, num=mp1.shape[0])
fakeMp1 = np.logspace(logMinMass , logMaxMass, num=1000)
fakeMp2 = np.logspace(logMinMass , logMaxMass*100.0, num=1000)
#print(fakeMp1)
#print(fakeMp1.shape)
# SPL
ngtm_spl = readerPlan.P_spl(fakeMp1, means_spl)
fac = 1.e5
plt.loglog(mp1, fac*ngtm, color=(0,0,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_spl, color=(0,0,0,1), label='SPL', linestyle='--')
# STPL
ngtm_stpl = readerPlan.P_stpl(fakeMp1, means_stpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(1,0,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(1,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_stpl, color=(1,0,0,1), label='STPL', linestyle='--')
arg = np.argmin(np.absolute(minMass*np.exp(means_stpl[1])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_stpl[arg], color=(1,0,0,1), marker="s", linewidth=0, markersize=7)
# VTPL
ngtm_vtpl = readerPlan.P_vtpl(fakeMp1, means_vtpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(0,1,0,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,1,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_vtpl, color=(0,1,0,1), label='VTPL', linestyle='--')
#arg = np.argmin(np.absolute(minMass*np.exp(means_vtpl[2])-fakeMp1))
#plt.loglog(fakeMp1[arg], nm*fac*ngtm_vtpl[arg], color=(0,1,0,1), marker="s", linewidth=0, markersize=7)
# BCPL
ngtm_bcpl = readerPlan.P_bcpl(fakeMp1, means_bcpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(0,0,1,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,0,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_bcpl, color=(0,0,1,1.0), label='BCPL', linestyle='--')
arg = np.argmin(np.absolute(minMass*np.exp(means_bcpl[2])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_stpl[arg], color=(0,0,1,1), marker="s", linewidth=0, markersize=7)
# TPL
ngtm_tpl = readerPlan.P_tpl(fakeMp1, means_tpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(1,0,1,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(1,0,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_tpl, color=(1,0,1,1), label='TPL', linestyle='--')
# BPL
ngtm_bpl = readerPlan.P_bpl(fakeMp1, means_bpl)
fac /=10
plt.loglog(mp1, fac*ngtm, color=(0,1,1,1), marker=".", linewidth=0, markersize=2)
plt.loglog(mp1, fac*ngtm, color=(0,1,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(fakeMp1, nm*fac*ngtm_bpl, color=(0,1,1,1), label='BPL', linestyle='--')
arg = np.argmin(np.absolute(minMass*np.exp(means_bpl[2])-fakeMp1))
plt.loglog(fakeMp1[arg], nm*fac*ngtm_bpl[arg], color=(0,1,1,1), marker="s", linewidth=0, markersize=7)
# plot labels etc.
plt.xlabel(r'$M_p [M_G]$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(6.e-1, 1.e8)
plt.xlim(1e-3, 1.e1)
plt.legend(prop={'size':10}, loc='best')
tools.saveAndClear(pathSave + "hist_cumulative_" + str(n) +".png", figNum=0)


































#np.savetxt("./mp_control_253.csv", mp1)
#np.savetxt("./ngtm_control_253.csv", ngtm)

#np.savetxt("./mp_weak_291.csv", mp1)
#np.savetxt("./ngtm_weak_291.csv", ngtm)

#np.savetxt("./mp_moderate_367.csv", mp1)
#np.savetxt("./ngtm_moderate_367.csv", ngtm)





















#
