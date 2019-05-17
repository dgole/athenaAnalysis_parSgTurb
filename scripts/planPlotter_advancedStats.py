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
pathPlan = pathBase + 'planOutput/'
pathSave = pathBase + 'plots/plan/'
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
# Fit the histogram
p_fit, err_fit, c_fit = readerPlan.get_p_fit(doPlan, n)
# Fit with spl
means_spl, errsPlus_spl, errsMinus_spl, maxLike_spl = readerPlan.bootstrap(mp1, readerPlan.fit_spl, 1, nb=nb)
# Fit with stpl
means_stpl, errsPlus_stpl, errsMinus_stpl, maxLike_stpl = readerPlan.bootstrap(mp1, readerPlan.fit_stpl, 2, nb=nb)
# Fit with vtpl
means_vtpl, errsPlus_vtpl, errsMinus_vtpl, maxLike_vtpl =	readerPlan.bootstrap(mp1, readerPlan.fit_vtpl, 3, nb=nb)
# Fit with bcpl
means_bcpl, errsPlus_bcpl, errsMinus_bcpl, maxLike_bcpl =	readerPlan.bootstrap(mp1, readerPlan.fit_bcpl, 3, nb=nb)
# Fit with tpl
means_tpl, errsPlus_tpl, errsMinus_tpl, maxLike_tpl     = readerPlan.bootstrap(mp1, readerPlan.fit_tpl, 2, nb=nb)
# Fit with bpl
means_bpl, errsPlus_bpl, errsMinus_bpl, maxLike_bpl     = readerPlan.bootstrap(mp1, readerPlan.fit_bpl, 3, nb=nb)

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
def report(name, paramNames, means, errsPlus, errsMinus, maxLike, dbic, daic):
	nParams = len(paramNames)
	print("########################################################")
	print(name)
	if nParams == 1:
		print(paramNames[0].ljust(8) + " = "
			+ str(np.round(means,3)) + " p "
			+ str(np.round(errsPlus[0],3)) + " m "
			+ str(np.round(errsMinus[0],3)))
	else:
		for n in range(nParams):
			print(paramNames[n].ljust(8) + " = "
				+ str(np.round(means[n],3)) + " p "
				+ str(np.round(errsPlus[n],3)) + " m "
				+ str(np.round(errsMinus[n],3)))
	print("ln(like) = " + str(np.round(maxLike,3)))
	print("DBIC     = " + str(np.round(dbic,3)))
	print("DAIC     = " + str(np.round(daic,3)))


print("########################################################")
print("Nclumps: " + str(doPlan.nClumpsList[n]))
print("########################################################")
print("Single PL MLE Slope")
print(str(np.round(p_mle,  3)) + " pm "
	+ str(np.round(err_mle,3)))
print("########################################################")
print("Single PL Direct Fit")
print(str(np.round(p_fit,  3)) + " pm "
	+ str(np.round(err_fit,3)))
report("SPL", ["alpha"],
 	   means_spl, errsPlus_spl, errsMinus_spl,
	   maxLike_spl, dbic_spl, daic_spl)
report("STPL", ["alpha", "x_exp"],
 	   means_stpl, errsPlus_stpl, errsMinus_stpl,
	   maxLike_stpl, dbic_stpl, daic_stpl)
report("VTPL", ["alpha", "beta", "x_exp"],
 	   means_vtpl, errsPlus_vtpl, errsMinus_vtpl,
	   maxLike_vtpl, dbic_vtpl, daic_vtpl)
report("BCPL", ["alpha1", "alpha2", "x_br"],
 	   means_bcpl, errsPlus_bcpl, errsMinus_bcpl,
	   maxLike_bcpl, dbic_bcpl, daic_bcpl)
report("TPL", ["alpha", "x_tr"],
 	   means_tpl, errsPlus_tpl, errsMinus_tpl,
	   maxLike_tpl, dbic_tpl, daic_tpl)
report("BPL", ["alpha1", "alpha2", "x_br"],
 	   means_bpl, errsPlus_bpl, errsMinus_bpl,
	   maxLike_bpl, dbic_bpl, daic_bpl)

################################################################################
# cumulative hist plot
# SPL
ngtm_spl = readerPlan.P_spl(mp1, means_spl)
fac = 1.0
plt.loglog(mp1, fac*ngtm, color=(0,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, nm*fac*ngtm_spl, color=(0,0,0,1), label='SPL', linestyle='--')
# STPL
ngtm_stpl = readerPlan.P_stpl(mp1, means_stpl)
fac *=10
plt.loglog(mp1, fac*ngtm, color=(1,0,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, nm*fac*ngtm_stpl, color=(1,0,0,1), label='STPL', linestyle='--')
# VTPL
ngtm_vtpl = readerPlan.P_vtpl(mp1, means_vtpl)
fac *=10
plt.loglog(mp1, fac*ngtm, color=(0,1,0,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, nm*fac*ngtm_vtpl, color=(0,1,0,1), label='VTPL', linestyle='--')
# BCPL
ngtm_bcpl = readerPlan.P_bcpl(mp1, means_bcpl)
fac *=10
plt.loglog(mp1, fac*ngtm, color=(0,0,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, nm*fac*ngtm_bcpl, color=(0,0,1,1.0), label='BCPL', linestyle='--')
# TPL
ngtm_tpl = readerPlan.P_tpl(mp1, means_tpl)
fac *=10
plt.loglog(mp1, fac*ngtm, color=(1,0,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, nm*fac*ngtm_tpl, color=(1,0,1,1), label='TPL', linestyle='--')
# BPL
ngtm_bpl = readerPlan.P_bpl(mp1, means_bpl)
fac *=10
plt.loglog(mp1, fac*ngtm, color=(0,1,1,0.2), marker=".", linewidth=5, markersize=2)
plt.loglog(mp1, nm*fac*ngtm_bpl, color=(0,1,1,1), label='BPL', linestyle='--')
# plot labels etc.
plt.xlabel(r'$M_p$')
plt.ylabel(r'$N(>M_p)$')
plt.ylim(1.e-1, 1.e8)
plt.xlim(1e-5, 1.e-1)
plt.legend(prop={'size':6})
tools.saveAndClear(pathSave + "hist_cumulative.png", figNum=0)
